import gene_selection
import insilico_deletion_model as ISD

from models import LR_BinarySearch_SampleSelectionModel
from models import LR_ChromatinModel

import configparser
import argparse

from collections.abc import Iterable
from sys import stderr
import numpy as np
import h5py as h5
from scipy import stats, sparse
import os
import itertools
import json
import multiprocessing

from utils import LoadingBar, Log, GeneSet, Gene

_config = configparser.ConfigParser()
_config.read('lisa_config.ini')


class LISA:

    factor_binding_technologies = dict(
        chipseq = 'ChIP-seq',
        motifs = 'Motifs'
    )

    def __init__(self, species, *, 
        background_strategy = 'regulatory',
        num_background_genes = 3000,
        num_datasets_selected_anova = 200,
        num_datasets_selected = 10,
        cores = -1,
        isd_method = 'chipseq',
        use_covariates = False,
        verbose = True,
    ):
        self.isd_options = _config.get('lisa_params', 'isd_methods').split(',')
        self.background_options = _config.get('lisa_params', 'background_strategies').split(',')
        self.max_query_genes = int(_config.get('lisa_params', 'max_user_genelist_len'))

        assert( isinstance(num_background_genes, int) )
        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )
        assert( isinstance(cores, int) )
        if cores <= -1:
            cores = multiprocessing.cpu_count()
        self.cores = min(cores, multiprocessing.cpu_count(), num_datasets_selected)
        self.num_background_genes = num_background_genes
        self.num_datasets_selected_anova = num_datasets_selected_anova
        self.num_datasets_selected = num_datasets_selected

        assert( num_datasets_selected_anova > num_datasets_selected ), 'Anova must select more datasets than the regression model'
        assert( num_datasets_selected_anova > 0 and num_datasets_selected > 0 ), 'Number of datasets selected must be positive'

        assert( background_strategy in  self.background_options ), 'Background strategy must be in \{{}\}'.format(', '.join(self.background_options))
        assert( isd_method in  self.isd_options ), 'ISD method must be \{{}\}'.format(', '.join(self.isd_options))
        assert( species in ['mm10','hg38'] ), "Species must be either hg38 or mm10"
        
        self.background_strategy = background_strategy
        self.isd_method = self.factor_binding_technologies[isd_method]
        self.species = species

        self.use_covariates = use_covariates
        self.data_source = _config.get('paths', 'h5_path').format(species = self.species)
        self.log = Log(stderr if verbose else os.devnull)
        self._is_loaded = False

        self.log.append(
"""
___________________________________________________________________________________________________________________________

Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6
X. Shirley Liu Lab, 2020
___________________________________________________________________________________________________________________________
""")


    def _preprocess_gene_list(self, genes):
        return [str(gene).strip() for gene in genes if len(gene) > 0]

    def _load_gene_info(self):

        self.all_genes = GeneSet()
        
        with open(_config.get('genes','master_gene_list').format(species = self.species), 'r') as genes:
            self.all_genes.from_str(genes.read())
        

    def _load_factor_binding_data(self, data_object):

        dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
            .format(technology = self.isd_method)][...].astype(str))
        
        binding_sparse = sparse.load_npz(_config.get('factor_binding','matrix').format(species = self.species, technology = self.isd_method))

        self.factor_binding = binding_sparse
        self.factor_dataset_ids = dataset_ids


    def _load_rp_map(self, data_object):

        self.rp_map = sparse.load_npz(_config.get('RP_map','matrix').format(species = self.species)).tocsr()

        #select genes of interest from rp_map
        with open(_config.get('genes','gene_locs').format(species = self.species), 'r') as f:
            self.rp_map_locs = np.array([line.strip() for line in f.readlines()])

    def _load_data(self):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')
        try:
            with self.log.section('Loading data into memory (only on first prediction):') as log:

                with h5.File(self.data_source, 'r') as data:
                    
                    log.append('Loading gene information ...')
                    self._load_gene_info()

                    with log.section('Loading binding data ...') as log:
                        self._load_factor_binding_data(data)

                    log.append('Loading regulatory potential map ...')
                    self._load_rp_map(data)

                    log.append('Loading direct binding peak-RP model ...')
                    self.direct_rp_matrix, _ = self.load_rp_matrix(data, self.isd_method, np.ones(len(self.rp_map_locs)).astype(np.bool))

                self.log.append('Loading metadata ...')

                with open(_config.get('paths','metadata'), 'r', encoding = 'latin') as metdata_file:
                    metadata = [[field.strip() for field in line.split('\t')] for line in metdata_file.readlines()]
                                        
                meta_headers, metadata = metadata[0], metadata[1:]

                self.metadict = {metaline[0] : dict(zip(meta_headers[1:], metaline[1:])) for metaline in metadata}

                self.log.append('Done!')
            self._is_loaded = True
            
        except OSError:
            self._download_data()

    def _download_data(self):
        with self.log.section('Downloading {} data:'.format(self.species)):
            self.log.append('Fetching data from cistrome.org ...')
            self.log.append('Decompressing data ...')
            self.log.append('Done!\n')
            raise NotImplementedError()


    def load_rp_matrix(self, data_object, technology, location_mask):

        #symbol_index = data_object[_config.get('accessibility_assay', 'reg_potential_gene_symbols').format(technology=technology)][...].astype(str)
        dataset_ids = data_object[_config.get('accessibility_assay', 'reg_potential_dataset_ids').format(technology = technology)][...].astype(str)
        #symbol_mask = np.isin(symbol_index, list(label_dict.keys()))

        rp_matrix = data_object[_config.get('accessibility_assay', 'reg_potential_matrix').format(technology=technology)][location_mask, :][...]

        #subset RP matrix here

        return rp_matrix, dataset_ids

    def _load_accessibility_profiles(self, data_object, selected_datasets, accessibility_assay):
        #dataset_ids = list(zip(*[key.split(',') for key in data_object[factor_binding_technology].keys()]))
        loadingbar = LoadingBar('Reading {} data'.format(accessibility_assay), len(selected_datasets), 20)

        accessibility_profiles = []
        for selected_dataset in selected_datasets:
            self.log.append(loadingbar, update_line = True)
            accessibility_profiles.append(
                data_object[_config.get('accessibility_assay','binned_reads')\
                .format(technology = accessibility_assay, dataset_id = selected_dataset)][...][:,np.newaxis]
            )
        accessibility_profiles = np.concatenate(accessibility_profiles, axis = 1)

        return accessibility_profiles


    def _train_chromatin_model(self, data_object, technology, location_mask, label_vector):

        rp_matrix, dataset_ids = self.load_rp_matrix(data_object, technology, location_mask)

        selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected)
        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')

        dataset_mask = selection_model.fit(rp_matrix, label_vector )

        subset_rp_matrix, selected_dataset_ids = rp_matrix[:, dataset_mask], dataset_ids[dataset_mask]

        chromatin_model.fit(subset_rp_matrix, label_vector )

        return dict(selected_datasets = selected_dataset_ids, selection_model = selection_model, chromatin_model = chromatin_model)


    def _get_metadata(self, sample_ids):
        return dict(
            factor_name = [self.metadict[_id]['factor'] for _id in sample_ids],
            sample_id = sample_ids,
            cell_line = [self.metadict[_id]['cell_line'] for _id in sample_ids],
            cell_type = [self.metadict[_id]['cell_type'] for _id in sample_ids],
            tissue = [self.metadict[_id]['tissue'] for _id in sample_ids]
        )

    def predict(self, query_list, background_list = []):

        self.log = Log(stderr)
        self.log.append('Initializing LISA ...')
        self.log.append('Using {} cores.'.format(str(self.cores)))

        if not self._is_loaded:
            self._load_data()

        assert( isinstance(query_list, Iterable) )

        if background_list is None:
            background_list = []
        
        assert( isinstance(background_list, Iterable))

        try:
            with h5.File(self.data_source, 'r') as data:
                
                with self.log.section('Matching genes and selecting background ...'):

                    if len(background_list) > 0 and self.background_strategy == 'provided':
                        self.log.append('User provided background genes!')
                    #process user-supplied text
                    query_list = self._preprocess_gene_list(query_list)
                    background_list = self._preprocess_gene_list(background_list)

                    #match user-supplied text with genes, and select background genes
                    label_dict, query_symbols, background_symbols = gene_selection.select_genes(query_list, self.all_genes, 
                        num_background_genes = self.num_background_genes, background_strategy = self.background_strategy, 
                        max_query_genes = self.max_query_genes, background_list = background_list)

                    #Show user number of query and background genes selected
                    self.log.append('Selected {} query genes and {} background genes.'\
                        .format(str(len(query_symbols)), str(len(background_symbols))))

                    #subset the rp map for those genes and make label vector
                    rp_map_mask = np.isin(self.rp_map_locs, list(label_dict.keys()))
                    gene_subset_rp_map = self.rp_map[rp_map_mask, : ]

                    rp_map_locs_subset = self.rp_map_locs[rp_map_mask]
                    label_vector = np.array([label_dict[gene_loc] for gene_loc in rp_map_locs_subset])

                    #find bins of gene-subsetted rp-map with RP > 0
                    bin_subset = np.squeeze(np.array(gene_subset_rp_map.tocsc().sum(axis = 0) > 0))

                    #subset rp_map and factor hits on bins with RP > 0
                    subset_factor_binding = self.factor_binding[bin_subset, :]
                    bin_gene_subset_rp_map = gene_subset_rp_map[:, bin_subset]
               
                with self.log.section('Modeling DNase purturbations:'):
                    #DNase model building and purturbation
                    self.log.append('Selecting discriminative datasets and training chromatin model ...')

                    dnase_model = self._train_chromatin_model(data, 'DNase', rp_map_mask, label_vector)
                    
                    with self.log.section('Calculating in-silico deletions:'):

                        dnase_accessibility_profiles = self._load_accessibility_profiles(data, dnase_model['selected_datasets'], 'DNase')[bin_subset, :]
                        
                        dnase_test_stats, dnase_pvals, dnase_deltaR_scores = ISD.get_insilico_deletion_significance(log = self.log, accessibility_profiles = dnase_accessibility_profiles, 
                            factor_binding=subset_factor_binding, rp_map = bin_gene_subset_rp_map, chromatin_model = dnase_model['chromatin_model'], labels = label_vector, cores = self.cores)
                    
                    self.log.append('Done!')

                with self.log.section('Modeling H3K27ac purturbations:'):
                    self.log.append('Selecting discriminative datasets and training chromatin model ...')

                    acetylation_model = self._train_chromatin_model(data, 'H3K27ac', rp_map_mask, label_vector)
                    
                    with self.log.section('Calculating in-silico deletions:'):

                        acetylation_accessibility_profiles = self._load_accessibility_profiles(data, acetylation_model['selected_datasets'], 'H3K27ac')[bin_subset, :]
                        
                        acetylation_test_stats, acetylation_pvals, acetyaltion_deltaR_scores = ISD.get_insilico_deletion_significance(log = self.log, accessibility_profiles = acetylation_accessibility_profiles, 
                            factor_binding=subset_factor_binding, rp_map = bin_gene_subset_rp_map, chromatin_model = acetylation_model['chromatin_model'], labels = label_vector, cores = self.cores)
                        
                    self.log.append('Done!')

                with self.log.section('Modeling direct RP purturbations ...'):

                    self.log.append('Calculating direct binding p-values ...')
                    direct_rp_matrix = np.log2( self.direct_rp_matrix[rp_map_mask, :] + 1 )
                    
                    direct_test_stats, direct_p_vals = ISD.get_delta_RP_p_value(direct_rp_matrix, label_vector, cores = self.cores)

                with self.log.section('Mixing effects using Cauchy combination ...'):

                    aggregate_pvals = np.array([dnase_pvals, acetylation_pvals, direct_p_vals]).T
                    
                    combined_test_stats, combined_p_values = ISD.cauchy_combination(aggregate_pvals)

        except OSError as err:
            self._download_data()
            print(err)

        self.log.append('Formatting output ...')

        results = dict(
            **self._get_metadata(self.factor_dataset_ids),
            combined_p_value = combined_p_values,
            combined_p_value_adjusted = list(np.array(combined_p_values) * len(combined_p_values)),
            combined_test_statistic = combined_test_stats,
            DNase_p_value = dnase_pvals,
            DNase_test_statistic = dnase_test_stats,
            H3K27ac_p_value = acetylation_pvals,
            H3K27ac_test_statistic = acetylation_test_stats,
            Direct_knockout_test_statistic = direct_test_stats,
            Direct_knockout_p_value = direct_p_vals,
        )

        #convert to row-indexed data and sort by combined p-value
        results_table = list(zip(*results.values()))
        results_table = sorted(results_table, key = lambda col : col[5])
        
        #eliminate duplicate TF samples for summary results
        encountered = {}
        summary = []
        for result_line in results_table:
            factor_name = result_line[0]
            if not factor_name in encountered:
                summary.append(result_line[:7])
                encountered[factor_name] = True

        #transpose to column-indexed data and add "Rank" column
        summary_columns = [list(range(1, len(summary) + 1)), * list(zip(*summary))]
        summary_headers = ['Rank', * results.keys()]
        #convert to dictionary
        summary = dict( zip(summary_headers, summary_columns) )
        
        #convert to dictionary
        results = dict( zip(results.keys(), list(zip(*results_table))) )
       
        self.log.append('Done!')

        #return model_metadata as big dictionary
        return summary, dict(
            query_genes = query_symbols,
            background_genes = background_symbols,
            query_gene_locs = [gene_loc for gene_loc, is_query in label_dict.items() if is_query],
            background_gene_locs = [gene_loc for gene_loc, is_query in label_dict.items() if not is_query],
            all_results = results,
            DNase_models= dict(
                selection_model = dnase_model['selection_model'].get_info(),
                chromatin_model = dnase_model['chromatin_model'].get_info(),
                selected_datasets = self._get_metadata(list(dnase_model['selected_datasets'])),
            ),
            H3K27ac_models = dict(
                selection_model = acetylation_model['selection_model'].get_info(),
                chromatin_model = acetylation_model['chromatin_model'].get_info(),
                selected_datasets = self._get_metadata(list(acetylation_model['selected_datasets'])),
            ),
        ), dnase_deltaR_scores, dnase_model['chromatin_model'].rp_0

    @staticmethod
    def pretty_print_results(results_dict, top_n = 200):
        
        headers = list(results_dict.keys())
        rows = list(zip(*[results_list[:top_n] for results_list in results_dict.values()]))
        for line in [headers, *rows]:
            print('\t'.join([str(value) for value in line]))


if __name__ == "__main__":

    #define command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description = 'Landscape In Silico deletion Analysis (LISA): Find transcription factors important to the regulation of your gene set using using ChIP-seq informed purturbations of chromatin landscape models.'
        )

    parser.add_argument('query_genes', type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene list. One gene per line in either symbol or refseqID format')

    parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')

    parser.add_argument('-o','--metadata_output', required = False, type = argparse.FileType('w', encoding = 'utf-8'), help = 'metadata file output prefix. If left empty, will not save.')

    parser.add_argument('-n','--num_results', type = int, default = 200, help = 'Number of predictions to return.')

    parser.add_argument('-s', '--background_strategy', choices = _config.get('lisa_params', 'background_strategies').split(','),
        default = 'regulatory',
        help = """Background genes selection strategy. LISA samples background genes to compare to user\'s genes-of-interest from a diverse
        regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
        """)

    parser.add_argument('--background_genes', type = argparse.FileType('r', encoding = 'utf-8'), required = False,
        help = 'user-supplied list of backgroung genes. Used when -s flag is set to "provided"')

    parser.add_argument('-b','--num_background_genes', type = int, default = _config.get('lisa_params', 
    'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes')

    parser.add_argument('-a','--num_datasets_selected_anova', type = int, default = _config.get('lisa_params', 'num_anova_selected') or 200,
        help = 'Number of datasets to select using anova as preliminary filter')

    parser.add_argument('-d','--num_datasets_selected', type = int, default = _config.get('lisa_params','num_datasets_selected') or 10,
        help = 'Number of discriminative datasets to select using logistic regression model. Number of datasets selected linearly impacts performance')
    parser.add_argument('-c','--cores', type = int, default = -1)

    parser.add_argument('-m','--isd_method', type = str, choices = _config.get('lisa_params','isd_methods').split(','), default = 'chipseq',
        help = 'Use ChIP-seq peaks (recommended) or motif hits to represent TF binding')

    args = parser.parse_args()

    query_list = args.query_genes.readlines()
    del args.query_genes

    background_list = args.background_genes
    del args.background_genes

    output_filestream = args.metadata_output
    del args.metadata_output

    num_results = args.num_results
    del args.num_results
    
    lisa = LISA(**vars(args))
    predictions, metadata = lisa.predict(query_list, background_list)

    if not output_filestream is None:
        metadata['args'] = dict(**vars(args))
        output_filestream.write(json.dumps(metadata, indent = 4))
        output_filestream.close()

    lisa.pretty_print_results(predictions, top_n = num_results)