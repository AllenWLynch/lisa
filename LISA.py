
import build_chromatin_model
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

from utils import LoadingBar, Log

_config = configparser.ConfigParser()
_config.read('lisa_config.ini')


class LISA:

    factor_binding_technologies = dict(
        chipseq = 'ChIP-seq_1000kb',
        motifs = 'motif_hits_1000kb'
    )

    def __init__(self, species, *, 
        background_strategy = 'regulatory',
        num_background_genes = 3000,
        num_datasets_selected_anova = 200,
        num_datasets_selected = 10,
        threads = -1,
        isd_method = 'chipseq',
        use_covariates = False,
    ):
        self.isd_options = _config.get('lisa_params', 'isd_methods').split(',')
        self.background_options = _config.get('lisa_params', 'background_strategies').split(',')
        self.max_query_genes = int(_config.get('lisa_params', 'max_user_genelist_len'))

        assert( isinstance(num_background_genes, int) )
        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )
        assert( isinstance(threads, int) )
        self.threads = max(-1, threads)
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
        self.log = Log(stderr)
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

    def _load_gene_info(self, data_object):
        return (
            data_object[_config.get('gene_info', 'gene_symbols')][...].astype(str),
            data_object[_config.get('gene_info', 'gene_refseqIDs')][...].astype(str),
            data_object[_config.get('gene_info', 'tad_domains')][...].astype(str)
        )

    def _load_factor_binding_data(self, data_object):

        if self.isd_method == 'ChIP-seq_1000kb':
            dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
                .format(technology = 'ChIP-seq_RP')][...].astype(str))
        else:
            dataset_ids = list(data_object[self.isd_method].keys())

        loadingbar = LoadingBar('Reading {} data'.format(self.isd_method), len(dataset_ids), 20)

        peaks, dataset_index,factor_names = [],[],[]

        for i, dataset_id in enumerate(dataset_ids):
            self.log.append(loadingbar, update_line = True)

            dataset = data_object[_config.get('factor_binding', 'tf_binding_data').format(technology = self.isd_method, dataset_id = dataset_id)]
            peaks.append(dataset[...].astype(np.int64))
            dataset_index.append(np.full(dataset.shape, i))
            factor_names.append(dataset.attrs['factor'])

        self.log.append('Formatting {} data ...'.format(self.isd_method))
        peaks = np.concatenate(peaks)
        dataset_index = np.concatenate(dataset_index)

        factor_binding = sparse.csc_matrix(
            (
                np.ones(len(peaks)),
                (peaks, dataset_index)
            ),
            shape = (int(_config.get('lisa_params','num_bins')), len(dataset_ids))
        )

        return factor_binding, factor_names, dataset_ids

    def _load_rp_map(self, data_object):

        self.log.append('Importing regulatory potential map ...')

        rp_map_symbols = data_object[_config.get('RP_map', 'gene_symbols')][...].astype(str)
        rp_map = sparse.load_npz(_config.get('RP_map','matrix').format(species = self.species)).tocsc()

        return rp_map, rp_map_symbols

    def _load_data(self):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')
        try:
            with self.log.section('Loading data into memory (only on first prediction):') as log:

                with h5.File(self.data_source, 'r') as data:
                    
                    log.append('Loading gene information ...')
                    self.gene_symbols, self.gene_ids, self.tad_domains = self._load_gene_info(data)

                    with log.section('Loading binding data:') as log:
                        self.factor_binding, self.factor_names, self.factor_dataset_ids = self._load_factor_binding_data(data)

                    self.rp_map, self.rp_map_symbols = self._load_rp_map(data)

                    self.log.append('Done!')
            self._is_loaded = True
            
        except OSError:
            self._download_data()


    def _download_data(self):
        with log.section('Downloading {} data:'.format(self.species)):
            log.append('Fetching data from cistrome.org ...')
            log.append('Decompressing data ...')
            log.append('Done!\n')
            raise NotImplementedError()


    def load_rp_matrix(self, data_object, technology, label_dict):

        symbol_index = data_object[_config.get('accessibility_assay', 'reg_potential_gene_symbols').format(technology=technology)][...].astype(str)
        dataset_ids = data_object[_config.get('accessibility_assay', 'reg_potential_dataset_ids').format(technology = technology)][...].astype(str)
        symbol_mask = np.isin(symbol_index, list(label_dict.keys()))

        rp_matrix = data_object[_config.get('accessibility_assay', 'reg_potential_matrix').format(technology=technology)][symbol_mask, :][...]

        rp_matrix = rp_matrix[self.get_indices_in_order(symbol_index[symbol_mask], label_dict.keys()), :]

        return rp_matrix, dataset_ids

    def _train_chromatin_model(self, data_object, technology, label_dict):

        rp_matrix, dataset_ids = self.load_rp_matrix(data_object, technology, label_dict)

        sample_selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected)
        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-4,4,0.5))})

        narrow_rp_matrix, selected_dataset_ids, selection_model, chromatin_model, normalization_fn = build_chromatin_model\
            .build_chromatin_model(
                rp_matrix, dataset_ids, label_dict.values(), sample_selection_model, chromatin_model,
                num_anova_features = self.num_datasets_selected_anova, use_covariates = self.use_covariates, n_jobs = self.threads
            )

        return dict(
            selected_datasets = selected_dataset_ids,
            selection_model = selection_model,
            chromatin_model = chromatin_model,
            norm_fn = normalization_fn
        )


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

    @staticmethod
    def get_indices_in_order(from_list, target_ordering):
        from_map = {gene : i for i, gene in enumerate(from_list)}
        reordered = np.array([from_map[symbol] for symbol in target_ordering])
        return reordered

    def predict(self, query_list, background_list = []):

        self.log = Log(stderr)
        self.log.append('Initializing LISA ...')

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
                    label_dict = gene_selection.select_genes(query_list, self.gene_symbols, self.gene_ids, self.tad_domains, 
                        num_background_genes = self.num_background_genes, background_strategy = self.background_strategy, 
                        max_query_genes = self.max_query_genes, background_list = background_list)

                    num_query_genes = len([is_query for is_query in label_dict.values() if is_query])
                    self.log.append('Selected {} query genes and {} background genes.'\
                        .format(str(num_query_genes), str(len(label_dict) - num_query_genes)))

                    rp_map_symbols_idx = self.get_indices_in_order(self.rp_map_symbols, label_dict.keys())

                    subset_rp_map = self.rp_map[:, rp_map_symbols_idx].transpose()
               
                with self.log.section('Modeling DNase purturbations:'):
                    #DNase model building and purturbation
                    self.log.append('Selecting discriminative datasets and training chromatin model ...')

                    dnase_model = self._train_chromatin_model(data, 'DNase', label_dict)
                    
                    with self.log.section('Calculating in-silico deletions:'):

                        dnase_accessibility_profiles = self._load_accessibility_profiles(data, dnase_model['selected_datasets'], 'DNase')

                        dnase_test_stats, dnase_pvals = ISD.get_insilico_deletion_significance(log = self.log, accessibility_profiles = dnase_accessibility_profiles, 
                            factor_binding=self.factor_binding, rp_map = subset_rp_map, chromatin_model = dnase_model['chromatin_model'], labels = label_dict.values(), normalizer = dnase_model['norm_fn'])
                    
                    self.log.append('Done!')

                with self.log.section('Modeling H3K27ac purturbations:'):
                    self.log.append('Selecting discriminative datasets and training chromatin model ...')

                    acetylation_model = self._train_chromatin_model(data, 'H3K27ac', label_dict)
                    
                    with self.log.section('Calculating in-silico deletions:'):

                        acetylation_accessibility_profiles = self._load_accessibility_profiles(data, acetylation_model['selected_datasets'], 'H3K27ac')
                        
                        acetylation_test_stats, acetylation_pvals = ISD.get_insilico_deletion_significance(log = self.log, accessibility_profiles = acetylation_accessibility_profiles, 
                            factor_binding=self.factor_binding, rp_map = subset_rp_map, chromatin_model = acetylation_model['chromatin_model'], labels = label_dict.values(), normalizer = acetylation_model['norm_fn'])
                        
                    self.log.append('Done!')

                if self.isd_method == 'ChIP-seq_1000kb':

                    with self.log.section('Modeling direct ChIP-seq purturbations ...'):

                        chip_rp_matrix, dataset_ids = self.load_rp_matrix(data, 'ChIP-seq_RP', label_dict)

                        chip_test_stats, chip_pvals = ISD.get_delta_RP_p_value(chip_rp_matrix, label_dict.values())

                    aggregate_pvals = np.array([dnase_pvals, acetylation_pvals, chip_pvals]).T

                else:
                    aggregate_pvals = np.array([dnase_pvals, acetylation_pvals]).T

                with self.log.section('Mixing effects using Cauchy combination ...'):
                    
                    combined_test_stats, combined_p_values = ISD.cauchy_combination(aggregate_pvals)

                self.log.append('Formatting output ...')

                results = dict(
                    factor = self.factor_names,
                    dataset_id = self.factor_dataset_ids,
                    combined_p_value = combined_p_values,
                    combined_test_statistic = combined_test_stats,
                    Dnase_p_value = dnase_pvals,
                    Dnase_test_statistic = dnase_test_stats,
                    H3K27ac_p_value = acetylation_pvals,
                    H3K27ac_test_statistic = acetylation_test_stats,
                )

                if self.isd_method == 'ChIP-seq_1000kb':
                    results['Chip-seq_p_value'] = chip_pvals
                    results['Chip-seq_test_statistic'] = chip_test_stats

                results_table = list(zip(*results.values()))
                results_table = sorted(results_table, key = lambda col : col[2])

                results = dict( zip(results.keys(), list(zip(*results_table ))) )
                
                self.log.append('Done!')

                #return model_metadata as big dictionary
                return results, dict(
                    query_genes = [symbol for symbol, is_query in label_dict.items() if is_query],
                    background_genes = [symbol for symbol, is_query in label_dict.items() if not is_query],
                    Dnase_models= dict(
                        selection_model = dnase_model['selection_model'].get_info(),
                        chromatin_model = dnase_model['chromatin_model'].get_info(),
                        selected_datasets = list(dnase_model['selected_datasets'])
                    ),
                    H3K27ac_models = dict(
                        selection_model = acetylation_model['selection_model'].get_info(),
                        chromatin_model = acetylation_model['chromatin_model'].get_info(),
                        selected_datasets = list(acetylation_model['selected_datasets'])
                    ),
                )
                

        except OSError as err:
            self._download_data()
            print(err)

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
    parser.add_argument('-t','--threads', type = int, default = -1)

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