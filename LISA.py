
import gene_selection
import LISA_assays

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
    '''
    The LISA object is the user's interface with the LISA algorithm. It holds the parameters specified by the user and 
    handles data loading from hdf5 and 

    '''

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
        oneshot = False,
    ):
        self.isd_options = _config.get('lisa_params', 'isd_methods').split(',')
        self.background_options = _config.get('lisa_params', 'background_strategies').split(',')
        self.max_query_genes = int(_config.get('lisa_params', 'max_user_genelist_len'))

        assert( isinstance(num_background_genes, int) )
        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )

        self.num_background_genes = num_background_genes
        self.num_datasets_selected_anova = num_datasets_selected_anova
        self.num_datasets_selected = num_datasets_selected

        self._set_cores(cores)

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

        self.assays = []
        self.oneshot = oneshot

        self.log.append(
"""
___________________________________________________________________________________________________________________________

Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6
X. Shirley Liu Lab, 2020
___________________________________________________________________________________________________________________________
""")

    def _set_cores(self, cores):
        assert( isinstance(cores, int) )

        max_cores = multiprocessing.cpu_count() - 1
        if cores <= -1:
            cores = max_cores
        self.cores = min(cores, max_cores, self.num_datasets_selected)

    def _preprocess_gene_list(self, genes):
        return [str(gene).strip() for gene in genes if len(gene) > 0]

    def _load_gene_info(self):

        all_genes = GeneSet()
        
        with open(_config.get('genes','master_gene_list').format(species = self.species), 'r') as genes:
            all_genes.from_str(genes.read())

        return all_genes
        

    def _load_factor_binding_data(self, data_object):

        self.factor_dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
            .format(technology = self.isd_method)][...].astype(str))
        
        return sparse.load_npz(_config.get('factor_binding','matrix').format(species = self.species, technology = self.isd_method))

    def _load_rp_map(self, data_object):

        rp_map = sparse.load_npz(_config.get('RP_map','matrix').format(species = self.species)).tocsr()

        #select genes of interest from rp_map
        with open(_config.get('genes','gene_locs').format(species = self.species), 'r') as f:
            rp_map_locs = np.array([line.strip() for line in f.readlines()])

        return rp_map, rp_map_locs


    def add_assay(self, data_object, assay):
        if not self.oneshot:
            assay.load(data_object)
        self.assays.append(assay)


    def _load_data(self):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')
        try:
            with self.log.section('Loading data into memory (only on first prediction):') as log:

                with h5.File(self.data_source, 'r') as data:
                    
                    log.append('Loading gene information ...')
                    self.all_genes = self._load_gene_info()

                    with log.section('Loading binding data ...') as log:
                        self.factor_binding = self._load_factor_binding_data(data)

                    log.append('Loading regulatory potential map ...')
                    self.rp_map, self.rp_map_locs = self._load_rp_map(data)

                    self.add_assay(data, 
                        LISA_assays.PeakRP_Assay(self.isd_method, _config, self.cores, self.log)
                    )

                    self.add_assay(data, 
                        LISA_assays.Accesibility_Assay('DNase', _config, self.cores, self.log,
                            rp_map = self.rp_map, factor_binding = self.factor_binding,
                            selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                            chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                        )
                    )

                    self.add_assay(data, 
                        LISA_assays.Accesibility_Assay('H3K27ac', _config, self.cores, self.log,
                            rp_map = self.rp_map, factor_binding = self.factor_binding,
                            selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                            chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                        )
                    )

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

    def _get_metadata(self, sample_ids):
        return dict(
            factor_name = [self.metadict[_id]['factor'] for _id in sample_ids],
            sample_id = sample_ids,
            cell_line = [self.metadict[_id]['cell_line'] for _id in sample_ids],
            cell_type = [self.metadict[_id]['cell_type'] for _id in sample_ids],
            tissue = [self.metadict[_id]['tissue'] for _id in sample_ids]
        )

    
    @staticmethod
    def combine_tests(p_vals, weights=None):
        #https://arxiv.org/abs/1808.09011
        p_vals = np.array(p_vals).astype(np.float64)

        assert(len(p_vals.shape) == 2 ), 'P-values must be provided as matrix of (samples, multiple p-values)'
        assert(p_vals.shape[1] > 1), 'Must have multiple p-values to combine'

        #clip p-values at minimum float64 value to prevent underflow
        p_vals = np.clip(p_vals, np.nextafter(0, 1, dtype = np.float64), 1.0)
        if weights is None:
            weights = np.ones((1, p_vals.shape[1])) / p_vals.shape[1]
        else:
            assert(p_vals.shape[1] == weights.shape[1])

        test_statistic = np.sum(weights * np.tan((0.5-p_vals) * np.pi), axis = 1)
        combined_p_value = 0.5 - np.arctan(test_statistic)/np.pi

        return combined_p_value

    #process user-supplied text
    def select_genes(self, query_list, background_list):
        '''
        query_list: list of any text that may be a query gene
        background_list : list of any text that may be a background gene
        These lists are lightly processed, then matched with possible genes in the lisa database

        returns:
        gene_mask: mask for which genes in lisa were selected for query and background
        label_vector: query or background set membership for each of those genes
        gene_info_dict: gene info for results printout
        '''
        query_list = self._preprocess_gene_list(query_list)
        background_list = self._preprocess_gene_list(background_list)
        
        #match user-supplied text with genes, and select background genes
        label_dict, gene_info_dict = gene_selection.select_genes(query_list, self.all_genes, 
            num_background_genes = self.num_background_genes, background_strategy = self.background_strategy, 
            max_query_genes = self.max_query_genes, background_list = background_list)

        #Show user number of query and background genes selected
        self.log.append('Selected {} query genes and {} background genes.'\
            .format(str(len(gene_info_dict['query_symbols'])), str(len(gene_info_dict['background_symbols']))))

        #subset the rp map for those genes and make label vector
        gene_mask = np.isin(self.rp_map_locs, list(label_dict.keys()))
        #make label vector from label_dict based of gene_loc ordering in rp_map data
        label_vector = np.array([label_dict[gene_loc] for gene_loc in self.rp_map_locs[gene_mask]])
        
        return gene_mask, label_vector, gene_info_dict

    def predict(self, query_list, background_list = []):

        self.log.append('Initializing LISA using {} cores ...'.format(str(self.cores)))

        if not self._is_loaded:
            self._load_data()

        assert( isinstance(query_list, Iterable) )

        if background_list is None:
            background_list = []
        
        assert( isinstance(background_list, Iterable))

        try:
                            
            with self.log.section('Matching genes and selecting background ...'):

                if len(background_list) > 0 and self.background_strategy == 'provided':
                    self.log.append('User provided background genes!')
                
            gene_mask, label_vector, gene_info_dict = self.select_genes(query_list, background_list)

            with h5.File(self.data_source, 'r') as data:

                assay_pvals, assay_info = {},{}
                for assay in self.assays:
                    assay_pvals[assay.technology + '_p_value'] = assay.predict(gene_mask, label_vector, data_object = data)
                    assay_info[assay.technology + '_model_info'] = assay.get_info(self._get_metadata)

            with self.log.section('Mixing effects using Cauchy combination ...'):

                aggregate_pvals = np.array(list(assay_pvals.values())).T
                
                combined_p_values = self.combine_tests(aggregate_pvals)

        except OSError as err:
            self._download_data()
            print(err)

        self.log.append('Formatting output ...')

        results = dict(
            **self._get_metadata(self.factor_dataset_ids),
            combined_p_value = combined_p_values,
            combined_p_value_adjusted = list(np.minimum(np.array(combined_p_values) * len(combined_p_values), 1.0)),
            **assay_pvals,
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
            **gene_info_dict,
            results = results,
            **assay_info
        )

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