
import gene_selection
import LISA_assays

from models import LR_BinarySearch_SampleSelectionModel
from models import LR_ChromatinModel

import configparser
import argparse

from collections.abc import Iterable
import sys
import numpy as np
import h5py as h5
from scipy import sparse
import os
import json
import multiprocessing
from utils import LoadingBar, Log, LISA_Results

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

    def __init__(self, species, 
        cores = 1,
        knockout_method = 'chipseq',
        num_datasets_selected_anova = 200,
        num_datasets_selected = 10,
        verbose = True,
        oneshot = False,
        log = None,
    ):
        self.knockout_options = _config.get('lisa_params', 'isd_methods').split(',')
        self.background_options = _config.get('lisa_params', 'background_strategies').split(',')
        self.max_query_genes = int(_config.get('lisa_params', 'max_user_genelist_len'))

        
        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )

        self.num_datasets_selected_anova = num_datasets_selected_anova
        self.num_datasets_selected = num_datasets_selected

        self._set_cores(cores)

        assert( num_datasets_selected_anova > num_datasets_selected ), 'Anova must select more datasets than the regression model'
        assert( num_datasets_selected_anova > 0 and num_datasets_selected > 0 ), 'Number of datasets selected must be positive'

        assert( knockout_method in  self.knockout_options ), 'ISD method must be \{{}\}'.format(', '.join(self.knockout_options))
        assert( species in ['mm10','hg38'] ), "Species must be either hg38 or mm10"
        
        self.isd_method = self.factor_binding_technologies[knockout_method]
        self.species = species

        self.data_source = _config.get('paths', 'h5_path').format(species = self.species)

        if log is None:
            self.log = Log(sys.stderr, verbose = verbose)
        else:
            self.log = log

        self._is_loaded = False

        self.assays = []
        self.oneshot = oneshot

        #load all gene information
        self.all_genes = self._load_gene_info()

        #load gene axis labels for all assays, both these loading steps are very fast, and allow for more finely-tuned loading from the hdf5 later
        with open(_config.get('genes','gene_locs').format(species = self.species), 'r') as f:
            self.rp_map_locs = np.array([line.strip() for line in f.readlines()])

        self.used_oneshot = False


    def _set_cores(self, cores):
        assert( isinstance(cores, int) )

        max_cores = multiprocessing.cpu_count() - 1
        if cores <= -1:
            cores = max_cores
        self.cores = min(cores, max_cores, self.num_datasets_selected)

    def _preprocess_gene_list(self, genes):
        return [str(gene).strip() for gene in genes if len(gene) > 0]

    def _load_gene_info(self):

        all_genes = gene_selection.GeneSet()
        
        with open(_config.get('genes','master_gene_list').format(species = self.species), 'r') as genes:
            all_genes.from_str(genes.read())

        return all_genes
        

    def _load_factor_binding_data(self, data_object):

        self.factor_dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
            .format(technology = self.isd_method)][...].astype(str))
        
        return sparse.load_npz(_config.get('factor_binding','matrix').format(species = self.species, technology = self.isd_method))

    def _load_rp_map(self, data_object):

        rp_map = sparse.load_npz(_config.get('RP_map','matrix').format(species = self.species)).tocsr()

        return rp_map


    def add_assay(self, data_object, gene_mask, assay):
        assay.load(data_object, gene_mask=gene_mask, oneshot=self.oneshot)
        self.assays.append(assay)


    def _load_data(self, gene_mask = None):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')
        try:
            with self.log.section('Loading data into memory (only on first prediction):') as log:

                with h5.File(self.data_source, 'r') as data:

                    with log.section('Loading binding data ...') as log:
                        self.factor_binding = self._load_factor_binding_data(data)

                    log.append('Loading regulatory potential map ...')
                    self.rp_map = self._load_rp_map(data)

                    self.add_assay(data, gene_mask,
                        LISA_assays.PeakRP_Assay(self.isd_method, _config, self.cores, self.log)
                    )

                    self.add_assay(data, gene_mask,
                        LISA_assays.Accesibility_Assay('DNase', _config, self.cores, self.log,
                            rp_map = self.rp_map, factor_binding = self.factor_binding,
                            selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                            chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                        )
                    )

                    self.add_assay(data, gene_mask,
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
    def make_gene_mask(self, query_list, background_list, num_background_genes, background_strategy):
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
            num_background_genes = num_background_genes, background_strategy = background_strategy, 
            max_query_genes = self.max_query_genes, background_list = background_list)

        #Show user number of query and background genes selected
        self.log.append('Selected {} query genes and {} background genes.'\
            .format(str(len(gene_info_dict['query_symbols'])), str(len(gene_info_dict['background_symbols']))))

        #subset the rp map for those genes and make label vector
        gene_mask = np.isin(self.rp_map_locs, list(label_dict.keys()))
        #make label vector from label_dict based of gene_loc ordering in rp_map data
        label_vector = np.array([label_dict[gene_loc] for gene_loc in self.rp_map_locs[gene_mask]])
        
        return gene_mask, label_vector, gene_info_dict

    def predict(self, query_list, 
            background_list = [], 
            background_strategy = 'regulatory',
            num_background_genes = 3000,
        ):

        if self.oneshot and self.used_oneshot:
            raise AssertionError('When instantiated in one-shot, model cannot be used for multiple predictions')

        if background_list is None:
            background_list = []
        assert( isinstance(num_background_genes, int) )
        assert( background_strategy in  self.background_options ), 'Background strategy must be in \{{}\}'.format(', '.join(self.background_options))
        assert( isinstance(query_list, Iterable) )
        assert( isinstance(background_list, Iterable))

        self.log.append('Initializing LISA using {} cores ...'.format(str(self.cores)))

        try:
                            
            with self.log.section('Matching genes and selecting background ...'):

                if len(background_list) > 0 and self.background_strategy == 'provided':
                    self.log.append('User provided background genes!')
                
            gene_mask, label_vector, gene_info_dict = self.make_gene_mask(query_list, background_list, num_background_genes, background_strategy)

            if not self._is_loaded:
                self._load_data(gene_mask)

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

        results = LISA_Results.fromdict(
            **self._get_metadata(self.factor_dataset_ids),
            combined_p_value = combined_p_values,
            combined_p_value_adjusted = list(np.minimum(np.array(combined_p_values) * len(combined_p_values), 1.0)),
            **assay_pvals
        )

        results = results.sortby('combined_p_value_adjusted', add_rank = True)
        
        self.log.append('Done!')
        
        self.used_oneshot = True
        #return model_metadata as big dictionary
        return results, dict(
            **gene_info_dict,
            **assay_info
        )


INSTANTIATION_KWARGS = ['cores','knockout_method','verbose','oneshot']
PREDICTION_KWARGS = ['background_list','num_background_genes','background_strategy']

def extract_kwargs(args, keywords):
    return {key : vars(args)[key] for key in keywords}

def is_valid_prefix(prefix):

    if os.path.isdir(prefix) or os.path.isfile(prefix) or os.path.isdir(os.path.dirname(prefix)):
        return prefix
        
    raise argparse.ArgumentTypeError('{}: Invalid file prefix.'.format(prefix))


def lisa_oneshot(args):

    results, metadata = LISA(args.species, **extract_kwargs(args, INSTANTIATION_KWARGS)).predict(args.query_list.readlines(), **extract_kwargs(args, PREDICTION_KWARGS))
    
    if args.save_metadata:
        if args.output_prefix:
            metadata_filename = args.output_prefix + '.metadata.json' 
        else:
            metadata_filename = os.path.basename(args.query_list.name) + '.metadata.json'

        with open(metadata_filename, 'w') as f:
            f.write(json.dumps(metadata, indent=4))

    if not args.output_prefix is None:
        with open(args.output_prefix + '.lisa.tsv', 'w') as f:
            f.write(results.to_tsv())
    else:
        print(results.to_tsv())
        
def lisa_multi(args):

    log = Log(target = sys.stderr, verbose = args.verbose)

    lisa = LISA(args.species, **extract_kwargs(args, INSTANTIATION_KWARGS), log = log)

    results_summary = []

    for query in args.query_lists:
        
        query_name = os.path.basename(query.name)

        with log.section('Modeling gene list {}:'.format(str(query_name))):

            results, metadata = lisa.predict(query.readlines(), **extract_kwargs(args, PREDICTION_KWARGS))

            with open(args.output_prefix + query_name + '.lisa.tsv', 'w') as f:
                f.write(results.to_tsv())

            if args.save_metadata:
                with open(args.output_prefix + query_name + '.metadata.json', 'w') as f:
                    f.write(json.dumps(metadata, indent=4))

            top_TFs = results.subset(range(0,20)).todict()['factor_name']
            
            top_TFs_unique, encountered = [], set()
            for TF in top_TFs:
                if not TF in encountered:
                    top_TFs_unique.append(TF)
                    encountered.add(TF)
        

            results_summary.append((query_name, top_TFs_unique))

    print('Sample\tTop Regulatory Factors')
    for result_line in results_summary:
        print(result_line[0], ', '.join(result_line[1]), sep = '\t')


def lisa_one_v_rest():

    raise NotImplementedError()


def build_common_args(parser):
    parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
    parser.add_argument('-c','--cores', required = True, type = int, default = 1)
    parser.add_argument('--seed', type = int, default = None, help = 'Random seed for gene selection. Allows for reproducing exact results.')
    parser.add_argument('--knockout_method', type = str, choices = _config.get('lisa_params','isd_methods').split(','), default = 'chipseq',
        help = 'Use ChIP-seq peaks (recommended) or motif hits to represent TF binding')
    parser.add_argument('--save_metadata', action = 'store_true', default = False, help = 'Save json-formatted metadata from processing each gene list.')


def build_multiple_lists_args(parser):
    parser.add_argument('query_lists', type = argparse.FileType('r', encoding = 'utf-8'), nargs = "+", help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
    parser.add_argument('-o','--output_prefix', required = True, type = is_valid_prefix, help = 'Output file prefix.')
    parser.add_argument('-b','--num_background_genes', type = int, default = _config.get('lisa_params', 'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes. These genes are selection from other gene lists.')
    parser.add_argument('--random_background', action = 'store_const', const = 'random', default = 'regulatory', dest = 'background_strategy', help = 'Use random background selection rather than "regulatory" selection.')
    parser.add_argument('-v','--verbose',type = int, default = 2)

if __name__ == "__main__":

    #define command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description =
"""
___________________________________________________________________________________________________________________________

Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6
X. Shirley Liu Lab, 2020
___________________________________________________________________________________________________________________________
""")

    subparsers = parser.add_subparsers(help = 'commands')

    #__ LISA oneshot command __#################

    oneshot_parser = subparsers.add_parser('oneshot', help = 'Use LISA to infer genes from one gene list. If you have multiple lists, this option will be slower than using "multi" due to data-loading time.\n')
    build_common_args(oneshot_parser)

    oneshot_parser.add_argument('query_list', type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
    oneshot_parser.add_argument('-o','--output_prefix', required = False, type = is_valid_prefix, help = 'Output file prefix. If left empty, will write results to stdout.')
    oneshot_parser.add_argument('--background_strategy', choices = _config.get('lisa_params', 'background_strategies').split(','),
        default = 'regulatory',
        help = """Background genes selection strategy. LISA samples background genes to compare to user\'s genes-of-interest from a diverse
        regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
        """)
    background_genes_group = oneshot_parser.add_mutually_exclusive_group()
    background_genes_group.add_argument('--background_list', type = argparse.FileType('r', encoding = 'utf-8'), required = False,
        help = 'user-supplied list of backgroung genes. Used when --background_strategy flag is set to "provided"')
    background_genes_group.add_argument('-b','--num_background_genes', type = int, default = _config.get('lisa_params', 'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes')
    oneshot_parser.add_argument('-v','--verbose',type = int, default = 4)
    oneshot_parser.set_defaults(func = lisa_oneshot, oneshot = True)
    
    #__ LISA multi command __#################

    multi_parser = subparsers.add_parser('multi', help = 'Process multiple genelists, provided one list per line, genes seperated by commas. This reduces data-loading time if using the same parameters for all lists.\n')
    build_common_args(multi_parser)
    build_multiple_lists_args(multi_parser)
    multi_parser.set_defaults(func = lisa_multi, oneshot = False, background_list = None)
    
    #__ LISA one-vs-rest command __#################

    one_v_rest_parser = subparsers.add_parser('one-vs-rest', help = 'Compare gene lists in a one-vs-rest fashion. Useful downstream of cluster analysis.\n')
    build_common_args(one_v_rest_parser)
    build_multiple_lists_args(one_v_rest_parser)
    one_v_rest_parser.set_defaults(func = lisa_one_v_rest, oneshot = False, background_list = None)

    args = parser.parse_args()

    try:
        args.func(args)
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
