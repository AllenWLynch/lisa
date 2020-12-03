
#lisa modules
import lisa.gene_selection as gene_selection
from lisa.utils import LoadingBar, Log, LISA_Results, Metadata
from lisa._version import __version__

#standard library
import configparser
from collections.abc import Iterable
import sys
import json
import multiprocessing
from urllib import request, error
import tarfile
import os
import sys

#required
import numpy as np
import h5py as h5
from scipy import sparse, stats
#from pyranges.pyranges import PyRange

PACKAGE_PATH = os.path.dirname(__file__)
CONFIG_PATH = os.path.join(PACKAGE_PATH, 'config.ini')
REQURED_DATASET_VERSION = '.'.join(__version__.split('.')[:2])
INSTALL_PATH = os.path.join(PACKAGE_PATH, 'data')

_config = configparser.ConfigParser()
_config.read(CONFIG_PATH)

class DownloadRequiredError(BaseException):
    pass

class LISA_Core:
    '''
    The LISA object is the user's interface with the LISA algorithm. It holds the parameters specified by the user and 
    handles data loading from hdf5
    '''
    #class dictionary, converts cmd line alias to standard symbols for knockout methods
    factor_binding_technologies = dict(
        chipseq = 'ChIP-seq',
        motifs = 'Motifs'
    )

    def __init__(self, species, 
        cores = 1,
        isd_method = 'chipseq',
        num_datasets_selected = 10,
        verbose = True,
        log = None,
        rp_map = 'basic'
    ):
        num_datasets_selected_anova = 200
        #all paramter checking is done on LISA instantiation to ensure consistency between python module and cmd line usage
        self.isd_options = _config.get('lisa_params', 'isd_methods').split(',')
        self.background_options = _config.get('lisa_params', 'background_strategies').split(',')
        self.max_query_genes = int(_config.get('lisa_params', 'max_user_genelist_len'))
        
        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )

        self.num_datasets_selected_anova = num_datasets_selected_anova
        self.num_datasets_selected = num_datasets_selected

        self._set_cores(cores)

        assert( num_datasets_selected_anova > num_datasets_selected ), 'Anova must select more datasets than the regression model'
        assert( num_datasets_selected_anova > 0 and num_datasets_selected > 0 ), 'Number of datasets selected must be positive'
        assert(num_datasets_selected_anova < 500)
        assert(num_datasets_selected < 25)

        assert( isd_method in  self.isd_options ), 'ISD method must be \{{}}'.format(', '.join(self.isd_options))
        assert( species in ['mm10','hg38'] ), "Species must be either hg38 or mm10"
        
        self.isd_method = self.factor_binding_technologies[isd_method]
        self.species = species

        self.data_source = _config.get('paths', 'h5_path').format(package_path = PACKAGE_PATH, species = self.species)

        self.data_path = INSTALL_PATH

        #use provided log object or create a new one. Applications which call LISA may want to pass in their own log object for control of indentation
        if log is None:
            self.log = Log(sys.stderr, verbose = verbose)
        else:
            assert( isinstance(log, Log) )
            self.log = log
        
        #data has not yet been loaded
        self._is_loaded = False

        self.assays = []

        if isinstance(rp_map, str):
            rp_map_styles = _config.get('lisa_params','rp_map_styles').split(',')
            assert(rp_map in rp_map_styles), 'RP map must be numpy/scipy.sparse array, or be one of provided maps: {}'.format(','.join(rp_map_styles))
        else:
            assert( isinstance(rp_map, np.ndarry) or isinstance(rp_map, scipy.sparse)), 'RP map must be either numpy ndarry or scipy.sparse matrix'
        self.rp_map = rp_map

        if self.num_datasets_selected % self.cores != 0:
            self.log.append('''WARNING: You\'ve allocated {} cores with {} datasets selected.
To ensure maximum speed, allocate as many cores as datasets selected.
For better efficiency, make #datasets a multiple of #cores.'''.format(self.cores, self.num_datasets_selected))

    def _set_cores(self, cores):
        assert( isinstance(cores, int) and cores >= -1)
        #leave one core out for the rest of us
        max_cores = multiprocessing.cpu_count() - 1
        if cores <= -1:
            cores = max_cores
        #use the minimum number of cores between the user's specificaion, the number of datasets to be processed in parallel, and the number of cores on the machine.
        #this prevents LISA from using more resources than required.
        self.cores = min(cores, max_cores, self.num_datasets_selected)

    #____ DATA LOADING FUNCTIONS _____

    def _load_gene_info(self):

        all_genes = gene_selection.GeneSet()
        
        with open(_config.get('genes','master_gene_list').format(package_path = PACKAGE_PATH, species = self.species), 'r') as genes:
            all_genes.from_str(genes.read())

        with open(_config.get('genes','gene_locs').format(package_path = PACKAGE_PATH, species = self.species), 'r') as f:
            rp_map_locs = np.array([line.strip() for line in f.readlines()])

        self.all_genes, self.rp_map_locs = all_genes, rp_map_locs
        
    def _load_factor_binding_data(self, data_object):

        self.factor_dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
            .format(technology = self.isd_method)][...].astype(str))
        
        self.factor_binding = sparse.load_npz(_config.get('factor_binding','matrix').format(package_path = PACKAGE_PATH, species = self.species, technology = self.isd_method))

    def _load_rp_map(self):
        if isinstance(self.rp_map, str):
            #self.rp_map = sparse.load_npz(_config.get('RP_map','matrix').format(package_path = PACKAGE_PATH, species = self.species, style = self.rp_map)).tocsr()
            rp_map_name_params = dict(species = self.species, version = REQURED_DATASET_VERSION, name = self.rp_map)
            
            rp_map_path = _config.get('paths','rp_map').format(package_path = PACKAGE_PATH, **rp_map_name_params)
            
            if not os.path.exists(rp_map_path):
                try:
                    self.fetch_from_cistrome(_config.get('downloads', 'rp_maps')\
                        .format(**rp_map_name_params), rp_map_path, is_tar=False)
                except error.URLError as err:
                    self.log.append('Cannot connect to cistrome server, or cistrome server does not have the RP map requested. \nDefaulting to "basic" RP map instead.')
                    rp_map_name_params['name'] = 'basic'
                    rp_map_path = _config.get('paths','rp_map').format(package_path = PACKAGE_PATH, **rp_map_name_params)
            
            self.rp_map = sparse.load_npz(rp_map_path).tocsr()

        return self.rp_map

    def _load_data(self):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')

        with self.log.section('Loading data into memory (only on first prediction):') as log:

            self._load_gene_info()

            log.append('Loading regulatory potential map ...')
            self._load_rp_map()

            with h5.File(self.data_source, 'r') as data:

                with log.section('Loading binding data ...') as log:
                    self._load_factor_binding_data(data)

            self.factor_metadata = self.link_metadata(self.isd_method).select(self.factor_dataset_ids)

            self.log.append('Done!')

        self._is_loaded = True

    def link_metadata(self, technology):
        return Metadata(_config.get('metadata', technology).format(package_path = PACKAGE_PATH, species = self.species, technology = technology),
                    _config.get('metadata',technology+'_headers').split(','))

    #____ DATA DOWNLOADING FUNCTIONS ____

    def _check_for_data(self):
        if not (os.path.isdir(self.data_path) and os.path.isdir(os.path.join(self.data_path, self.species))):
            self.log.append('Data not found, must download from CistromeDB ...')
            raise DownloadRequiredError()

        else:
            try:
                with open(_config.get('paths','dataset_version').format(package_path = PACKAGE_PATH, species = self.species), 'r') as v:
                    dataset_version = v.read().strip()

                if not REQURED_DATASET_VERSION == dataset_version:
                    self.log.append('Dataset version mismatch, must download new dataset from CistromeDB ...')
                    raise DownloadRequiredError()

            except FileNotFoundError:
                raise DownloadRequiredError()

    @staticmethod
    def extract_tar(tarball, write_name, rm = False):
        with tarfile.open(tarball) as tar:
            tar.extractall(path = write_name)
        if rm:
            os.remove(tarball)

    def fetch_from_cistrome(self, url, write_name, is_tar = True):

        write_path = os.path.dirname(write_name)

        if not os.path.isdir(write_path):
                os.mkdir(write_path)

        filename, _ = request.urlretrieve(url)
        
        if is_tar:
            self.log.append('Extracting data ...')
            self.extract_tar(filename, write_name, rm = True)
        else:
            os.rename(filename, write_name)
        
        return True

    def get_dataset_url(self):
        return _config.get('downloads','dataset').format(species = self.species, version = REQURED_DATASET_VERSION)

    def download_data(self):

        with self.log.section('Grabbing {} data (~15 minutes):'.format(self.species)):
            
            self.log.append('Downloading from database ...')
                        
            dataset_url = self.get_dataset_url()

            try:
                self.fetch_from_cistrome(dataset_url, self.data_path, is_tar=True)
            except error.URLError as err:
                self.log.append('ERROR: Cannot connect to cistrome.org for data (usually due to security settings on some servers)!\nView github pages for manual dataset install instructions.')
                self.log.append(err)
                sys.exit()

    #____ GENE SELECTION FUNCTIONS ____
    def _preprocess_gene_list(self, genes):
        #make sure gene lists are composed of stripped strings
        return [str(gene).strip() for gene in genes if len(gene) > 0]
        
    def _sample_background_genes(self, query_genes, background_strategy = 'regulatory', num_background_genes = 3000, seed = None):

        assert( num_background_genes >= len(query_genes) and len(query_genes) <= 19000//2 ), "More query genes selected than background genes"

        background_candidates = self.all_genes.get_distinct_genes_by_symbol(excluding = query_genes.get_symbols())

        assert(len(background_candidates) >= num_background_genes), 'Number of background candidates must be greater than or equal number of genes to select as background.'
        #if no down-sampling is needed:
        if len(background_candidates) == num_background_genes:
            background_genes = background_candidates

        if background_strategy == 'regulatory':

            background_genes = background_candidates.sample_by_TAD(num_background_genes, seed = seed)

            if len(background_genes) > num_background_genes:
                background_genes = background_genes.random_sample(num_background_genes, seed = seed)

        elif background_strategy == 'random':
            background_genes = background_candidates.random_sample(num_background_genes, seed = seed)
        else:
            raise AssertionError('Background selection strategy {} not supported'.format(background_strategy))

        return background_genes

    def _get_query_and_background_genes(self, query_list, *, background_list = [], num_background_genes = 3000, background_strategy =
        'regulatory', seed = None):

        query_list = self._preprocess_gene_list(query_list)
        background_list = self._preprocess_gene_list(background_list)

        query_genes = self.all_genes.match_user_provided_genes(query_list)

        assert(20 <= len(query_genes) <= self.max_query_genes), 'User must provide list of 20 to {} unique genes. Provided {}'\
            .format(str(self.max_query_genes), str(len(query_genes)))

        if background_strategy == 'provided':
            
            background_matches = self.all_genes.match_user_provided_genes(background_list)
            background_genes = background_matches.get_distinct_genes_by_symbol(excluding = query_genes.get_symbols())
            assert( len(background_genes) > len(query_genes) ), 'Number of background genes must exceed number of query genes provided'
        
        else:
            background_genes = self._sample_background_genes(query_genes, background_strategy = background_strategy, 
                num_background_genes = num_background_genes, seed = seed)

        return query_genes, background_genes

    @staticmethod
    def create_label_dictionary(query_list, background_list):
        return { #this dict is used to create label vector for LISA algo
            gene.get_location() : int(i < len(query_list))
            for i, gene in enumerate(list(query_list) + list(background_list))
        }, dict( #return this dict for results (purely informational)
            query_symbols = query_list.get_symbols(),
            background_symbols = background_list.get_symbols(),
            query_locations = query_list.get_locations(),
            background_locations = background_list.get_locations(),
        )

    def _make_gene_mask(self, query_genes, background_genes):

        label_dict, gene_info_dict = self.create_label_dictionary(query_genes, background_genes)

        #Show user number of query and background genes selected
        self.log.append('Selected {} query genes and {} background genes.'\
            .format(str(len(gene_info_dict['query_symbols'])), str(len(gene_info_dict['background_symbols']))))

        #subset the rp map for those genes and make label vector
        gene_mask = np.isin(self.rp_map_locs, list(label_dict.keys()))
        #make label vector from label_dict based of gene_loc ordering in rp_map data
        label_vector = np.array([label_dict[gene_loc] for gene_loc in self.rp_map_locs[gene_mask]])
        
        return gene_mask, label_vector, gene_info_dict

    #KEYSTONE FUNCTION FOR GENE SELECTION
    def _choose_genes(self, query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = None):
        #validate user args
        if background_list is None:
            background_list = []
        assert( isinstance(num_background_genes, int) )
        assert( background_strategy in  self.background_options ), 'Background strategy must be in \{{}\}'.format(', '.join(self.background_options))
        assert( isinstance(query_list, Iterable) )
        assert( isinstance(background_list, Iterable))
        #check to make sure data is downloaded

        with self.log.section('Matching genes and selecting background ...'):

            if len(background_list) > 0 and background_strategy == 'provided':
                self.log.append('User provided background genes!')

            query_genes, background_genes = self._get_query_and_background_genes(query_list, background_list = background_list, 
                    num_background_genes = num_background_genes, background_strategy = background_strategy, seed = seed)
            
            gene_mask, label_vector, gene_info_dict = self._make_gene_mask(query_genes, background_genes)

        return gene_mask, label_vector, gene_info_dict

    #____ ASSAY INSTANTIATION AND CALLING
    def add_assay(self, assay):
        self.assays.append(assay)

    def _run_assays(self, gene_mask, label_vector):

        with h5.File(self.data_source, 'r') as data:
            #run each assay specified in the _load_data step. Each assay returns a p_value for each factor binding sample, and returns metadata
            #each assay takes the same arguments and makes the same returns but can perform different tests
            assay_pvals, assay_info = {},{}
            for assay in self.assays:
                assay_pvals[assay.technology + '_p_value'] = assay.predict(gene_mask, label_vector, data)
                assay_info[assay.technology] = assay.get_info()

        return assay_pvals, assay_info

    def _initialize_assays(self):
        raise NotImplementedError('_initialize_assays function is not implemented on base LISA object') 

    def _initialize_resources(self):
        try:
            self._check_for_data()
        except DownloadRequiredError:
            self.download_data()
            self._check_for_data()

        self.log.append('Using {} cores ...'.format(str(self.cores)))

        if not self._is_loaded:
            self._load_data()

    #___ RESULTS SUMMARIZATION AND FORMATTING____
    @staticmethod
    def _combine_tests(p_vals, weights=None):
        #https://arxiv.org/abs/1808.09011
        #combine p-values from all assays
        p_vals = np.array(p_vals).astype(np.float64)

        assert(len(p_vals.shape) == 2 ), 'P-values must be provided as matrix of (samples, multiple p-values)'
        assert(p_vals.shape[1] > 1), 'Must have multiple p-values to combine'

        if weights is None:
            weights = 1 / p_vals.shape[1]
        else:
            assert(p_vals.shape[1] == weights.shape[1])

        individual_statistics = np.where(p_vals > 1e-15, np.tan((0.5-p_vals) * np.pi), 1/p_vals/np.pi)
        test_statistic = np.sum(weights * individual_statistics, axis = 1)
        
        # lim x->inf of arctan(x) = pi/2
        # so for large x, arctan(x) ~ pi/2 - 1/x, as 1/x -> 0
        # substitute this arctan approx when x > 1e15

        too_large = test_statistic > 1e15
        combined_p_value = np.where(~too_large, 0.5 - np.arctan(test_statistic)/np.pi, 1/test_statistic/np.pi)

        return test_statistic, combined_p_value

    def _format_results(self, assay_pvals, assay_info, gene_info_dict, **kwargs):
        self.log.append('Formatting output ...')
        #instantiate a lisa results object
        results = LISA_Results.fromdict(
            **self.factor_metadata,
            **assay_pvals,
            **kwargs
        )

        results = results.sortby('summary_p_value', add_rank = True, reverse = False)
        
        self.log.append('Done!')
        
        #return model_metadata as big dictionary
        return results, dict(
            **gene_info_dict,
            **assay_info
        )

    def _summarize_p_values(self, assay_pvals):
        
        if len(assay_pvals) > 1:
            #combine p-values from all assays on a per-sample basis
            with self.log.section('Mixing effects using Cauchy combination ...'):

                aggregate_pvals = np.array(list(assay_pvals.values())).T
                
                combined_statistic, summary_p_value = self._combine_tests(aggregate_pvals)
        else:
            summary_p_value = np.array(list(assay_pvals.values())).reshape(-1)

        return summary_p_value


    #___ CLASS PREDICTION FUNCTION (USER INTERFACE) ___
    def predict(self, query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556):
        
        try:
            self._initialize_resources()

            if len(self.assays) == 0:
                #the initialize assay section is the mutable section for each interface
                self._initialize_assays()

            gene_mask, label_vector, gene_info_dict = self._choose_genes(query_list, background_list, background_strategy, num_background_genes, seed)

            assay_pvals, assay_info = self._run_assays(gene_mask, label_vector)

            summary_p_value = self._summarize_p_values(assay_pvals)

            return self._format_results(assay_pvals, assay_info, gene_info_dict, summary_p_value = summary_p_value)

        except (FileNotFoundError, OSError):
            raise DownloadRequiredError('Data is malformed or incomplete, run "lisa download [species]" to redownload dataset')