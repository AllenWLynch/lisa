
#lisa modules
import lisa.gene_selection as gene_selection
import lisa.assays as assays
from lisa.utils import LoadingBar, Log, LISA_Results
from lisa.models import LR_BinarySearch_SampleSelectionModel
from lisa.models import LR_ChromatinModel
from lisa._version import __version__

#standard library
import configparser
from collections.abc import Iterable
import sys
import json
import multiprocessing
from urllib import request
import tarfile
import os

#required
import numpy as np
import h5py as h5
from scipy import sparse

PACKAGE_PATH = os.path.dirname(__file__)
CONFIG_PATH = os.path.join(PACKAGE_PATH, 'config.ini')
REQURED_DATASET_VERSION = '.'.join(__version__.split('.')[:2])

_config = configparser.ConfigParser()
_config.read(CONFIG_PATH)

class DownloadRequiredError(BaseException):
    pass

class LISA:
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
        num_datasets_selected_anova = 200,
        num_datasets_selected = 10,
        verbose = True,
        oneshot = False,
        log = None,
    ):
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

        self.data_path = os.path.join(PACKAGE_PATH, 'data')

        #use provided log object or create a new one. Applications which call LISA may want to pass in their own log object for control of indentation
        if log is None:
            self.log = Log(sys.stderr, verbose = verbose)
        else:
            assert( isinstance(log, Log) )
            self.log = log
        
        #data has not yet been loaded
        self._is_loaded = False

        self.assays = []

        assert( isinstance(oneshot, bool ))
        self.oneshot = oneshot

        self.used_oneshot = False


    def _set_cores(self, cores):
        assert( isinstance(cores, int) and cores >= -1)
        #leave one core out for the rest of us
        max_cores = multiprocessing.cpu_count() - 1
        if cores <= -1:
            cores = max_cores
        #use the minimum number of cores between the user's specificaion, the number of datasets to be processed in parallel, and the number of cores on the machine.
        #this prevents LISA from using more resources than required.
        self.cores = min(cores, max_cores, self.num_datasets_selected)

    def _preprocess_gene_list(self, genes):
        #make sure gene lists are composed of stripped strings
        return [str(gene).strip() for gene in genes if len(gene) > 0]

    def _load_gene_info(self):

        all_genes = gene_selection.GeneSet()
        
        with open(_config.get('genes','master_gene_list').format(package_path = PACKAGE_PATH, species = self.species), 'r') as genes:
            all_genes.from_str(genes.read())

        with open(_config.get('genes','gene_locs').format(package_path = PACKAGE_PATH, species = self.species), 'r') as f:
            rp_map_locs = np.array([line.strip() for line in f.readlines()])

        return all_genes, rp_map_locs
        

    def _load_factor_binding_data(self, data_object):

        self.factor_dataset_ids = list(data_object[_config.get('accessibility_assay','reg_potential_dataset_ids')\
            .format(technology = self.isd_method)][...].astype(str))
        
        return sparse.load_npz(_config.get('factor_binding','matrix').format(package_path = PACKAGE_PATH, species = self.species, technology = self.isd_method))

    def _load_rp_map(self, data_object):

        return sparse.load_npz(_config.get('RP_map','matrix').format(package_path = PACKAGE_PATH, species = self.species)).tocsr()


    def add_assay(self, data_object, gene_mask, assay):
        assay.load(data_object, gene_mask=gene_mask, oneshot=self.oneshot)
        self.assays.append(assay)


    def _load_data(self, gene_mask = None):

        if self._is_loaded:
            raise AssertionError('Data is already loaded')
        
        #if not already set, load this data
        try:
            self.rp_map_locs
        except AttributeError:
            self.all_genes, self.rp_map_locs = self._load_gene_info()
        
        with self.log.section('Loading data into memory (only on first prediction):') as log:

            with h5.File(self.data_source, 'r') as data:

                with log.section('Loading binding data ...') as log:
                    self.factor_binding = self._load_factor_binding_data(data)

                log.append('Loading regulatory potential map ...')
                self.rp_map = self._load_rp_map(data)

                #Add assays to LISA's steps. Each assay follows the same instantiation and prediction calling, making them modular and substitutable.
                #Adding an assay loads the required data for that assay

                #The first assay to be loaded must be the ChIP-seq or Motif direct knockout, because this assay supplies the 
                #metadata for the final output table since it holds information on all factor binding samples
                self.add_assay(data, gene_mask,
                    assays.PeakRP_Assay(self.isd_method, _config, self.cores, self.log,
                    assays.LISA_RP_Assay.invitro_metadata if self.isd_method == 'ChIP-seq' else assays.LISA_RP_Assay.insilico_metadata, 
                    _config.get('metadata',self.isd_method).format(package_path = PACKAGE_PATH, species = self.species, technology =self.isd_method),
                    )
                )

                self.add_assay(data, gene_mask,
                    assays.Accesibility_Assay('DNase', _config, self.cores, self.log, 
                        _config.get('metadata','DNase').format(package_path = PACKAGE_PATH, species = self.species, technology = 'DNase'),
                        rp_map = self.rp_map, factor_binding = self.factor_binding,
                        selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                    )
                )

                self.add_assay(data, gene_mask,
                    assays.Accesibility_Assay('H3K27ac', _config, self.cores, self.log,
                        _config.get('metadata', 'H3K27ac').format(package_path = PACKAGE_PATH, species = self.species, technology = 'H3K27ac'),
                        rp_map = self.rp_map, factor_binding = self.factor_binding,
                        selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                    )
                )

            self.log.append('Done!')

        self._is_loaded = True

    
    def _check_for_data(self):
        if not (os.path.isdir(self.data_path) and os.path.isdir(os.path.join(self.data_path, self.species))):
            self.log.append('Data not found, must download from CistromeDB ...')
            raise DownloadRequiredError()

        else:
            with open(_config.get('paths','dataset_version').format(package_path = PACKAGE_PATH, species = self.species), 'r') as v:
                dataset_version = v.read().strip()

            if not REQURED_DATASET_VERSION == dataset_version:
                self.log.append('Dataset version mismatch, must download new dataset from CistromeDB ...')


    def _download_data(self):

        print(__file__)
        print('here')

        assert(False)

        with self.log.section('Grabbing {} data (~15 minutes):'.format(self.species)):
            
            self.log.append('Downloading from database ...')
            #make data directory if does not exist
            if not os.path.isdir(self.data_path):
                os.mkdir(self.data_path)
            
            #http://cistrome.org/~alynch/data/genes.tar.gz
            download_dataset = _config.get('downloads','{species}_{version}'.format(species = self.species, version = REQURED_DATASET_VERSION))

            filename, _ = request.urlretrieve(
                download_dataset, 
                os.path.join(self.data_path, self.species + '_data.tar.gz')
            )
            
            self.log.append('Extracting data ...')

            with tarfile.open(filename) as tar:
                tar.extractall(path = self.data_path)

            os.remove(filename)

            with open(_config.get('paths','dataset_version').format(package_path = PACKAGE_PATH, species = self.species), 'w') as v:
                v.write(REQURED_DATASET_VERSION)

            self.log.append('Done!\n')

    @staticmethod
    def _combine_tests(p_vals, weights=None):
        #https://arxiv.org/abs/1808.09011
        #combine p-values from all assays
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
    def _make_gene_mask(self, query_list, background_list, num_background_genes, background_strategy, seed = None):
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
            max_query_genes = self.max_query_genes, background_list = background_list, seed = seed)

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
            seed = None,
        ):

        if self.oneshot and self.used_oneshot:
            raise AssertionError('When instantiated in one-shot, model cannot be used for multiple predictions')

        if background_list is None:
            background_list = []
        assert( isinstance(num_background_genes, int) )
        assert( background_strategy in  self.background_options ), 'Background strategy must be in \{{}\}'.format(', '.join(self.background_options))
        assert( isinstance(query_list, Iterable) )
        assert( isinstance(background_list, Iterable))
 
        try:
            self._check_for_data()
        except DownloadRequiredError:
            self._download_data()
            self._check_for_data()

        try:
            self.log.append('Using {} cores ...'.format(str(self.cores)))

            if not self._is_loaded:
                self.all_genes, self.rp_map_locs = self._load_gene_info()

            with self.log.section('Matching genes and selecting background ...'):

                if len(background_list) > 0 and background_strategy == 'provided':
                    self.log.append('User provided background genes!')
                
                gene_mask, label_vector, gene_info_dict = self._make_gene_mask(query_list, background_list, num_background_genes, background_strategy, seed = seed)

            #based on genes selected, only subset of data can be loaded
            if not self._is_loaded:
                self._load_data(gene_mask)
                #this method sets _is_loaded to True

            with h5.File(self.data_source, 'r') as data:
                
                #run each assay specified in the _load_data step. Each assay returns a p_value for each factor binding sample, and returns metadata
                #each assay takes the same arguments and makes the same returns but can perform different tests
                assay_pvals, assay_info = {},{}
                for assay in self.assays:
                    assay_pvals[assay.technology + '_p_value'] = assay.predict(gene_mask, label_vector, data_object = data)
                    assay_info[assay.technology + '_model_info'] = assay.get_info()
        
        except (FileNotFoundError, OSError):
            raise DownloadRequiredError('Data is malformed or incomplete, run "lisa download [species]" to redownload dataset')

        #combine p-values from all assays on a per-sample basis
        with self.log.section('Mixing effects using Cauchy combination ...'):

            aggregate_pvals = np.array(list(assay_pvals.values())).T
            
            combined_p_values = self._combine_tests(aggregate_pvals)

        self.log.append('Formatting output ...')
        #instantiate a lisa results object
        results = LISA_Results.fromdict(
            **self.assays[0].get_metadata(),
            combined_p_value = combined_p_values,
            combined_p_value_adjusted = list(np.minimum(np.array(combined_p_values) * len(combined_p_values), 1.0)),
            **assay_pvals
        )

        results = results.sortby('combined_p_value_adjusted', add_rank = True)
        
        self.log.append('Done!')
        
        #if using in one-shot mode, prevent user from calling "predict" on this object again
        self.used_oneshot = True
        #return model_metadata as big dictionary
        return results, dict(
            **gene_info_dict,
            **assay_info
        )