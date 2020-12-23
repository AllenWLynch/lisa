
from lisa.core.lisa_core import LISA_Core, PACKAGE_PATH, REQURED_DATASET_VERSION
from lisa.core.lisa_core import CONFIG_PATH as base_config_path
import lisa.lisa_public_data.assays as assays
from lisa.lisa_public_data.models import LR_BinarySearch_SampleSelectionModel, LR_ChromatinModel
import numpy as np
from scipy import sparse
import multiprocessing
import os
import configparser
"""
LISA_Core implements the main methods for results formatting and data loading that make
up a standard LISA application. Extensions of the core method may instantiate different assays
depending on the data available to asses TF influence in different ways.

If only expression data is available, the base LISA interface can be instantiated, which 
will use public data to estimate TF influence.
"""

CONFIG_PATH = os.path.join(os.path.dirname(__file__),'config.ini')
_config = configparser.ConfigParser()
_config.read([base_config_path, CONFIG_PATH])

class LISA(LISA_Core):

    def __init__(self, species, config = _config, rp_map = 'basic', assays = ['Direct','H3K27ac','DNase'], num_datasets_selected = int(_config.get('lisa_params', 'num_datasets_selected')), cores = 1, **kwargs):
        super().__init__(species, config, **kwargs)
        assert(len(assays) > 0), 'Must provide at least one assay to run.'
        assert(all([assay in self._config.get('lisa_params','assays').split(',') for assay in assays])), 'An assay chosen by the user is not a valid choice: \{{}}'.format(self._config.get('lisa_params','assays'))
        #assert(self.rp_map == 'basic'), 'For base LISA predictor, rp map must be "basic".'

        num_datasets_selected_anova = int(self._config.get('lisa_params','num_datasets_selected_anova'))

        assert( isinstance(num_datasets_selected, int) )
        assert( isinstance(num_datasets_selected_anova, int) )

        self.num_datasets_selected_anova = num_datasets_selected_anova
        self.num_datasets_selected = num_datasets_selected

        self.cores = self._set_cores(cores)

        assert( num_datasets_selected_anova > num_datasets_selected ), 'Anova must select more datasets than the regression model'
        assert( num_datasets_selected_anova > 0 and num_datasets_selected > 0 ), 'Number of datasets selected must be positive'
        assert(num_datasets_selected_anova < 500)
        assert(num_datasets_selected < 25)

        if self.num_datasets_selected % self.cores != 0:
            self.log.append('''WARNING: You\'ve allocated {} cores with {} datasets selected.
To ensure maximum speed, allocate as many cores as datasets selected.
For better efficiency, make #datasets a multiple of #cores.'''.format(self.cores, self.num_datasets_selected))

        self.data_source = _config.get('paths', 'h5_path').format(package_path = PACKAGE_PATH, species = self.species)

        self.schedule_assays = sorted(list(set(assays)))

        if isinstance(rp_map, str):
            rp_map_styles = self._config.get('lisa_params','rp_map_styles').split(',')
            assert(rp_map in rp_map_styles), 'RP map must be numpy/scipy.sparse array, or be one of provided maps: {}'.format(','.join(rp_map_styles))
        else:
            assert( isinstance(rp_map, np.ndarry) or isinstance(rp_map, scipy.sparse)), 'RP map must be either numpy ndarry or scipy.sparse matrix'
        self.rp_map = rp_map
        self.generate_rp_matrix = rp_map != 'basic'

    def _set_cores(self, cores):
        assert( isinstance(cores, int) and cores >= -1)
        #leave one core out for the rest of us
        max_cores = multiprocessing.cpu_count() - 1
        if cores <= -1:
            cores = max_cores
        #use the minimum number of cores between the user's specificaion, the number of datasets to be processed in parallel, and the number of cores on the machine.
        #this prevents LISA from using more resources than required.
        self.cores = min(cores, max_cores, self.num_datasets_selected)
        return self.cores

    # Change this part
    def _load_factor_binding_data(self):
        super()._load_factor_binding_data('1000')
        
    def _load_rp_map(self):
        if isinstance(self.rp_map, str):
            #self.rp_map = sparse.load_npz(self._config.get('RP_map','matrix').format(package_path = PACKAGE_PATH, species = self.species, style = self.rp_map)).tocsr()
            rp_map_name_params = dict(species = self.species, version = REQURED_DATASET_VERSION, name = self.rp_map)
            
            rp_map_path = self._config.get('paths','rp_map').format(package_path = PACKAGE_PATH, **rp_map_name_params)
            
            if not os.path.exists(rp_map_path):
                try:
                    self.fetch_from_cistrome(self._config.get('downloads', 'rp_maps')\
                        .format(**rp_map_name_params), rp_map_path, is_tar=False)
                except error.URLError as err:
                    self.log.append('Cannot connect to cistrome server, or cistrome server does not have the RP map requested. \nDefaulting to "basic" RP map instead.')
                    rp_map_name_params['name'] = 'basic'
                    rp_map_path = self._config.get('paths','rp_map').format(package_path = PACKAGE_PATH, **rp_map_name_params)
            
            self.rp_map = sparse.load_npz(rp_map_path).tocsr()

        return self.rp_map


    def get_factor_gene_mask(self):

        factor_genes = self.all_genes.match_user_provided_genes(self.factor_metadata['factor'])

        loc_symbol_dict = dict(zip(factor_genes.get_locations(), factor_genes.get_symbols()))
        self.factor_gene_mask = np.isin(self.rp_map_locs, factor_genes.get_locations())
        #make label vector from label_dict based of gene_loc ordering in rp_map data
        self.factor_mask_keys = [loc_symbol_dict[gene_loc] for gene_loc in self.rp_map_locs[self.factor_gene_mask]]

        return self.factor_gene_mask, self.factor_mask_keys

    def _format_results(self, assay_pvals, assay_info, gene_info_dict, **kwargs):
        
        new_columns = {}
        for assay, info in assay_info.items():
            if 'factor_acc_z_scores' in info:
                z_scores = info.pop('factor_acc_z_scores')
                z_dict = dict(zip(self.factor_mask_keys, z_scores))
                new_columns[assay + 'factor_accessibility_z_score'] = [z_dict.get(factor, 'NA') for factor in self.factor_metadata['factor']]

        return super()._format_results(assay_pvals, assay_info, gene_info_dict, **new_columns, **kwargs)

    def _initialize_assays(self, **assay_kwargs):
        #Add assays to LISA's steps. Each assay follows the same instantiation and prediction calling, making them modular and substitutable.
        #Adding an assay loads the required data for that assay

        try: # if factor gene mask not instantiated
            self.factor_gene_mask
        except AttributeError:
            self.get_factor_gene_mask()

        for assay in self.schedule_assays:
            if assay == 'Direct':
                self.add_assay(
                    assays.PeakRP_Assay(
                        technology = self.isd_method, **assay_kwargs, 
                        metadata = self.link_metadata(self.isd_method),
                        generate_rp_matrix = self.generate_rp_matrix
                    )
                )
            elif assay == 'DNase' or assay == 'H3K27ac':
                self.add_assay(
                    assays.Accesibility_Assay(technology = assay, **assay_kwargs,
                        metadata = self.link_metadata(assay), factor_gene_mask = self.factor_gene_mask,
                        selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2'), 
                        cores = self.cores, generate_rp_matrix = self.generate_rp_matrix
                    )
                )
            else:
                raise AssertionError('Invalid assay encountered: {}'.format(str(assay)))