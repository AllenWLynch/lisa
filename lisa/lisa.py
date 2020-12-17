
from lisa.lisa_core import LISA_Core, _config
import lisa.assays as assays
from lisa.models import LR_BinarySearch_SampleSelectionModel
from lisa.models import LR_ChromatinModel
import numpy as np
import multiprocessing
"""
LISA_Core implements the main methods for results formatting and data loading that make
up a standard LISA application. Extensions of the core method may instantiate different assays
depending on the data available to asses TF influence in different ways.

If only expression data is available, the base LISA interface can be instantiated, which 
will use public data to estimate TF influence.
"""
class LISA(LISA_Core):

    def __init__(self, *args, assays = ['Direct','H3K27ac','DNase'], num_datasets_selected = 10, cores = 1, **kwargs):
        super().__init__(*args, **kwargs)
        assert(len(assays) > 0), 'Must provide at least one assay to run.'
        assert(all([assay in _config.get('lisa_params','assays').split(',') for assay in assays])), 'An assay chosen by the user is not a valid choice: \{{}}'.format(_config.get('lisa_params','assays'))
        assert(self.rp_map == 'basic'), 'For base LISA predictor, rp map must be "basic".'

        num_datasets_selected_anova = 200

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

        self.schedule_assays = sorted(list(set(assays)))


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

    def _initialize_assays(self):
        #Add assays to LISA's steps. Each assay follows the same instantiation and prediction calling, making them modular and substitutable.
        #Adding an assay loads the required data for that assay

        #The first assay to be loaded must be the ChIP-seq or Motif direct knockout, because this assay supplies the 
        #metadata for the final output table since it holds information on all factor binding samples
        assay_kwargs = dict(config = _config, cores = self.cores, log = self.log)

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
                    )
                )
            elif assay == 'DNase' or assay == 'H3K27ac':
                self.add_assay(
                    assays.Accesibility_Assay(technology = assay, **assay_kwargs,
                        metadata = self.link_metadata(assay), factor_gene_mask = self.factor_gene_mask,
                        rp_map = self.rp_map, factor_binding = self.factor_binding,
                        selection_model = LR_BinarySearch_SampleSelectionModel(self.num_datasets_selected_anova, self.num_datasets_selected),
                        chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-2,4.1,0.5))}, penalty = 'l2')
                    )
                )
            else:
                raise AssertionError('Invalid assay encountered: {}'.format(str(assay)))