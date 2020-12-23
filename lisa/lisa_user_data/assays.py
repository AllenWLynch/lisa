from lisa.core.utils import LoadingBar
from lisa.core.assays import LISA_RP_Assay, get_delta_RP, get_deltaRP_activation
from multiprocessing import Pool
import numpy as np
from scipy import stats, sparse
import os

class ISD_Assay(LISA_RP_Assay):

    def __init__(self, region_scores, **kwargs):
        super().__init__(**kwargs, technology = 'Regions', metadata = [])
        self.region_scores = region_scores

    def predict(self, gene_mask, label_vector, *args,**kwargs):

        with self.log.section('Modeling insilico deletions:'):

            profile = self.region_scores[:, np.newaxis]

            factor_sums = self.factor_binding.sum(axis = 0)

            subset_rp_map = self.rp_map[gene_mask, :]
            self.log.append('Calcuating null RP model ...')
            null_rp_matrix = self.make_rp_matrix(profiles = profile, rp_map = subset_rp_map)

            self.log.append('Performing knockouts ...')
            delta_RP = np.squeeze(get_delta_RP(profile, self.factor_binding, subset_rp_map).transpose(0,2,1), axis=-1)

            self.log.append('Calculating Δ regulatory score ...')
            delta_reg_scores = get_deltaRP_activation(null_rp_matrix, delta_RP)

            self.log.append('Calculating p-values ...')
            p_vals = self.get_delta_RP_p_value(delta_reg_scores, label_vector)

        self.log.append('Done!')

        return p_vals

    def get_info(self):
        return dict()