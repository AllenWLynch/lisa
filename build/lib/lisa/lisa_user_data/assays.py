from lisa.core.utils import LoadingBar
from lisa.core.assays import LISA_RP_Assay, get_delta_RP, get_deltaRP_activation
from multiprocessing import Pool
import numpy as np
from scipy import stats, sparse
import os

class ISD_Assay(LISA_RP_Assay):

    def predict(self, gene_mask, label_vector,*,region_scores,**kwargs):

        with self.log.section('Modeling insilico deletions:'):

            profile = (region_scores/region_scores.sum() * 1e5)[:,np.newaxis]

            subset_profile, subset_factor_binding, subset_rp_map = \
                self.make_subset_rp_map(rp_map = self.rp_map, factor_binding = self.factor_binding,
                    gene_mask = gene_mask, accessibility = profile)

            self.log.append('Calcuating null RP model ...')

            null_rp_matrix = self.make_rp_matrix(profiles = subset_profile, rp_map = subset_rp_map)

            self.log.append('Performing knockouts ...')
            delta_RP = get_delta_RP(subset_profile, subset_factor_binding, subset_rp_map)

            self.log.append('Calculating Δ regulatory score ...')
            delta_reg_scores = get_deltaRP_activation(null_rp_matrix, delta_RP)

            self.log.append('Calculating p-values ...')
            p_vals = self.get_delta_RP_p_value(delta_reg_scores, label_vector)

            query_reg_score_matrix, background_reg_score_matrix, top_factor_idx = self.get_delta_reg_score_matrix(p_vals, delta_reg_scores, label_vector)

            '''self.debug = dict(
                profile = profile,
                gene_mask = gene_mask,
                label_vector = label_vector,
                subset_rp_map = subset_rp_map,
                null_rp_matrix = null_rp_matrix,
                delta_RP_matrix = delta_RP,
                delta_reg_scores = delta_reg_scores,
            )'''

        self.log.append('Done!')

        return p_vals, dict(
            dataset_ids = list(np.array(self.factor_dataset_ids)[top_factor_idx]),
            query_reg_scores = query_reg_score_matrix.tolist(),
            background_reg_scores = background_reg_score_matrix.tolist()
        )
        