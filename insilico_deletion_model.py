
import numpy as np
from scipy import sparse, stats
from models import ChromatinModel
from utils import LoadingBar

def get_delta_RP(profile, binding_data, rp_map):

    return np.array(rp_map.dot(binding_data.multiply(profile.reshape((-1,1)))).todense())[:,:,np.newaxis]


def get_insilico_deletion_significance(*, log, accessibility_profiles, factor_binding, rp_map, chromatin_model, labels, normalizer, alternative = 'greater'): #enforced kwargs
    """
    log: log object for printing
    accessibility_profiles: (bins, num_datasets): list of chromatin-accessiblity datasets on which to analyze the effects of ISD
    factor_binding: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
    rp_map: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)
    chromatin_model: a ChromatinModel object that implements fit and get_deltaRP_activation. Purturbation of this model summarizes effect of ISD
    labels: query vs. background labels
    alternative: alternative hypothesis type. Usually, testing to see if delta regulatory score is greater for query genes

    returns:
    p-value of wilcoxon test for each TF
    """
    #get the delta reg score for each TF on each gene, summarized for each dataset
    assert( isinstance(chromatin_model, ChromatinModel) )

    loadingbar = LoadingBar('Performing in-silico knockouts', accessibility_profiles.shape[1], 20, cold_start = True)

    log.append(loadingbar, update_line = True)
    
    knockouts = []
    for profile in accessibility_profiles.T:
        knockouts.append(get_delta_RP(profile, factor_binding, rp_map))
        log.append(loadingbar, update_line = True)
     
    log.append('Calculating Δ regulatory score ...')
    datacube = np.concatenate(knockouts, axis = 2)

    num_genes, num_TFs, num_samples = datacube.shape

    isd_conditions = datacube.reshape((-1, num_samples)) # (genes + TFs) x samples

    isd_conditions = normalizer(isd_conditions)

    delta_regulation_score = chromatin_model.get_deltaRP_activation(isd_conditions).reshape(num_genes, num_TFs)
    
    log.append('Caculating p-values ...')
    #apply the wilcoxon test on each TF
    query_delta = delta_regulation_score[labels.astype(np.bool)]
    background_delta = delta_regulation_score[~labels.astype(np.bool)]
    
    p_vals = []
    for q, b in zip(query_delta.T, background_delta.T):
        try:
            p_vals.append(
                stats.mannwhitneyu(q, b, alternative = alternative)
            )
        except ValueError:
            p_vals.append(None)

    test_statistic, p_values = list(zip(*p_vals))

    return test_statistic, p_values


def cauchy_combination(p_vals, wi=None):
    #https://arxiv.org/abs/1808.09011
    p_vals = np.array(p_vals, np.float64)

    assert(len(p_vals.shape) == 2 ), 'P-values must be provided as matrix of (samples, multiple p-values)'
    assert(p_vals.shape[1] > 1), 'Must have multiple p-values to combine'

    if np.any(p_vals <= 1e-15): # np.finfo(np.float64) 1e-15
        raise NotImplementedError('P-values underflowing float64 datatype, install mpmath to get accurate results.')
    else:
        if wi is None:
            wi = np.ones((1, p_vals.shape[1])) / p_vals.shape[1]
        else:
            assert(p_vals.shape[1] == weights.shape[1])

        test_statistic = np.sum(wi * np.tan((0.5-p_vals) * np.pi), axis = 1)
        combined_p_value = 0.5 - np.arctan(test_statistic)/np.pi

        return test_statistic, combined_p_value

