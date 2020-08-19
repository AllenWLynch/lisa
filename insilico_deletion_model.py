
import numpy as np
from scipy import sparse, stats
from models import ChromatinModel
from utils import LoadingBar

def get_delta_RP(profile, binding_data, rp_map):
    '''
    profile: chromatin profile, bin x 1 array of accessibility at each genomic region
    rp_map: sparse bin x gene matrix mapping a genomic region to a gene based on RP proximity
    binding_data: bin x dataset matrix, with binary hits to show TF binding occurs in a location

    the line below defines the mapping of chromatin accessibility at TF binding sites to genes.

    returns: gene x TF x 1 matrix of delta RPs
    '''
    return np.array(rp_map.dot(binding_data.astype(np.bool).multiply(profile.reshape((-1,1)))).todense())[:,np.newaxis,:]


def get_delta_RP_p_value(delta_regulation_score, labels):
    '''
    delta_regulation_score: gene x TF, model output of delta-RP matrix. more purturbation of genes of interest correspond with higher delta regulation score
    labels: query vs. background labels for each gene in delta-RP matrix
    '''

    labels = np.array(list(labels))

    #seperate matrix into query and background sets
    query_delta = delta_regulation_score[labels.astype(np.bool)]
    background_delta = delta_regulation_score[~labels.astype(np.bool)]

    #for each TF, calculate p-value of difference in delta-R distributions between query and background genes
    p_vals = []
    for i, (q, b) in enumerate(zip(query_delta.T, background_delta.T)):
        try:
            p_vals.append(
                stats.mannwhitneyu(q, b, alternative = 'greater')
            )
        except ValueError as err:
            # Occurs when the values of query and background genes are the same. This means that a chip-seq sample did not affect any genes.
            #print('Null on dataset: {}'.format(str(i)))
            p_vals.append((None, 1.0))

    test_statistic, p_values = list(zip(*p_vals))

    return test_statistic, p_values


def get_insilico_deletion_significance(*, log, accessibility_profiles, factor_binding, rp_map, chromatin_model, labels, normalizer, full_return = False): #enforced kwargs
    """
    log: log object for printing
    accessibility_profiles: (bins, num_datasets): list of chromatin-accessiblity datasets on which to analyze the effects of ISD
    factor_binding: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
    rp_map: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)
    chromatin_model: a ChromatinModel object that implements fit and get_deltaRP_activation. Purturbation of this model summarizes effect of ISD
    labels: query vs. background labels

    returns:
    p-value of wilcoxon test for each TF
    """
    assert( isinstance(chromatin_model, ChromatinModel) )

    loadingbar = LoadingBar('Performing in-silico knockouts', accessibility_profiles.shape[1], 20, cold_start = True)

    log.append(loadingbar, update_line = True)
    
    #loop through datasets and find delta-RP from all knockouts
    knockouts = []
    for profile in accessibility_profiles.T:
        knockouts.append(get_delta_RP(profile, factor_binding, rp_map))
        log.append(loadingbar, update_line = True)
     
    #concatenate datasets to for gene x TF x datasets shaped matrix
    log.append('Calculating Δ regulatory score ...')
    datacube = np.concatenate(knockouts, axis = 1)

    num_genes, num_samples, num_TFs = datacube.shape

    datacube = normalizer(datacube)

    #isd_conditions = datacube.transpose(0,2,1).reshape((-1, num_samples))

    #use chromatin model to summarize change in delta-RP across all datasets
    delta_regulation_score = chromatin_model.get_deltaRP_activation(datacube)
    
    assert(delta_regulation_score.shape == (num_genes, num_TFs))
    
    log.append('Caculating p-values ...')
    #apply the wilcoxon test on each TF
    statistic, p_vals = get_delta_RP_p_value(delta_regulation_score, labels)

    if full_return:
        return statistic, p_vals, isd_conditions, delta_regulation_score
    
    return statistic, p_vals
    

def cauchy_combination(p_vals, wi=None):
    #https://arxiv.org/abs/1808.09011
    p_vals = np.array(p_vals).astype(np.float64)

    assert(len(p_vals.shape) == 2 ), 'P-values must be provided as matrix of (samples, multiple p-values)'
    assert(p_vals.shape[1] > 1), 'Must have multiple p-values to combine'

    #clip p-values at minimum float64 value to prevent underflow
    p_vals = np.clip(p_vals, np.nextafter(0, 1, dtype = np.float64), 1.0)
    if wi is None:
        wi = np.ones((1, p_vals.shape[1])) / p_vals.shape[1]
    else:
        assert(p_vals.shape[1] == weights.shape[1])

    test_statistic = np.sum(wi * np.tan((0.5-p_vals) * np.pi), axis = 1)
    combined_p_value = 0.5 - np.arctan(test_statistic)/np.pi

    return test_statistic, combined_p_value

