
import numpy as np
from scipy import sparse, stats
from models import ChromatinModel
from utils import LoadingBar
from multiprocessing import Pool

def mannu_test_function(x):
    query, background = x
    try:
        return stats.mannwhitneyu(query, background, alternative = 'greater')
    #if all values in query and background are equal (no knockouts), throws value error
    except ValueError:
        #catch, return none for test statistic, 1.0 for p-val
        return (None, 1.0)

def get_delta_RP(profile, binding_data, rp_map):
    '''
    profile: chromatin profile, bin x 1 array of accessibility at each genomic region
    rp_map: sparse bin x gene matrix mapping a genomic region to a gene based on RP proximity
    binding_data: bin x dataset matrix, with binary hits to show TF binding occurs in a location

    the line below defines the mapping of chromatin accessibility at TF binding sites to genes.

    returns: gene x TF x 1 matrix of delta RPs
    '''
    return np.array(rp_map.dot(binding_data.astype(np.bool).multiply(profile.reshape((-1,1)))).todense())[:,np.newaxis,:]

#distributes arguments for get_delta_RP function, to be used by multiprocessing module
def delta_RP_wrapper(x):
    return get_delta_RP(*x)

#Gotta fix this up
def get_delta_RP_p_value(delta_regulation_score, labels, cores = 1):
    '''
    delta_regulation_score: gene x TF, model output of delta-RP matrix. more purturbation of genes of interest correspond with higher delta regulation score
    labels: query vs. background labels for each gene in delta-RP matrix
    '''

    labels = np.array(list(labels))

    #seperate matrix into query and background sets
    query_delta = delta_regulation_score[labels.astype(np.bool)]
    background_delta = delta_regulation_score[~labels.astype(np.bool)]

    #for each TF, calculate p-value of difference in delta-R distributions between query and background genes
    test_parameters = list(zip(query_delta.T, background_delta.T))

    if cores == 1:
        p_vals = [
            mannu_test_function((q,b)) for q,b in test_parameters
        ]

    else:
        with Pool(cores) as p:
            p_vals = p.map(mannu_test_function, test_parameters)

    test_statistic, p_values = list(zip(*p_vals))

    return test_statistic, p_values


class KnockoutGenerator:

    def __init__(self, profiles, factor_binding, rp_map):
        self.factor_binding = factor_binding
        self.rp_map = rp_map
        self.profiles = profiles

    def __iter__(self):

        for profile in self.profiles.T:
            yield profile, self.factor_binding, self.rp_map

def get_insilico_deletion_significance(*, log, accessibility_profiles, factor_binding, rp_map, chromatin_model, labels, cores): #enforced kwargs
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

    log.append('Performing in-silico knockouts ...')
    #loop through datasets and find delta-RP from all knockouts
    '''knockouts = []
    loadingbar = LoadingBar('Performing in-silico knockouts', accessibility_profiles.shape[1], 20, cold_start = True)
    log.append(loadingbar, update_line = True)
    for profile in accessibility_profiles.T:
        knockouts.append(get_delta_RP(profile, factor_binding, rp_map))
        log.append(loadingbar, update_line = True)'''
    with Pool(cores) as p:
        knockouts = p.map(delta_RP_wrapper, iter(KnockoutGenerator(accessibility_profiles, factor_binding, rp_map)))
     
    #concatenate datasets to for gene x TF x datasets shaped matrix
    log.append('Calculating Δ regulatory score ...')
    datacube = np.concatenate(knockouts, axis = 1)

    num_genes, num_samples, num_TFs = datacube.shape

    delta_regulation_score = chromatin_model.get_deltaRP_activation(datacube)
    
    assert(delta_regulation_score.shape == (num_genes, num_TFs))
    
    log.append('Caculating p-values ...')
    #apply the wilcoxon test on each TF
    statistic, p_vals = get_delta_RP_p_value(delta_regulation_score, labels, cores=cores)

    return statistic, p_vals, datacube

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

