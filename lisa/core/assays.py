import numpy as np
from scipy import stats, sparse
import os
from multiprocessing import Pool

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

def mannu_test_function(x):
    query, background = x
    try:
        return stats.mannwhitneyu(query, background, alternative = 'greater')
    #if all values in query and background are equal (no knockouts), throws value error
    except ValueError:
        #catch, return none for test statistic, 1.0 for p-val
        return (None, 1.0)

def get_deltaRP_activation(rp_0, rp_knockout):
    """
    rp_knockout: is a datacube of shape (genes, samples, TFs),
    this method must implement a transformation into a genes x TFs matrix, sumarrizing the dataset-axis effects
    """
    #subtract define deltaX to be the log2 of the fraction of knocked-out RP
    deltaX = np.log2(rp_0 - rp_knockout + 1) - np.log2(rp_0 + 1)
    
    #flip sign so that more knockout = more deltaR
    return -deltaX


class LISA_RP_Assay:

    def __init__(self, *, technology, config, log, metadata, rp_map, factor_binding):
        self.config = config
        self.technology = technology
        self.log = log
        self.loaded = False
        self.metadata = metadata
        self.rp_map = rp_map
        self.factor_binding = factor_binding

    @staticmethod
    def make_subset_rp_map(*,rp_map, factor_binding, gene_mask, accessibility):
        
        bin_mask = np.squeeze(np.array(rp_map[gene_mask, : ].tocsc().sum(axis = 0) > 0))
        #subset rp_map and factor hits on bins with RP > 0
        subset_factor_binding = factor_binding[bin_mask, :]
        
        subset_rp_map = rp_map[gene_mask, :][:, bin_mask]

        subset_accessibility = accessibility[bin_mask, :]
        
        return subset_accessibility, subset_factor_binding, subset_rp_map

    def load_rp_matrix(self, data_object):

        self.log.append('Loading {} RP matrix ...'.format(self.technology))
        
        dataset_ids = data_object[self.config.get('accessibility_assay', 'reg_potential_dataset_ids').format(technology = self.technology)][...].astype(str)

        rp_matrix = data_object[self.config.get('accessibility_assay', 'reg_potential_matrix').format(technology = self.technology)][...]
        
        self.loaded = True
        return rp_matrix, dataset_ids

    @staticmethod
    def make_rp_matrix(*, profiles, rp_map):

        rp_matrix = rp_map.dot(profiles)

        return rp_matrix

    def get_delta_RP_p_value(self, gene_TF_scores, label_vector):
        '''
        gene_TF_scores: gene x TF, model output of delta-RP matrix. more purturbation of genes of interest correspond with higher delta regulation score
        '''
        #seperate matrix into query and background sets
        query_delta = gene_TF_scores[label_vector.astype(np.bool)]
        background_delta = gene_TF_scores[~label_vector.astype(np.bool)]

        #for each TF, calculate p-value of difference in delta-R distributions between query and background genes
        test_parameters = list(zip(query_delta.T, background_delta.T))

        try:
            self.cores
        except AttributeError:
            p_vals = [
                mannu_test_function((q,b)) for q,b in test_parameters
            ]
        else:
            with Pool(self.cores) as p:
                p_vals = p.map(mannu_test_function, test_parameters)

        _, p_values = list(zip(*p_vals))

        return p_values

    def get_info(self):
        raise NotImplementedError()

    def predict(self, gene_mask, label_vector, *args, data_object = None, debug = False, **kwargs):
        raise NotImplementedError()