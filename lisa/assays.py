
from lisa.utils import LoadingBar
from multiprocessing import Pool
import numpy as np
from scipy import stats, sparse
import os

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


class LISA_RP_Assay:

    invitro_metadata = ['factor','cell_line','cell_type','tissue']
    insilico_metadata = ['factor','dbd','description','cell_line']

    def __init__(self, *, technology, config, cores, log, metadata_headers, metadata_path, oneshot):
        self.config = config
        self.technology = technology
        self.cores = cores
        self.log = log
        self.loaded = False
        self.oneshot = False
        self.metadata_headers = metadata_headers
        self.metadata_path = metadata_path

    def load_rp_matrix(self, data_object, gene_mask = None):

        self.log.append('Loading {} RP matrix ...'.format(self.technology))
        
        dataset_ids = data_object[self.config.get('accessibility_assay', 'reg_potential_dataset_ids').format(technology = self.technology)][...].astype(str)

        if self.oneshot and not gene_mask is None:
            rp_matrix = data_object[self.config.get('accessibility_assay', 'reg_potential_matrix').format(technology = self.technology)][gene_mask, :]
        else:
            self.oneshot = False
            rp_matrix = data_object[self.config.get('accessibility_assay', 'reg_potential_matrix').format(technology = self.technology)][...]
        
        self.loaded = True
        return rp_matrix, dataset_ids


    def get_delta_RP_p_value(self, gene_TF_scores, label_vector):
        '''
        gene_TF_scores: gene x TF, model output of delta-RP matrix. more purturbation of genes of interest correspond with higher delta regulation score
        '''
        #seperate matrix into query and background sets
        query_delta = gene_TF_scores[label_vector.astype(np.bool)]
        background_delta = gene_TF_scores[~label_vector.astype(np.bool)]

        #for each TF, calculate p-value of difference in delta-R distributions between query and background genes
        test_parameters = list(zip(query_delta.T, background_delta.T))

        if self.cores == 1:
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

    def predict(self, gene_mask, label_vector, data_object = None, debug = False):
        raise NotImplementedError()

    def load_metadata(self):
        #reformat metadata tsv into conventiant id-indexed dictionary
        with open(self.metadata_path, 'r', encoding = 'latin') as metdata_file:
            metadata = [[field.strip() for field in line.split('\t')] for line in metdata_file.readlines()]
                                
        meta_headers, metadata = metadata[0], metadata[1:]

        return {metaline[0] : dict(zip(meta_headers[1:], metaline[1:])) for metaline in metadata}

    def load(self, data_object, gene_mask = None):
        self.rp_matrix, self.dataset_ids = self.load_rp_matrix(data_object, gene_mask = gene_mask)
        self.metadata = self.load_metadata()

    def get_metadata(self, sample_ids):
        return dict(
            sample_id = sample_ids,
            **{
                header : [self.metadata[_id][header] for _id in sample_ids]
                for header in self.metadata_headers 
            }
        )

class PeakRP_Assay(LISA_RP_Assay):

    def predict(self, gene_mask, label_vector, data_object, debug = False):

        try:
            self.rp_matrix
        except AttributeError:
            self.load(data_object, gene_mask=gene_mask)
            
        self.log.append('Calculating {} peak-RP p-values ...'.format(self.technology))
        
        #calculate p-values by directly applying mannu-test on RP matrix. Subset the RP matrix for genes-of-interest if required
        p_vals = self.get_delta_RP_p_value(self.rp_matrix[gene_mask, :] if not self.oneshot else self.rp_matrix, label_vector)

        return p_vals

    def get_info(self):
        return dict()

    def get_metadata(self):
        return super().get_metadata(self.dataset_ids)

#Generator that repeatedly yields the same factor_binding and rp_map matrices with a new accessibility profile 
#to the multiprocessing pool creator. This reduces data redundancy.
class KnockoutGenerator:

    def __init__(self, accessibility_profiles, factor_binding, rp_map):
        self.accessibility_profiles = accessibility_profiles
        self.rp_map = rp_map
        self.factor_binding = factor_binding

    def __iter__(self):
        for profile in self.accessibility_profiles.T:
            yield profile, self.factor_binding, self.rp_map
    
class Accesibility_Assay(LISA_RP_Assay):

    def __init__(self, *, technology, config, cores, log, metadata_path, oneshot, rp_map, factor_binding, selection_model, chromatin_model,
        metadata_headers = LISA_RP_Assay.invitro_metadata):
        super().__init__(technology = technology, config=config,cores=cores, log=log, metadata_headers=metadata_headers,
            metadata_path=metadata_path, oneshot=oneshot)
        self.rp_map = rp_map
        self.factor_binding = factor_binding
        self.selection_model = selection_model
        self.chromatin_model = chromatin_model

    def get_info(self):
        return dict(
            selection_model = self.selection_model.get_info(),
            chromatin_model = self.chromatin_model.get_info(),
            selected_datasets = self.get_metadata(list(self.selected_dataset_ids)),
        )

    def load_accessibility_profiles(self, data_object, selected_dataset_ids):

        loadingbar = LoadingBar('Reading {} data'.format(self.technology), len(selected_dataset_ids), 20)

        accessibility_profiles = []
        for selected_dataset in selected_dataset_ids:
            self.log.append(loadingbar, update_line = True)
            accessibility_profiles.append(
                data_object[self.config.get('accessibility_assay','binned_reads')\
                .format(technology = self.technology, dataset_id = selected_dataset)][...][:,np.newaxis]
            )

        accessibility_profiles = np.concatenate(accessibility_profiles, axis = 1)

        return accessibility_profiles

    #iterator for distributing these matrices to multiple cores
    

    def calculate_ISDs(self, accessibility_profiles, factor_binding, rp_map): #enforced kwargs
        """
        log: log object for printing
        accessibility_profiles: (bins, num_datasets): list of chromatin-accessiblity datasets on which to analyze the effects of ISD
        factor_binding: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
        rp_map: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)

        returns:
        delta_regulatory_score (genes x TFs)
        """
        '''
        # Code for non-multiprocessing implementation of this ISD algorithm. Much slower.
        knockouts = []
        loadingbar = LoadingBar('Performing in-silico knockouts', accessibility_profiles.shape[1], 20, cold_start = True)
        log.append(loadingbar, update_line = True)
        for profile in accessibility_profiles.T:
            knockouts.append(get_delta_RP(profile, factor_binding, rp_map))
            log.append(loadingbar, update_line = True)'''

        self.log.append('Performing in-silico knockouts ...')

        with Pool(self.cores) as p:
            knockouts = p.map(delta_RP_wrapper, iter(KnockoutGenerator(accessibility_profiles, factor_binding, rp_map)))
        
        #concatenate datasets to for gene x TF x datasets shaped matrix
        self.log.append('Calculating Δ regulatory score ...')
        datacube = np.concatenate(knockouts, axis = 1)

        num_genes, _, num_TFs = datacube.shape

        delta_regulation_score = self.chromatin_model.get_deltaRP_activation(datacube)
        
        assert(delta_regulation_score.shape == (num_genes, num_TFs))

        return delta_regulation_score, datacube

    
    def predict(self, gene_mask, label_vector, data_object, debug = False):

        with self.log.section('Modeling {} purturbations:'.format(self.technology)):

            try:
                self.rp_matrix
            except AttributeError:
                self.load(data_object, gene_mask= gene_mask)

            #find bins of gene-subsetted rp-map with RP > 0
            bin_mask = np.squeeze(np.array(self.rp_map[gene_mask, : ].tocsc().sum(axis = 0) > 0))
            #subset rp_map and factor hits on bins with RP > 0
            subset_factor_binding = self.factor_binding[bin_mask, :]
            
            subset_rp_map = self.rp_map[gene_mask, :][:, bin_mask]

            subset_rp_matrix = self.rp_matrix[gene_mask, :] if not self.oneshot else self.rp_matrix

            #DNase model building and purturbation
            self.log.append('Selecting discriminative datasets and training chromatin model ...')

            #select the most discriminative datasets
            dataset_mask = self.selection_model.fit(subset_rp_matrix, label_vector)
            #subset the best datasets
            subset_rp_matrix, self.selected_dataset_ids = subset_rp_matrix[:, dataset_mask], self.dataset_ids[dataset_mask]
            #fit the chromatin model to these datasets
            self.chromatin_model.fit(subset_rp_matrix, label_vector, self.cores)
                    
            with self.log.section('Calculating in-silico deletions:'):

                accesibility_profiles = self.load_accessibility_profiles(data_object, self.selected_dataset_ids)

                accesibility_profiles = accesibility_profiles[bin_mask, :]
                
                delta_reg_scores, datacube = self.calculate_ISDs(accesibility_profiles, subset_factor_binding, subset_rp_map)

                self.log.append('Calculating p-values ...')
                
                p_vals = self.get_delta_RP_p_value(delta_reg_scores, label_vector)
        
            self.log.append('Done!')

        if debug:
            return p_vals, dict(
                gene_mask = gene_mask,
                label_vector = label_vector,
                subset_rp_matrix = subset_rp_matrix,
                subset_rp_map = subset_rp_map,
                subset_factor_binding = subset_factor_binding,
                datacube = datacube,
                delta_reg_scores = delta_reg_scores,
                dataset_mask = dataset_mask,
                bin_mask = bin_mask,
            )
        else:
            return p_vals