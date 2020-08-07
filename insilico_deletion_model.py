
import numpy as np
from scipy import sparse
from chromatin_model import ChromatinModel

#investigate density of everything
def calculate_dataset_deltaRP(reads_per_bin, ChIP_per_bin, rp_matrix):

    #validate sparse types
    assert( isinstance(reads_per_bin, sparse.csc_matrix) )
    assert( isinstance(ChIP_per_bin, sparse.csc_matrix) )
    assert( isinstance(rp_matrix, sparse.csr_matrix) )

    #validate shapes
    num_bins = rp_matrix.shape[0]
    assert( reads_per_bin.shape == (num_bins, 1)) # bins x 1
    assert( ChIP_per_bin.shape[0] == num_bins ) # bins x TFs

    # reads_per_bin (*) ChIP_per_bin = reads_per_TF (bins x TFs)
    # reads_per_TF^T \dot rp_matrix = deleted RP (TF x gene)
    # deleted RP^T (*) dataset_coeffcient = weighted delta RP (gene x TFs)
    # weighted delta RP + running_regression_total (just add this to a total and crunch one dataset at a time rather than build huge 3d datacube)
    delta_rp = reads_per_bin.multiply(ChIP_per_bin).transpose().dot(rp_matrix).transpose().todense()[:, :, np.newaxis]

    return delta_rp

def get_delta_regulation(datasets, ChIP_per_bin, rp_matrix, chromatin_model):

    assert( isinstance(chromatin_model, ChromatinModel))

    datacube = np.concatenate(
        [calculate_dataset_deltaRP(dataset, ChIP_per_bin, rp_matrix) # genes x TFs
            for dataset in datasets],
        axis = 2
    )

    num_genes, num_TFs, num_samples = datacube.shape

    isd_conditions = datacube.reshape((-1, num_samples)) # (genes + TFs) x samples

    delta_regulation_score = chromatin_model.get_delta_regulation(isd_conditions).reshape(num_genes, num_TFs) # genes x TFs
    
    return delta_regulation_score

def get_delta_regulation_significance(delta_regulation_score, labels):
    pass





