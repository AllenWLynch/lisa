
import numpy as np
from scipy import sparse, stats
from chromatin_model import ChromatinModel


def calculate_dataset_deltaRP(reads_per_bin, ChIP_per_bin, rp_matrix):
    """
    For a given chromatin accessibility dataset, calculated the change in RP from knocking out each TF
    reads_per_bin: (num_bins, ), an array describing the reads in each 1kb bin of the genome for a chrom accessibility assay
    ChIP_per_bin: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
    rp_matrix: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)
    """
    #validate sparse types
    assert( isinstance(reads_per_bin, np.array )
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
    """
    Loops through datasets-of-interest, performing ISD for each. Then, compiles the results and gets the effective change with respect to the chromatin model trained in the previous step
    datasets: (bins, num_datasets): list of chromatin-accessiblity datasets on which to analyze the effects of ISD
    ChIP_per_bin: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
    rp_matrix: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)
    chromatin_model: a ChromatinModel object that implements fit and get_deltaRP_activation. Purturbation of this model summarizes effect of ISD

    returns: delta_regulation_score: (genes x TFs)
    """

    assert( isinstance(chromatin_model, ChromatinModel) )

    #potentially multi-thread this
    datacube = np.concatenate(
        [calculate_dataset_deltaRP(dataset, ChIP_per_bin, rp_matrix) # genes x TFs
            for dataset in datasets],
        axis = 2
    )

    num_genes, num_TFs, num_samples = datacube.shape

    isd_conditions = datacube.reshape((-1, num_samples)) # (genes + TFs) x samples

    delta_regulation_score = chromatin_model.get_delta_regulation(isd_conditions).reshape(num_genes, num_TFs) # genes x TFs
    
    return delta_regulation_score

def get_delta_regulation_significance(tf_delta_regulation_score, labels, alternative = 'greater'):
    """
    Performs wilcoxon test on the results of one TF's change in regulatory score on query vs. background genes
    tf_delta_regulation_score: (num genes, ): 0D array for one TF's effecst on reguatory score
    labels: query vs. background labels
    alternative: alternative hypothesis type. Usually, testing to see if delta regulatory score is greater for query genes

    returns:
    p-value of wilcoxon test
    """
    query_delta_rp = tf_delta_regulation_score[labels.astype(np.bool)]
    background_delta_rp = tf_delta_regulation_score[ ~labels.astype(np.bool) ]

    return stats.wilcoxon(query_delta_rp, background_delta_rp, alternative = alternative)[1]


def get_insilico_deletion_significance(*, datasets, ChIP_per_bin, rp_matrix, chromatin_model, labels, alternative = 'greater'): #enforced kwargs
    """
    datasets: (bins, num_datasets): list of chromatin-accessiblity datasets on which to analyze the effects of ISD
    ChIP_per_bin: (num_bins, TFs), a sparse binary map showing bins that contain a chip-seq or motif peak to knock out
    rp_matrix: (num_bins, genes), precalculated matrix mapping the reads in a bin to a regulatory score for each gene (this is a huge matrix)
    chromatin_model: a ChromatinModel object that implements fit and get_deltaRP_activation. Purturbation of this model summarizes effect of ISD
    labels: query vs. background labels
    alternative: alternative hypothesis type. Usually, testing to see if delta regulatory score is greater for query genes

    returns:
    p-value of wilcoxon test for each TF
    """
    #get the delta reg score for each TF on each gene, summarized for each dataset
    delta_regulation_score = get_delta_regulation(datasets, ChIP_per_bin, rp_matrix, chromatin_model) # (genes, TFs)

    #apply the wilcoxon test on each TF
    tf_p_vals = np.apply_along_axis(get_delta_regulation_significance, 0, delta_regulation_score) # (TFs, )

    return tf_p_vals

