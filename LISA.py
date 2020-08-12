
import build_chromatin_model
import background_genes_selection
import insilico_deletion_model as ISD

from models import LR_BinarySearch_SampleSelectionModel
from models import LR_ChromatinModel

if __name__ == "__main__":

#load gene lists
with h5.File('data/hg38/lisa_v2_hg38.h5', 'r') as store:
    gene_symbols = store['gene_info']['gene_symbols'][...].astype(str)
    gene_ids = store['gene_info']['gene_refseqIDs'][...].astype(str)
    gene_tag_groups = store['gene_info']['gene_tad_info'][...].astype(str)

#open user-supplied gene-list
with open('data/gene_list.txt', 'r') as f:
    user_gene_list = f.readlines()

gene_list = [gene.strip().upper() for gene in user_gene_list]

#select background genes from TAD distributions
selected_genes, labels = background_genes_selection\
.select_genes_for_chromatin_model(gene_list, gene_symbols, gene_ids, gene_tag_groups, 
                                  num_selected = 3000, user_background_genes = None, method = 'TAD')

#Load gene x dataset RP matrix
with h5.File('data/hg38/lisa_v2_hg38.h5', 'r') as data:
    
    gene_symbols = data['DNase_data']['gene_ids'][...].astype(str)
    intersected_genes = set(selected_genes).intersection(set(gene_symbols))
    intersected_ids = np.isin(gene_symbols, list(intersected_genes))
    
    rp_matrix = data['DNase_data']['regulatory_potential'][intersected_ids, :][...]
    dataset_ids = data['DNase_data']['sample_ids'][...].astype(str)

labels = labels[np.isin(selected_genes, list(intersected_genes))]

#instantiate models
sample_selection_model = LR_BinarySearch_SampleSelectionModel()

chromatin_model = LR_ChromatinModel({'C' : list(10.0**np.arange(-4,4,1))})

#train models
selected_datasets, selected_dataset_ids, selection_model, chromatin_model, normalization_fn\
    = build_chromatin_model.build_chromatin_model(rp_matrix, dataset_ids, 
                                            labels, sample_selection_model, chromatin_model, n_jobs = -1)

#read in chipseq data and build sparse matrix
with h5.File('data/hg38/lisa_v2_hg38.h5', 'r') as data:

    tf_ids = {tf : i for i, tf in enumerate(data['ChIP_data'].keys())}

    peaks = np.concatenate([
        data['ChIP_data'][factor][...].astype(np.int32)
        for factor in tf_ids.keys()
    ])

    data['ChIP_data']['SMAD4'].shape

    tf_index = np.concatenate([
        np.full(data['ChIP_data'][factor].shape, factor_id)
        for factor, factor_id in tf_ids.items()    ])

#read in DNase samples, filter by important to ChIP model
with h5.File('data/hg38/lisa_v2_hg38.h5', 'r') as data:
    
    num_peaks = peaks.shape
    num_bins = data['DNase_data']['reads_data'][selected_dataset_ids[0]].shape[0]

    chip_sparse = sparse.coo_matrix(
        (
            np.ones(num_peaks), 
            (peaks, tf_index)
        ), 
        shape = (num_bins, len(tf_ids))
    )

    chip_filter = np.array(chip_sparse.sum(axis = 1) >= 1).reshape(-1)

    chip_sparse = sparse.csr_matrix(chip_sparse)[np.array(chip_filter).reshape(-1), :]

    dnnase_reads = np.concatenate([
        data['DNase_data']['reads_data'][dataset_id][chip_filter][...].reshape((-1,1))
        for dataset_id in selected_dataset_ids
    ], axis = 1)

#Import RP matrix
rp_shape = (dnnase_reads.shape[0], labels.shape[0])
num_datapoints = 1e5

random_gene = np.random.choice(rp_shape[1], int(num_datapoints))
random_bin = np.random.choice(rp_shape[0], int(num_datapoints))

rp_matrix = sparse.coo_matrix(
    (
        np.ones(int(num_datapoints)),
        (random_bin, random_gene)
    ),
    shape = rp_shape
)

ISD.get_insilico_deletion_significance(datasets = dnnase_reads, ChIP_per_bin=chip_sparse, 
                rp_matrix= rp_matrix, chromatin_model = chromatin_model, labels = labels)
