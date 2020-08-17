import configparser
import h5py as h5
import numpy as np
import pandas as pd
from utils import LoadingBar
from collections import defaultdict

'''Steps to reformatting a lisa_v1 to lisa_v2 species-specific data repo:
1. ~~Get gene lists from tad dataframe, h3k27ac data, and Dnase data~~
2. ~~Intersect gene lists to make master list, subset the TAD data and define gene_info folder~~

For each accessibility assay:
3. ~~Get QC from metadata, and use qc to subset samples in Dnase and acetylation data. Also subset by used genes~~
4. ~~Reformat binned reads array~~

For each TF-binding assay:
5. ~~Aggregate by factor, save as factor-indexed arrays of binned hits~~'''


def subset_accessibility_rp_data(path_to_h5, gene_set, qc_status):
    
    with h5.File(path_to_h5, 'r') as acc_data:
        
        symbols = [line.decode('utf-8').split(':')[-1] for line in acc_data['RefSeq'][...]]

        symbol_dict = {}
        for i, symbol in enumerate(symbols):
            if not symbol in symbol_dict:
                symbol_dict[symbol] = i

        gene_subset = np.zeros(len(symbols), dtype = np.bool)
        gene_subset[list(symbol_dict.values())] = True
        
        dataset_ids = acc_data['IDs'][...].astype(str)
        
        qc_subset = qc_status[dataset_ids].values.astype(np.bool)
        
        rp_data = acc_data['RP'][:, qc_subset][gene_subset, :]
        
    return rp_data, np.array(symbols)[gene_subset], dataset_ids[qc_subset]

def reformat_accessibility_assay_RP(h5_object, technology, config, path_to_h5, gene_set, qc_status):

    rp_data, genes, sid = subset_accessibility_rp_data(path_to_h5, gene_set, qc_status)
    
    h5_object.create_dataset(
        config['accessibility_assay']['reg_potential_gene_symbols'].format(technology = technology),
        data = np.array(genes).astype('S')
    )
    
    h5_object.create_dataset(
        config['accessibility_assay']['reg_potential_dataset_ids'].format(technology = technology),
        data = np.array(sid).astype('S')
    )
    
    h5_object.create_dataset(
        config['accessibility_assay']['reg_potential_matrix'].format(technology = technology),
        data = rp_data
    )

    
def reformat_TF_hits(h5_object, config, technology, path_to_h5, window_converter, metadata):
    
    with h5.File(path_to_h5, 'r') as tf_hits:

        num_samples = len(list(tf_hits.keys()))
        
        loading_bar = LoadingBar('\tCollecting binding data', num_samples, 20)
        
        for sample in tf_hits.keys():
            print(loading_bar, end = '')
            
            if not sample == 'IDs':
                
                try:
                    sample_metadata = metadata.loc[str(sample)]
                
                    if sample_metadata.qc == 1:
                        
                        factor_name = sample_metadata.factor.replace('/', '-').replace('_','-')

                        peaks = tf_hits[sample][...]
                        peaks = window_converter[peaks]
                        
                        h5_object.create_dataset(
                            config['factor_binding']['tf_binding_data']\
                                .format(
                                    technology = technology, dataset_id = "{}_{}".format(factor_name, sample)
                                ),
                            data = np.array(list(peaks)).astype(np.int64)
                        )
                        
                except OSError:
                    print('\n\tError saving data for sample {}, factor: {}'\
                          .format(str(sample), sample_metadata.factor))
                except KeyError:
                    print('\n\tError: No metadata for sample {}'.format(str(sample)))
                    

def index_binned_reads(h5_object, technology, config, path_to_h5, qc_status):
    
    print('\tReading binned reads data ...')

    with h5.File(path_to_h5, 'r') as binned_reads:

        dataset_ids = binned_reads['IDs'][...].astype(str)

        dataset_subset = qc_status[dataset_ids].astype(np.bool).values

        reads_matrix = binned_reads['OrderCount'][:, dataset_subset]

    loading_bar = LoadingBar('\tWriting ID-indexed data', dataset_subset.sum(), 20)

    for i, dataset_id in enumerate(dataset_ids[dataset_subset]):
        print(loading_bar, end = '')
        h5_object.create_dataset(
            config['accessibility_assay']['binned_reads'].format(technology = technology, dataset_id = dataset_id),
            data = reads_matrix[:, i]
        )

def convert_lisa_data_formatting(new_dataset_name, 
    create_new = True, *
    metadata_file,
    tad_data,
    window_conversion_file,
    dnase_reads,
    acetylation_reads,
    chip_peaks,
    motif_peaks,
    dnase_rp,
    acetylation_rp,
    config_file,
    motif_metadata,
): 

    config = configparser.ConfigParser()
    config.read(config_file)

    print('Reading metadata ...')

    tad_data = pd.read_csv(tad_data, encoding = "ISO-8859-1", sep = '\t')

    with h5.File(dnase_rp, 'r') as dnase_rp_data:
        dnase_gene_symbols = [line.split(":")[-1] for line in dnase_rp_data['RefSeq'][...].astype(str)]
        
    with h5.File(acetylation_rp, 'r') as acetylation_rp_data:
        acetyalation_gene_symbols = [line.split(":")[-1] for line in acetylation_rp_data['RefSeq'][...].astype(str)]

    print('Intersecting gene lists. Compiling gene TAD information ...')
    combined_sybmols = set(tad_data.geneName) & set(dnase_gene_symbols) & set(acetyalation_gene_symbols)    

    tad_data = tad_data[tad_data.geneName.isin(combined_sybmols)]

    gene_symbols = tad_data.geneName.values
    gene_ids = tad_data.geneID.str.split('.').str.get(0).values
    tad_domain = (tad_data.k4me3_order_cluster.astype(str) + ',' + tad_data.tad_order_cluster.astype(str)).values

    metadata = pd.read_csv(metadata_file, sep = '\t', encoding = "ISO-8859-1").set_index('id')
    metadata.index = metadata.index.astype(np.str)

    window_converter = np.load(window_conversion_file)

    motif_metadata = pd.read_csv(motif_metadata, sep = '\t').set_index('id')
    motif_metadata['qc'] = 1
    motif_metadata['factor'] = motif_metadata.index

    with h5.File(new_dataset_name, 'w' if create_new else 'a') as data:
        
        
        data.create_dataset(config['gene_info']['gene_symbols'], data = gene_symbols.astype('S'))
        data.create_dataset(config['gene_info']['gene_refseqids'], data = gene_ids.astype('S'))
        data.create_dataset(config['gene_info']['tad_domains'], data = tad_domain.astype('S'))
        
        print('DNase:')
        print('\tReformatting DNase regulatory potential matrix ...')
        reformat_accessibility_assay_RP(data, 'DNase', config, dnase_rp, combined_sybmols, 
                                    metadata.qc)
        
        index_binned_reads(data, 'DNase', config, dnase_reads, metadata.qc)
        
        print('H3K27ac:')
        print('\tReformatting h3k27ac regulatory potential matrix ...')
        reformat_accessibility_assay_RP(data, 'H3K27ac', config, acetylation_rp,
                                    combined_sybmols, metadata.qc)
        
        index_binned_reads(data, 'H3K27ac', config, acetylation_reads, metadata.qc)
        
        print('Motifs:')    
        reformat_TF_hits(data, config, 'motif_hits_1000kb', motif_peaks, window_converter, motif_metadata)
        
        print('Chip-seq:')
        reformat_TF_hits(data, config, 'ChIP-seq_1000kb', chip_peaks, window_converter, metadata)