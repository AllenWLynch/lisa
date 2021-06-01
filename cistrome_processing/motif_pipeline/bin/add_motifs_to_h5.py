
from lisa.core.data_interface import DataInterface
import argparse
import os
import numpy as np
import pandas as pd
from scipy import sparse

TECHNOLOGY = 'Motifs'

def main(species, window_size, motif_metadata, bin_sorted_hits, group_loci = 100000):

    motif_metadata = pd.read_csv(motif_metadata, sep = '\t', header = None)
    motif_metadata.columns = ['dataset_id', 'factor', 'source']
    motif_metadata = motif_metadata.set_index('dataset_id')
    motif_metadata = motif_metadata.drop_duplicates()

    data = DataInterface(species, window_size= window_size, download_if_not_exists=False,
        make_new=False, load_genes=False)

    print(data.path)
    raise Exception()

    data.create_binding_dataset(TECHNOLOGY, motif_metadata.index.values, **motif_metadata.to_dict('list'))

    id_to_idx_map = dict(zip(data.list_binding_datasets(TECHNOLOGY), np.arange(len(data.list_binding_datasets(TECHNOLOGY)))))

    current_pos = 0
    last_added_chunk = 0
    i = 0
    rows,cols,scores=[],[],[]

    with open(bin_sorted_hits, 'r') as f:

        for line in f:
            motif_id, bin_num, score = line.strip().split()

            bin_num = int(bin_num)
            
            if bin_num < current_pos:
                raise Exception('Input file not sorted!')
            elif bin_num > current_pos and i >= group_loci:
                print('Adding matrix segment ...')
                matrix_form = sparse.coo_matrix((scores, (rows, cols))).tocsr()
                data.append_csr(TECHNOLOGY, matrix_form)
                last_added_chunk = bin_num
                i=0
                rows,cols,scores=[],[],[]

            tf_idx = id_to_idx_map[motif_id]
            rows.append(bin_num - last_added_chunk)
            cols.append(tf_idx)
            scores.append(int(score))
            current_pos = bin_num
            i+=1
        
        if len(rows) > 0:
            matrix_form = sparse.coo_matrix((scores, (rows, cols))).tocsr()
            data.append_csr(TECHNOLOGY, matrix_form)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Saves indices lists to factor binding h5. Filename specifies metadata: \{technology\}-\{dataset_id\}-\{metadata1_metadata2_...\}')
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('window_size', type = int)
    parser.add_argument('motif_metadata', type = str)
    parser.add_argument('hits', type = str)

    args = parser.parse_args()

    main(args.species, int(args.window_size), args.motif_metadata, args.hits)