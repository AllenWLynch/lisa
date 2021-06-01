
from lisa.core.data_interface import DataInterface
import argparse
import os
import numpy as np
import pandas as pd

def main(species, window_size, cistrome_metadata, motif_metadata, index_files):

    cistrome_metadata = pd.read_csv(cistrome_metadata, sep = '\t').set_index('DCid')
    cistrome_metadata.index = cistrome_metadata.index.astype(str)
    motif_metadata = pd.read_csv(motif_metadata, sep = '\t', header = None)
    motif_metadata.columns = ['dataset_id', 'factor', 'source']
    motif_metadata = motif_metadata.set_index('dataset_id')
    motif_metadata = motif_metadata.drop_duplicates()

    data = DataInterface(species, window_size= window_size, download_if_not_exists=False,
        make_new=False, load_genes=False)

    dataset_ids = motif_metadata.index.values

    data.create_binding_dataset('Motif', dataset_ids)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Saves indices lists to factor binding h5. Filename specifies metadata: \{technology\}-\{dataset_id\}-\{metadata1_metadata2_...\}')
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('window_size', type = int)
    parser.add_argument('cistrome_metadata', type = str)
    parser.add_argument('motif_metadata', type = str)
    parser.add_argument('index_files', type = str, nargs='+')

    args = parser.parse_args()

    main(args.species, int(args.window_size), args.cistrome_metadata, args.motif_metadata, args.index_files)