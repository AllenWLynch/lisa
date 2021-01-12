
import argparse
from lisa.core.data_interface import DataInterface
import numpy as np
import pandas as pd
import os

def main(args):

    cistrome_metadata = pd.read_csv(args.cistrome_metadata, sep = '\t').set_index('DCid')
    cistrome_metadata.index = cistrome_metadata.index.astype(str)

    data = DataInterface(args.species, window_size= args.window_size, download_if_not_exists=False,
        make_new=False, load_genes=True)

    rp_map_styles = data.get_rp_maps()

    if len(rp_map_styles) == 0:

        basic_rp_map, enhanced_rp_map = data.build_binned_rp_map('basic',10000), data.build_binned_rp_map('enhanced', 10000)

        data.add_rp_map('basic_10K', basic_rp_map)
        data.add_rp_map('enhanced_10K', enhanced_rp_map)

    else:

        basic_rp_map = data.get_rp_map('basic_10K')
        enhanced_rp_map = data.get_rp_map('enhanced_10K')


    for arr_name in args.coverage_arrays:
        
        coverage_array = np.load(arr_name)
        
        technology, dataset_id = os.path.basename(arr_name).split('_')

        dataset_id = '.'.join(dataset_id.split('.')[:-1])

        metadata_headers = data.get_metadata_headers(technology)

        meta_dict = cistrome_metadata.loc[dataset_id, metadata_headers].to_dict()

        data.add_profile_data(technology, dataset_id, coverage_array, [basic_rp_map, enhanced_rp_map],
            ['basic_10K','enhanced_10K'], **meta_dict)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('window_size', type = int)
    parser.add_argument('cistrome_metadata', type = str)
    parser.add_argument('coverage_arrays', type = str, nargs = '+')

    args = parser.parse_args()

    main(args)