from lisa.core.utils import indices_list_to_sparse_array
from lisa.core.utils import LoadingBar
import fire
import h5py as h5
from scipy import sparse

def reformat_TF_hits(*,h5_object, save_name, factor_ids, offset, num_bins):

    with open(factor_ids, 'r') as f:
        factor_ids = [x.strip() for x in f.readlines()]
    
    with h5.File(h5_object, 'r') as tf_hits:

        num_samples = len(factor_ids)
        
        loading_bar = LoadingBar('\tCollecting binding data', num_samples, 20)
        
        peaks_list = []
        for sample in factor_ids:
            print('\r',loading_bar, end = '')
            #try:
            peaks = tf_hits[sample][...] - offset
            
            peaks_list.append(peaks)
                    
            #except OSError:
            #    print('\n\tError saving data for sample {}, factor: {}'\
            #            .format(str(sample), sample_metadata.factor))
            #except KeyError:
            #    print('\n\tError: No metadata for sample {}'.format(str(sample)))
                    
    tf_hits = indices_list_to_sparse_array(peaks_list, int(num_bins))
        
    sparse.save_npz(save_name, tf_hits)

if __name__ == "__main__":
    fire.Fire(reformat_TF_hits)