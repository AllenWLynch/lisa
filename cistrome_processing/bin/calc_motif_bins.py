
from urllib import request
import gzip
from lisa.core.data_interface import DataInterface
from lisa.core.genome_tools import Region, BadRegionError
from lisa.core.utils import Log, LoadingBar
import argparse
from scipy import sparse
import os
from scipy.stats import gamma
import numpy as np
from sys import stdout, stderr

def main(species, motif_bed, window_size, gamma_threshold = 0.95):

    genome = DataInterface.load_genome(species, window_size)

    log = Log(target= stderr)
    
    factor_name = None
    window_nums, scores = [],[]

    with gzip.open(motif_bed, 'rb') as f:

        bed = f.readlines()

        bar = LoadingBar('Binning {} motif hits'.format(str(len(bed))), len(bed), cold_start=True)

        for i, line in enumerate(bed):
            
            chrom, start, end, factor, relscore, log_pval, strand = line.decode('utf-8').strip().split('\t')
            
            if i == 0:
                factor_name = factor

            try:
                hit_windows = genome.get_region_windows(Region(chrom, start, end))
                window_nums.extend(hit_windows)

                scores.extend([float(log_pval)/100]*len(hit_windows))

            except BadRegionError:
                pass

            log.append(bar, update_line=True)
    
    log.append('')
        
    log.append('Done')

    hits = sparse.csc_matrix((scores, window_nums, [0, len(window_nums)]), shape = (len(genome), 1)).tocoo().tocsc()

    sample_hit_scores = np.random.choice(np.array(hits.todense()).reshape(-1), size = 10000)

    min_bin_score = gamma(*gamma.fit(sample_hit_scores)).ppf(gamma_threshold)

    hit_indices = hits.indices[(hits.data >= min_bin_score) & (hits.data > 0)]
    
    return hit_indices, factor_name


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('motif_bed', type = str)
    parser.add_argument('-w', '--window_size', type = int, required=True)
    parser.add_argument('-o', '--output', type = str, required=True)
    args = parser.parse_args()

    hit_indices, factor_name = main(args.species, args.motif_bed, args.window_size)

    motif_id = '.'.join(os.path.basename(args.motif_bed).split('.')[:2])

    print(motif_id, factor_name, 'JASPAR', sep = '\t', file = stdout)

    with open(args.output, 'w') as f:
        print('\n'.join(map(str, hit_indices)), file = f)    