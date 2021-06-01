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

def main(*,species, motif_bed, window_size, dataset_id, output):

    genome = DataInterface.load_genome(species, window_size)

    factor_name = None
    window_nums, scores = [],[]

    #adjust p-val cutoff based on the filesize (only affects p-val if file is huge)
    pval_cutoff = 430

    with open(output, 'w') as o:

        with gzip.open(motif_bed, 'rb') as bed:

            for i, line in enumerate(bed):
                
                chrom, start, end, factor, relscore, log_pval, strand = line.decode('utf-8').strip().split('\t')
                
                if i == 0:
                    factor_name = factor
                    print('Binning {} motifs with pval cutoff of {} ...'.format(factor_name.upper(), str(pval_cutoff)), file = stderr)

                neg_log10_pval = int(log_pval)

                if neg_log10_pval >= pval_cutoff:

                    try:
                        hit_windows = genome.get_region_windows(Region(chrom, start, end))
                        
                        for hit_window in hit_windows:
                            print(dataset_id, hit_window, neg_log10_pval, sep = '\t', file = o)                    

                    except BadRegionError:
                        pass
    
    print(dataset_id, factor_name.upper(), 'JASPAR', sep = '\t')

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('motif_bed', type = str)
    parser.add_argument('-w', '--window_size', type = int, required=True)
    parser.add_argument('-o', '--output', type = str, required=True)
    parser.add_argument('-i','--id', type = str, required=True)
    args = parser.parse_args()

    main(species = args.species, motif_bed = args.motif_bed, 
        window_size = args.window_size, dataset_id = args.id, output = args.output)