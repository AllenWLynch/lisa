
from lisa import FromCoverage
import numpy as np
import argparse

def main(args):
    
    coverage_array = FromCoverage.convert_bigwig(args.bigwig_file, args.species)
    
    np.save(args.name, coverage_array)


if __name__ == "__main__":

    bw_parser = argparse.ArgumentParser()
    bw_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Pileup over hg38 or mm10 genome')
    bw_parser.add_argument('bigwig_file', type = str, help = 'bigwig file to convert into coverage array')
    bw_parser.add_argument('-n', '--name', type = str, required = True, help = 'Prefix under which to save coverage array file')

    args = bw_parser.parse_args()

    main(args)