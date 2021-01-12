
from lisa.core.io import parse_bedfile
from lisa.core.data_interface import DataInterface
from lisa.core.genome_tools import Region, BadRegionError
import argparse

def main(species, window_size, path):

    region_fields = parse_bedfile(path, header = False)

    regions = [Region(*r) for r in region_fields]

    genome = DataInterface.load_genome(species, window_size)

    indices = []

    for region in regions:
        try:
            windows = genome.get_region_windows(region)
            indices.extend(windows)
        except BadRegionError:
            pass

    return list(set(indices))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Using summit file or motif hits file, find genome indices where regions overlap')
    parser.add_argument('species', type = str, choices = ['hg38','mm10'])
    parser.add_argument('window_size', type = int)
    parser.add_argument('input', type = str)
    parser.add_argument('output', type = str)

    args = parser.parse_args()

    indices = main(args.species, int(args.window_size), args.input)

    with open(args.output, 'w') as f:
        print('\n'.join(map(str, indices)), file = f)