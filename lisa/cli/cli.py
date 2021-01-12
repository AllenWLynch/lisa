'''
********
Lisa CLI
********

Installing LISA using pip or conda adds the "lisa" command to your path. LISA's functionality is divided into three main subcommands:

* `lisa oneshot`_ : one genelist
* `lisa multi`_ : multiple genelists
* `lisa regions`_ : one genelist and a list of regions

Which are used depending on the evidence you have on hand. 

See the `User Guide <user_guide.rst>`_ for more usage information.
See the `Python API <python_api.rst>`_ for more in-depth description of tests and parameters.

'''

from lisa import FromRegions, FromGenes, FromCoverage
from lisa.core.utils import Log
from lisa.core.lisa_core import DownloadRequiredError
from lisa.core.data_interface import DatasetNotFoundError
from lisa._version import __version__
import configparser
import argparse
import os
import sys
import json
from collections import defaultdict
from shutil import copyfile
import lisa.cli.test_cli as tests
from shutil import copyfile
import numpy as np

from lisa.lisa_public_data.genes_test import _config as public_config
from lisa.lisa_user_data.regions_test import _config as user_config

#____COMMAND LINE INTERFACE________

INSTANTIATION_KWARGS = ['isd_method','verbose','assays', 'rp_map']
PREDICTION_KWARGS = ['background_list','num_background_genes','background_strategy', 'seed']

def extract_kwargs(args, keywords):
    return {key : vars(args)[key] for key in keywords}

def is_valid_prefix(prefix):

    if '/' in prefix:
        if os.path.isdir(prefix) or os.path.isfile(prefix) or os.path.isdir(os.path.dirname(prefix)):
            return prefix
        else:
            raise argparse.ArgumentTypeError('{}: Invalid file prefix.'.format(prefix))
    else:
        return prefix

def save_results(args, results, metadata):

    if args.save_metadata:
        if args.output_prefix:
            metadata_filename = args.output_prefix + '.metadata.json' 
        else:
            metadata_filename = os.path.basename(args.query_list.name) + '.metadata.json'

        with open(metadata_filename, 'w') as f:
            f.write(json.dumps(metadata, indent=4))

    if not args.output_prefix is None:
        with open(args.output_prefix + '.lisa.tsv', 'w') as f:
            f.write(results.to_tsv())
    else:
        print(results.to_tsv())

def lisa_oneshot(args):

    try:
        args.background_list = args.background_list.readlines()
    except AttributeError:
        pass

    results, metadata = FromGenes(args.species, **extract_kwargs(args, INSTANTIATION_KWARGS)).predict(args.query_list.readlines(), **extract_kwargs(args, PREDICTION_KWARGS))
    
    save_results(args, results, metadata)


def lisa_regions(args):

    try:
        args.background_list = args.background_list.readlines()
    except AttributeError:
        pass

    if not args.macs_xls:
        fn = FromRegions.using_bedfile
    else:
        fn = FromRegions.using_macs_output

    results, metadata = fn(args.species, args.query_genes, args.regions, rp_map = args.rp_map,
        rp_decay=args.rp_decay, isd_method=args.isd_method, background_list=args.background_list,
        background_strategy=args.background_strategy, num_background_genes = args.num_background_genes,
        seed=args.seed, header = args.header)
    
    save_results(args, results, metadata)

def lisa_coverage(args):

    try:
        args.background_list = args.background_list.readlines()
    except AttributeError:
        pass

    results, metadata = FromCoverage.using_bigwig(args.species, args.query_genes, args.bigwig_path, rp_map = args.rp_map,
        isd_method=args.isd_method, background_list=args.background_list,
        background_strategy=args.background_strategy, num_background_genes = args.num_background_genes,
        seed=args.seed)
    
    save_results(args, results, metadata)


def save_and_get_top_TFs(args, query_name, results, metadata):

    with open(args.output_prefix + query_name + '.lisa.tsv', 'w') as f:
        f.write(results.to_tsv())

    if args.save_metadata:
        with open(args.output_prefix + query_name + '.metadata.json', 'w') as f:
            f.write(json.dumps(metadata, indent=4))

    top_TFs = results.to_dict()['factor']

    return list(set(top_TFs[:10]))

def print_results_multi(results_summary):
    print('Sample\tTop Regulatory Factors:')
    for result_line in results_summary:
        print(result_line[0], ', '.join(result_line[1]), sep = '\t')

def lisa_multi(args):

    log = Log(target = sys.stderr, verbose = args.verbose)
    lisa = FromGenes(args.species, **extract_kwargs(args, INSTANTIATION_KWARGS), log = log)

    query_dict = {os.path.basename(query.name) : query.readlines() for query in args.query_lists}

    results_summary = []
    for query_name, query_list in query_dict.items():
    
        with log.section('Modeling {}:'.format(str(query_name))):
            try: 
                results, metadata = lisa.predict(query_list, **extract_kwargs(args, PREDICTION_KWARGS))

                top_TFs_unique = save_and_get_top_TFs(args, query_name, results, metadata)
            
                results_summary.append((query_name, top_TFs_unique))
            
            except AssertionError as err:
                log.append('ERROR: ' + str(err))

    print_results_multi(results_summary)

    #make_summary_table(args.prefix + '.combined.' + original_isd_names[args.isd_method] + '.csv'
def run_tests(args):

    if not args.skip_oneshot:
        tests.test_oneshot(args.test_genelist, args.background_genelist)

    tests.test_multi(args.genelists)
        
def build_common_args(parser):
    parser.add_argument('--seed', type = int, default = 2556, help = 'Random seed for gene selection. Allows for reproducing exact results.')
    parser.add_argument('--use_motifs', action = 'store_const', const = 'motifs', default='chipseq',
        dest = 'isd_method', help = 'Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database).')
    parser.add_argument('--save_metadata', action = 'store_true', default = False, help = 'Save json-formatted metadata from processing each gene list.')

def build_from_genes_args(parser, add_assays = True):
    #parser.add_argument('-c','--cores', required = True, type = int)
    if add_assays:
        parser.add_argument('-a','--assays',nargs='+',default=['Direct','H3K27ac','DNase'], choices=['Direct','H3K27ac','DNase'], help = 'Which set of insilico-deletion assays to run.')
    parser.add_argument('--rp_map_style', dest = 'rp_map', choices=public_config.get('lisa_params','rp_map_styles').split(','),
        default= public_config.get('lisa_params','rp_map_styles').split(',')[0], help = 'Which style of rp_map to assess influence of regions on genes. "basic" is stricly distance-based, while "enhanced" masks the exon and promoter regions of nearby genes.')

def build_multiple_lists_args(parser):
    parser.add_argument('query_lists', type = argparse.FileType('r', encoding = 'utf-8'), nargs = "+", help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
    parser.add_argument('-o','--output_prefix', required = True, type = is_valid_prefix, help = 'Output file prefix.')
    parser.add_argument('-v','--verbose',type = int, default = 2)
    parser.add_argument('-b','--num_background_genes', type = int, default = public_config.get('lisa_params', 'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes. These genes are selection from other gene lists.')
    parser.add_argument('--random_background', action = 'store_const', const = 'random', default = 'regulatory', dest = 'background_strategy', help = 'Use random background selection rather than "regulatory" selection.')
    
def build_one_list_args(parser, default_background_strategy = 'regulatory'):
    parser.add_argument('-o','--output_prefix', required = False, type = is_valid_prefix, help = 'Output file prefix. If left empty, will write results to stdout.')
    parser.add_argument('--background_strategy', choices = public_config.get('lisa_params', 'background_strategies').split(','),
        default = default_background_strategy,
        help = """Background genes selection strategy. LISA samples background genes to compare to user\'s genes-of-interest from a diverse
        regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
        """)
    background_genes_group = parser.add_mutually_exclusive_group()
    background_genes_group.add_argument('--background_list', type = argparse.FileType('r', encoding = 'utf-8'), required = False,
        help = 'user-supplied list of backgroung genes. Used when --background_strategy flag is set to "provided"')
    background_genes_group.add_argument('-b','--num_background_genes', type = int, default = public_config.get('lisa_params', 'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes')
    parser.add_argument('-v','--verbose',type = int, default = 4)

def confirm_file(arg):
    if os.path.isfile(arg):
        return arg
    else:
        raise argparse.ArgumentTypeError('ERROR: {} is not a valid file'.format(str(arg)))

class RstFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
        formatter_class=RstFormatter,
        description =
"""
Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6
X. Shirley Liu Lab, 2020
""")
parser.add_argument('--version', action = 'version', version = __version__)
subparsers = parser.add_subparsers(help = 'commands')

#__ LISA oneshot command __
oneshot_parser = subparsers.add_parser('oneshot', formatter_class=RstFormatter, description = '''
lisa oneshot
------------

You have:

* one genelist

Use LISA to infer influential TFs from one gene list, with background epigenetic landscape modeled using public data. 
If you have multiple lists, this option will be slower than using "multi" due to data-loading time. \n

Example::

    $ lisa oneshot hg38 ./genelist.txt -b 501 --seed=2556 --save_metadata > results.tsv

''')
oneshot_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
oneshot_parser.add_argument('query_list', type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
build_one_list_args(oneshot_parser)
build_from_genes_args(oneshot_parser)
build_common_args(oneshot_parser)
oneshot_parser.set_defaults(func = lisa_oneshot)

#__ LISA multi command __
multi_parser = subparsers.add_parser('multi', formatter_class=RstFormatter, description = '''
lisa multi
----------

You have:

* multiple genelists

Use LISA to infer influential TFs from multiple lists. This function processes each genelist independently in the same manner as the "oneshot" command, but reduces data loading time. Useful when performing 
the test on up and down-regulated genes from multiple RNA-seq clusters.

Example::

    $ lisa multi hg38 ./genelists/*.txt -b 501 -o ./results/

''')
multi_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
build_multiple_lists_args(multi_parser)
build_from_genes_args(multi_parser)
build_common_args(multi_parser)
multi_parser.set_defaults(func = lisa_multi, background_list = None)

from argparse import SUPPRESS
#____ LISA regions command ____
regions_parser = subparsers.add_parser('regions', formatter_class=RstFormatter, add_help = False, description = '''
lisa regions
------------

You have:

* one genelist
* regions (250 - 1000 bp wide) of interest related to that list
* optional: a positive score/weight associated with each region (you may pass zero-weight regions, but they do not affect the test and will be filtered out)

Use LISA to infer TF influence on your geneset, but provide your regions-of-interest rather than building a background epigenetic model using public data. When providing 
your own regions, LISA uses higher resolution, more precise binding data to increase the power of the test. Your regions should be between ~250 and 1000 bp in width, and the 
associated score should be positive. Scores are often read-depth at those regions, but can be any metic you think may influence gene regulation.

Example::

    $ lisa regions -r ./regions.bed -q ./genelist.txt --save_metadata > results.tsv
    $ lisa regions -r ./macs_peaks.xls -q ./genelist.txt --macs_xls > results.tsv

''')
regions_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
regions_required = regions_parser.add_argument_group('required arguments')
regions_required.add_argument('-q', '--query_genes', required = True, type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene list. One gene per line in either symbol or refseqID format')
regions_required.add_argument('-r', '--regions', type = confirm_file, required = True, help = 'Tad-delineated bed file with columns: chr, start, end[, score]. The score column is optional. If not provided, LISA will assign each region a uniform weight.')
regions_optional = regions_parser.add_argument_group('optional arguments')
regions_optional.add_argument('--header', action = 'store_true', default=False, help = 'Bed file has header row as first row. The header row may contain ')
regions_optional.add_argument('--macs_xls', action = 'store_true', default=False, help='If provided, regions file is a MACS2 .xls output file, and the "pileup" field is taken to be the region score.')
regions_optional.add_argument('--rp_map_style', dest = 'rp_map', choices=user_config.get('lisa_params','rp_map_styles').split(','),
    default=user_config.get('lisa_params','rp_map_styles').split(',')[0])
regions_optional.add_argument('--rp_decay', type = int, default = user_config.get('lisa_params','rp_decay'),
    help = 'Distance in base-pairs in which the influence of a region on a gene decays by half. Increase for more weight on distal elements, decrease for more weight on promoter elements.')
build_one_list_args(regions_optional, default_background_strategy='all')
build_common_args(regions_optional)
regions_optional.add_argument('-h', '--help', action = 'help', default=SUPPRESS)
regions_parser.set_defaults(func = lisa_regions)

#___ LISA coverage commands _____

coverage_parser = subparsers.add_parser('coverage', formatter_class = RstFormatter, add_help = False, description = '''
lisa coverage
------------

You have:

* one genelist
* bigwig of coverage over the genome

Use LISA to infer TF influence on your geneset using your own coverage data. This test is better suited than the "regions" test when your measure produces wide peaks/areas of influence.
An example of this is H3K27ac data, which correlates with gene expression similarly to accessibility, but produces wide peaks that may span many distinct TF binding locations.

Example::

    $ lisa coverage -bw ./sample.bigwig -q ./genelist.txt --save_metadata > results.tsv

''')
coverage_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
coverage_parser.add_argument('-q', '--query_genes', required = True, type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene list. One gene per line in either symbol or refseqID format')
coverage_parser.add_argument('-bw', '--bigwig_path', type = confirm_file, required = True, help = 'Bigwig file describing coverage over the genome.')
coverage_optional = coverage_parser.add_argument_group('optional arguments')
build_from_genes_args(coverage_optional, False)
build_one_list_args(coverage_optional, default_background_strategy='all')
build_common_args(coverage_optional)
coverage_optional.add_argument('-h', '--help', action = 'help', default=SUPPRESS)
coverage_parser.set_defaults(func = lisa_coverage)

#__ download command ___

def lisa_download(args):
    if args.command in ['oneshot','multi','coverage']:
        _class = FromGenes
    elif args.command == 'regions':
        _class = FromRegions
    else:
        raise AssertionError('Command {} not recognized'.format(args.command))
    if args.url:
        print(_class.get_dataset_url(args.species))
    else:
        _class.download_dataset(args.species)

download_data_parser = subparsers.add_parser('download', description = 'Download data from CistromeDB. Use if data recieved is incomplete or malformed.')
download_data_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Download data associated with human (hg38) or mouse (mm10) genes')   
download_data_parser.add_argument('command', choices=['oneshot', 'multi', 'regions', 'coverage'], help = 'For which command to download data')
download_data_parser.add_argument('--url', action = 'store_true', help = 'Get url for data download. Does not install data.')
download_data_parser.set_defaults(func = lisa_download)

#__ install command ___

def install_data(args):
    if args.command in ['oneshot','multi','coverage']:
        _class = FromGenes
    elif args.command == 'regions':
        _class = FromRegions
    else:
        raise AssertionError('Command {} not recognized'.format(args.command))
    if args.remove:
        os.rename(args.dataset, _class.get_dataset_path(args.species))
    else:
        copyfile(args.dataset, _class.get_dataset_path(args.species))

install_data_parser = subparsers.add_parser('install', description = 'Helper command for manually installing Lisa\'s data')
install_data_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Install data associated with human (hg38) or mouse (mm10) genes')   
install_data_parser.add_argument('command', choices=['oneshot', 'multi', 'regions', 'coverage'], help = 'For which command to install data')
install_data_parser.add_argument('dataset', type = confirm_file, help = 'Path to downloaded h5 dataset')
install_data_parser.add_argument('--remove', action = 'store_true', help = 'Delete dataset after installation is complete.')
install_data_parser.set_defaults(func = install_data)

#____ LISA run tests command ___
test_parser = subparsers.add_parser('run-tests')
test_parser.add_argument('species', type = str, choices=['hg38','mm10'])
test_parser.add_argument('test_genelist', type = confirm_file, help = 'test genelist for oneshot command')
test_parser.add_argument('background_genelist', type = confirm_file, help = 'background genelist for oneshot command')
test_parser.add_argument('genelists', nargs = '+', type = str, help = 'genelists for testing multi and one-vs-rest commands')
test_parser.add_argument('--skip_oneshot', action='store_true')
args = parser.parse_args()
test_parser.set_defaults(func = run_tests)
    
def main():
    #____ Execute commands ___
    args = parser.parse_args()

    try:
        args.func #first try accessing the .func attribute, which is empty if user tries ">>>lisa". In this case, don't throw error, display help!
    except AttributeError:
        print(parser.print_help(), file = sys.stderr)
    else:
        try:
            args.func(args)
        except (AssertionError, DownloadRequiredError, DatasetNotFoundError) as err:
            print(err, file = sys.stderr)