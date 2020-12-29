'''
********
Lisa CLI
********

Installing LISA using pip or conda adds the "lisa" command to your path. LISA's functionality is divided into three main subcommands:

* `lisa oneshot`_ : one genelist
* `lisa multi`_ : multiple genelists
* `lisa regions`_ : one genelist and a list of regions

Which are used depending on the evidence you have on hand. 

See the `User Guide <docs/user_guide.rst>`_ for more usage information.
See the `Python API <docs/python_api.rst>`_ for more in-depth description of tests and parameters.

'''

from lisa import FromRegions, FromGenes
from lisa.core.utils import Log
from lisa._version import __version__
from lisa.core.lisa_core import INSTALL_PATH, DownloadRequiredError
import configparser
import argparse
import os
import sys
import json
from collections import defaultdict
from shutil import copyfile
import lisa.cli.test_cli as tests

from lisa.lisa_public_data.lisa import _config as public_config
from lisa.lisa_user_data.lisa import _config as user_config

#____COMMAND LINE INTERFACE________

INSTANTIATION_KWARGS = ['cores','isd_method','verbose','assays']
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

def lisa_download(args):
    lisa = FromGenes(args.species)
    if args.url:
        print(lisa.get_dataset_url())
    else:
        lisa.download_data()

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


def save_and_get_top_TFs(args, query_name, results, metadata):

    with open(args.output_prefix + query_name + '.lisa.tsv', 'w') as f:
        f.write(results.to_tsv())

    if args.save_metadata:
        with open(args.output_prefix + query_name + '.metadata.json', 'w') as f:
            f.write(json.dumps(metadata, indent=4))

    try:
        top_TFs = results.filter_rows(lambda x : x <= 0.01, 'summary_p_value').todict()['factor']
    except KeyError:
        top_TFs = ['None']

    top_TFs_unique, encountered = [], set()
    for TF in top_TFs:
        if not TF in encountered:
            top_TFs_unique.append(TF)
            encountered.add(TF)

    return top_TFs_unique


def print_results_multi(results_summary):
    print('Sample\tTop Regulatory Factors (p < 0.01)')
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


def lisa_backcompatible(args):

    def make_original_csv_format(filename, results, p_val_col):

        with open(filename, 'w') as f:

            if args.isd_method == 'chipseq':
                print(',pval,species,factor,factor_type,cell_line,cell_type,tissue,qc', file = f)
                index_with = ['sample_id','factor',p_val_col,'species','factor','factor_type','cell_line','cell_type','tissue','qc']
            else:
                print(',pval,edition,source,sourcefile,status,numseqs,pmid', file = f)
                index_with = ['sample_id','factor',p_val_col,'edition','source','sourcefile','status','numseqs','pmid']
            
            for x in zip(*results.get_column(index_with)):
                print(x[0]+'|'+x[1],*x[2:],sep = ',', file = f)

    def make_model_dataframe(filename, model_dir, metdata):
        sample_id = metadata[model_dir]['selected_datasets']['sample_id']
        coefs = metadata[model_dir]['chromatin_model']['coefs']
        cell_type = metadata[model_dir]['selected_datasets']['cell_type']
        cell_line = metadata[model_dir]['selected_datasets']['cell_line']
        tissue = metadata[model_dir]['selected_datasets']['tissue']

        with open(filename, 'w') as f:
            print(',coefficients,cell_line,cell_type,tissue', file = f)
            for line in zip(sample_id, coefs, cell_type, cell_line, tissue):
                print(*line, sep = ',', file = f)

    def make_summary_table(filename, results):
        factor_pvals = defaultdict(list)
        for factor, pval in zip(*results.get_column(['factor','summary_p_value'])):
            factor_pvals[factor].append(pval)

        with open(filename, 'w') as f:
            print('Transcription Factor,1st Sample p-value,2nd Sample p-value,3rd Sample p-value,4th Sample p-value,5th Sample p-value', file = f)
            for factor, pval_list in factor_pvals.items():
                print(factor, *pval_list[:5],sep = ',',file = f)

    if args.isd_method == 'knockout':
        args.isd_method = 'chipseq'

    lisa = FromGenes(args.species, **extract_kwargs(args, INSTANTIATION_KWARGS))
    results, metadata = lisa.predict(args.query_list.readlines(), **extract_kwargs(args, PREDICTION_KWARGS))
    
    if args.save_metadata:
        with open(args.prefix + '.metadata.json', 'w') as f:
                f.write(json.dumps(metadata, indent=4))
        with open(args.prefix + '.lisa.tsv', 'w') as f:
            f.write(results.to_tsv())

    results = results.add_column('factor_type',['tf']*len(results))\
                .add_column('species',[args.species]*len(results)).add_column('qc',['1']*len(results))

    original_isd_names = {'chipseq' : 'chipseq', 'motifs' : 'motif99'}
    
    make_original_csv_format(args.prefix + '.DNase.' + original_isd_names[args.isd_method] + '.p_value.csv', results, 'DNase_p_value')
    make_original_csv_format(args.prefix + '.H3K27ac.' + original_isd_names[args.isd_method] + '.p_value.csv', results, 'H3K27ac_p_value')
    
    original_isd_names = {'chipseq' : 'chipseq', 'motifs' : 'motif'}

    make_original_csv_format(args.prefix + '.' + original_isd_names[args.isd_method] + '_direct.p_value.csv', results, lisa.isd_method+'_p_value')
    make_original_csv_format(args.prefix + '_' + original_isd_names[args.isd_method] + '_cauchy_combine_raw.csv', results, 'summary_p_value')

    make_model_dataframe(args.prefix + '.DNase.coefs.csv','DNase',metadata)
    make_model_dataframe(args.prefix + '.H3K27ac.coefs.csv','H3K27ac',metadata)

    #make_summary_table(args.prefix + '.combined.' + original_isd_names[args.isd_method] + '.csv', results)

def run_tests(args):

    if not args.skip_oneshot:
        tests.test_oneshot(args.test_genelist, args.background_genelist)

    tests.test_multi(args.genelists)
        
def build_common_args(parser):
    parser.add_argument('--seed', type = int, default = 2556, help = 'Random seed for gene selection. Allows for reproducing exact results.')
    parser.add_argument('--use_motifs', action = 'store_const', const = 'motifs', default='chipseq',
        dest = 'isd_method', help = 'Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database).')
    parser.add_argument('--save_metadata', action = 'store_true', default = False, help = 'Save json-formatted metadata from processing each gene list.')

def build_from_genes_args(parser):
    parser.add_argument('-c','--cores', required = True, type = int)
    parser.add_argument('-a','--assays',nargs='+',default=['Direct','H3K27ac','DNase'], choices=['Direct','H3K27ac','DNase'])

def build_multiple_lists_args(parser):
    parser.add_argument('query_lists', type = argparse.FileType('r', encoding = 'utf-8'), nargs = "+", help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
    parser.add_argument('-o','--output_prefix', required = True, type = is_valid_prefix, help = 'Output file prefix.')
    parser.add_argument('-v','--verbose',type = int, default = 2)
    parser.add_argument('-b','--num_background_genes', type = int, default = public_config.get('lisa_params', 'background_genes'),
        help = 'Number of sampled background genes to compare to user-supplied genes. These genes are selection from other gene lists.')
    parser.add_argument('--random_background', action = 'store_const', const = 'random', default = 'regulatory', dest = 'background_strategy', help = 'Use random background selection rather than "regulatory" selection.')
    

def build_one_list_args(parser):
    parser.add_argument('-o','--output_prefix', required = False, type = is_valid_prefix, help = 'Output file prefix. If left empty, will write results to stdout.')
    parser.add_argument('--background_strategy', choices = public_config.get('lisa_params', 'background_strategies').split(','),
        default = 'regulatory',
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

def extract_tar_func(args):
    FromGenes.extract_tar(args.tarball, INSTALL_PATH, rm = args.rm)

def lisa_regions(args):

    try:
        args.background_list = args.background_list.readlines()
    except AttributeError:
        pass

    results, metadata = FromRegions(args.species, args.regions, isd_method = args.isd_method)\
        .predict(args.query_list.readlines(), **extract_kwargs(args, PREDICTION_KWARGS))
    
    save_results(args, results, metadata)

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

    $ lisa oneshot hg38 ./genelist.txt -b 501 -c 5 --seed=2556 --save_metadata > results.tsv

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

    $ lisa multi hg38 ./genelists/*.txt -b 501 -c 5 -o ./results/

''')
multi_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
build_multiple_lists_args(multi_parser)
build_from_genes_args(multi_parser)
build_common_args(multi_parser)
multi_parser.set_defaults(func = lisa_multi, background_list = None)

#__ download command ___
download_data_parser = subparsers.add_parser('download', help = 'Download data from CistromeDB. Use if data recieved is incomplete or malformed.')
download_data_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Download data associated with human (hg38) or mouse (mm10) genes')   
download_data_parser.add_argument('--url', action = 'store_true', help = 'Get url for data download. Does not install data.')
download_data_parser.set_defaults(func = lisa_download)

#__ install command __
install_data_parser = subparsers.add_parser('unpack', help = 'Helper command for manually installing Lisa\'s data')
install_data_parser.add_argument('tarball', type = confirm_file, help = 'Path to downloaded data tarball')
install_data_parser.add_argument('--remove', dest = 'rm', action = 'store_true', help = 'Remove data tarball after unpacking is complete.')
install_data_parser.set_defaults(func = extract_tar_func)

#__ LISA back-compatible interface __
def get_background_option(arg):
    if os.path.isfile(arg):
        return argparse.FileType('r')(arg)
    elif arg == 'dynamic_auto_tad':
        return []
    else:
        raise argparse.ArgumentTypeError('ERROR: {} is neither a background genelist file, nor an allowed background selection option.'.format(arg)) 

backcompat_parser = subparsers.add_parser('backcompat', help = 'Interface with LISA using the version 1 command line constructs')
backcompat_parser.add_argument('query_list', type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene lists. One gene per line in either symbol or refseqID format')
backcompat_parser.add_argument('--species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
backcompat_parser.add_argument('--threads', dest = 'cores', required = True, type = int, default = 1)
backcompat_parser.add_argument('--seed', type = int, default = 2556, help = 'Random seed for gene selection. Allows for reproducing exact results. LISA may be more sensitive to seed when background count is low.')
backcompat_parser.add_argument('--prefix', required=True, type = is_valid_prefix, help = 'Output file prefix.')
backcompat_parser.add_argument('--stat_background_number', dest='num_background_genes', type = int, default = public_config.get('lisa_params', 'background_genes'),
    help = 'Number of sampled background genes to compare to user-supplied genes')
backcompat_parser.add_argument('--background', dest = 'background_list', type = get_background_option, help = 'background genes provided as file, or type of background selection')
backcompat_parser.add_argument('--method', choices=['knockout', 'chipseq', 'motifs'], dest = 'isd_method', required=True,
    help = 'Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database).')
#absorb all old and deprecated flags
for flag in 'clean,web,new_rp_h5,new_count,cluster,random,covariates,new_count_h5'.split(','):
    backcompat_parser.add_argument('--' + flag, required=False, help = 'deprecated option')
backcompat_parser.add_argument('-v','--verbose',type = int, default = 4)
backcompat_parser.set_defaults(func = lisa_backcompatible, background_strategy = 'regulatory')
backcompat_parser.add_argument('--save_metadata', action = 'store_true', default = False, help = 'Save json-formatted metadata from processing each gene list.')
parser.add_argument('--epigenome',nargs='+',default=['Direct','H3K27ac','DNase'],dest='assays',choices=['Direct','H3K27ac','DNase'])

#____ LISA run tests command ___
test_parser = subparsers.add_parser('run-tests')
test_parser.add_argument('species', type = str, choices=['hg38','mm10'])
test_parser.add_argument('test_genelist', type = confirm_file, help = 'test genelist for oneshot command')
test_parser.add_argument('background_genelist', type = confirm_file, help = 'background genelist for oneshot command')
test_parser.add_argument('genelists', nargs = '+', type = str, help = 'genelists for testing multi and one-vs-rest commands')
test_parser.add_argument('--skip_oneshot', action='store_true')
args = parser.parse_args()
test_parser.set_defaults(func = run_tests)

#____ LISA regions command ____
regions_parser = subparsers.add_parser('regions', formatter_class=RstFormatter, description = '''
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

    $ lisa regions -r ./regions.bed -q ./genelist.txt -b 501 --save_metadata > results.tsv

''')
regions_parser.add_argument('species', choices = ['hg38','mm10'], help = 'Find TFs associated with human (hg38) or mouse (mm10) genes')
regions_parser.add_argument('-r', '--regions', type = str, required = True, help = 'Bed file with columns: chr, start, end[, score]')
regions_parser.add_argument('-q' '--query_list', required = True, type = argparse.FileType('r', encoding = 'utf-8'), help = 'user-supplied gene list. One gene per line in either symbol or refseqID format')
build_one_list_args(regions_parser)
build_common_args(regions_parser)
regions_parser.set_defaults(func = lisa_regions)
    
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
        except (AssertionError, DownloadRequiredError) as err:
            print(err, file = sys.stderr)