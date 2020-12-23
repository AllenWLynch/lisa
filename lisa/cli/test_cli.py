from itertools import product
import subprocess
import os
import sys
import argparse

TEST_PATH = os.path.join('blah', 'test_files')

common_options = dict(
        use_motifs = ['','--use_motifs'],
        save_metadata = ['', '--save_metadata'],
        assays = ['--assays ' + a for a in ['DNase H3K27ac Direct', 'Direct', 'DNase']],
        verbose = '-v 3',
        cores = '-c 10',
    )

new_test_badge = '''
***********************************
Test #{i}
{command}
***********************************
'''

def permute_options(options):
    all_permutations = []
    for permute_option, values in options.items():
        if isinstance(values, list) and len(values) > 1:
            permutations = values[1:]
            for option, values in options.items():
                if option != permute_option and values[0] != '':
                    if isinstance(values, list):
                        add_option = values[0]
                    else: 
                        add_option = values
                    permutations = [permutation + ' ' + add_option for permutation in permutations]
            all_permutations.extend(permutations)
    return all_permutations

def exec_commands(commands, run = 'all'):

    failed_commands = 0
    for i, command in enumerate(commands):
        print(new_test_badge.format(command = command, i = str(i + 1)))
        try:
            process = subprocess.run(command.split(' '), capture_output=True)
            if process.returncode!=0:
                raise AssertionError()
        except AssertionError:
            print(process.stderr.decode('utf-8'))
            failed_commands += 1
            pass
        else:
            print('Passed!')

    return failed_commands

def run_tests(base_command, options):
    
    cmd_combos = [base_command + ' ' + option for option in permute_options(options)]

    failed = exec_commands(cmd_combos)

    print('Passed {}/{} tests.'\
        .format(str(len(cmd_combos) - failed), str(len(cmd_combos))))

def test_oneshot(test_genelist, background_list):

    oneshot_options = dict(
        background_strategy = ['--background_strategy regulatory -b 501', 
            '--background_strategy random -b 300', 
            '--background_strategy provided --background_list {}'.format(background_list)],
        output_prefix = ['', '-o results'],
        **common_options
    )

    run_tests('lisa oneshot hg38 ' + test_genelist, oneshot_options)
    run_tests('lisa oneshot mm10 ' + test_genelist, oneshot_options)

def test_multi(genelists):

    mulit_options = dict(
        num_background_genes = '-b 501',
        use_random_background = ['', '--random_background'],
        output = '-o results/',
        **common_options
    )

    run_tests('lisa multi hg38 ' + ' '.join(genelists), mulit_options)
    run_tests('lisa multi mm10 ' + ' '.join(genelists), mulit_options)