
import pandas as pd
from lisa import FromGenes
import requests
import fire
import argparse
import numpy as np

URL = 'http://lisa.cistrome.org/gallery/{id}_{dir}.txt'

TEST = [
    ('189', ['NFE2']),
    ('172',['PPARG']),
    ('187',['YY1']),
    ('220',['RUNX3']),
    ('135', ['HNF4A']),
    ('154',['MYOD1']),
    ('173',['EBF1']),
    ('182',['NR1D1']),
    ('263',['ESR1']),
    ('264',['GATA3']),
    ('183',['NKX2-5']),
    ('163',['GATA4']),
    ('280',['MED1']),
    ('37',['CEBPG']),
    ('304',['RUNX1']),
    ('275', ['BMI1']),
    ('222', ['BCL6']),
    ('292',['SOX2']),
    ('20',['IRF4']),
    ('181', ['FOXP3'])
]

def get_ranks(factor_name, results, rank_column = 'summary_p_value', max_rank = 500):

    results = pd.DataFrame(results.to_dict())

    results['col_rank'] = results[rank_column].rank()

    ranks = np.array(results[results.factor == factor_name].col_rank.values)

    ranks = np.clip(ranks, None, max_rank + 1)

    return list(ranks)

def get_genelist(cistrome_id, direction):

    r = requests.get(URL.format(id = cistrome_id, dir = direction))

    return [x.strip() for x in r.text.split('\n')]


def run_test(lisa_obj, factor_name, cistrome_id, direction, max_rank = 500, rank_column = 'summary_p_value'):
    
    genelist = get_genelist(cistrome_id, direction)

    try:

        results, _ = lisa_obj.predict(genelist, num_background_genes = 750, seed=2556)

        ranks = get_ranks(factor_name, results, rank_column=rank_column, max_rank=max_rank)

        print(factor_name + ': ', str(ranks))

        return ranks
    except Exception as err:
        print(repr(err))
        return []

def run_tests(lisa_obj, tests, max_rank = 500, rank_column = 'summary_p_value'):

    all_ranks = []
    for test in tests:

        for direction in ['up','down']:
            ranks = run_test(lisa_obj, test[1][0], test[0], direction, max_rank= max_rank, rank_column=rank_column)
            all_ranks.extend(ranks)

    return all_ranks


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('output', type = argparse.FileType('w'))
    args = parser.parse_args()

    lisa_obj = FromGenes('mm10', assays = ['DNase'], verbose=1)

    all_ranks = run_tests(lisa_obj, TEST, max_rank=500, rank_column='DNase_p_value')

    print('\n'.join([str(r) for r in all_ranks]), end = '', file = args.output)