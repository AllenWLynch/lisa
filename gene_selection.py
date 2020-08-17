
import numpy as np
import os
import random
from collections import Counter, defaultdict

#enforce sampling according to the TAD distribution of the genome.
def sampling_by_tad(background_gene_list, gene_names, TAD_data, num_selected=3000):
    '''
    background_gene_list: list, genes from user-supplied list of background candidates, or genes that are background candidates because they are not in query
    num_selected: number of genes to sample for background
    '''
    assert(num_selected < len(TAD_data)), 'Number of genes to sample must be less than the number of background gene candidates'

    #intersect background gene list with TAD_list to eliminate reserved genes
    TAD_data = [(gene_name, tad_group) for (gene_name, tad_group) in zip(gene_names, TAD_data) if gene_name in background_gene_list]

    #collect list of genes in each bin
    genes_in_TAD = defaultdict(list)
    for gene_name, tad_group in TAD_data:
        genes_in_TAD[tad_group].append(gene_name)
    
    #calculates number of genes in each TAD group
    num_genes_in_TAD = {tad_group : len(genes) for tad_group, genes in genes_in_TAD.items()}

    #calculates the number of genes expected to be sampled from each TAD, with a bit of an enrichment (1.1x). Ensure all TADs have >= 1 expected genes
    expected_samples = {
        tad_group : max(1, int(num_genes / len(TAD_data) * num_selected * 1.1))
        for tad_group, num_genes in num_genes_in_TAD.items()
    }

    #samples the expected number of genes from the TAD (or the number of genes if expected is > actual)
    selected_genes = []
    for tad_group, num_expected in expected_samples.items():
        sampled_genes = np.random.choice( genes_in_TAD[tad_group], min(num_genes_in_TAD[tad_group] - 1, num_expected), replace = False)
        selected_genes.extend(sampled_genes)

    return selected_genes


#Randomly sample background genes from the whole genome
def random_background_genes(background_gene_list, num_selected = 3000):
    '''
    background_gene_list: list, genes from user-supplied list of background candidates, or genes that are background candidates because they are not in query
    num_selected: number of background genes to sample
    '''
    assert(len(background_gene_list) > num_selected), 'Background gene list provided must contain more than {} genes'.format(str(num_selected))
    
    #if same number of genes, skip sample
    if len(background_gene_list) == num_selected:
        return background_gene_list
    
    else:
        return np.random.choice(background_gene_list, num_selected, replace = False)


#recieve user gene-of-interest list, process text, validate, and convert to gene symbols
def match_user_provided_genes(user_genelist, gene_symbols, gene_refseqIDs):
    """
    user_genelist: list of user-supplied text, supposedly containing some gene-like entries
    gene_symbols: list of gene symbols
    gene_refseqIDs: list of refseqIDs corresponding to their respective gene symbol
    """
    user_genelist = [user_symbol.upper() for user_symbol in user_genelist]

    gene_symbols, gene_refseqIDs = np.array([symbol.upper() for symbol in gene_symbols]), np.array([ref_id.upper() for ref_id in gene_refseqIDs])

    supplied_symbols = np.isin(gene_symbols, user_genelist)
    supplied_ids = np.isin(gene_refseqIDs, user_genelist)

    user_supplied_genes = np.logical_or(supplied_symbols, supplied_ids)

    return list(gene_symbols[user_supplied_genes])

#creates a label vector, indexed by gene symbol, with binary labels: {0 = background, 1 = query}
def create_label_dictionary(query_list, background_list):
    return { gene_symbol : i < len(query_list) for i, gene_symbol in enumerate(list(query_list) + list(background_list)) }


#converts user-entered data to two gene lists: a query list, and a background list
def select_genes(query_list, gene_symbols, gene_ids, tad_domains, num_background_genes = 3000, 
    background_strategy = 'regulatory', max_query_genes = 200, background_list = []):
    '''
    query_list: text-processed list of symbols or refseqIDs (or both) that user supplies
    gene_symbols: list of recognized gene_symbols
    gene_refseqIDs: list of corresponding refseqIDs
    tad_domains: tad and k4me3 cluster data for each gene. Determines distinct regulatory type
    num_background_genes: number of background genes to select
    background_strategy: regulatory, provided, or random background selection strategy
    user_background_genes: text-processed list of symbols or refseqIDs (or both) that user supplies as candidate background genes
    max_query_genes: config setting
    background_list: text-processed list of symbols or refseqIDs (or both)
    '''
                
    query_genes = match_user_provided_genes(query_list, gene_symbols, gene_ids)

    assert(20 <= len(query_genes) <= max_query_genes), 'User must provide list of 20 to {} unique genes. Provided {}'\
        .format(str(max_query_genes), str(len(query_genes)))

    if background_strategy == 'provided':

        background_genes = match_user_provided_genes(background_list, gene_symbols, gene_ids)

        assert( len(set(background_genes).intersection(set(query_genes))) == 0), 'Some genes in your query are also in your background genes'
        assert( len(background_genes) > len(query_genes) )

    else:

        assert( num_background_genes >= len(query_genes) and len(query_genes) <= 19000//2 ), "More query genes selected than background genes"
        background_candidates = list(set(gene_symbols).difference(query_genes))
        assert(len(background_candidates) >= num_background_genes), 'Number of background candidates must be greater than or equal number of genes to select as background.'

        #if no down-sampling is needed:
        if background_candidates == num_background_genes:
            background_genes = background_candidates

        elif background_strategy == 'regulatory':
            background_genes = sampling_by_tad(background_candidates, gene_symbols, tad_domains, num_selected=num_background_genes)

            if len(background_genes) > num_background_genes:
                background_genes = random_background_genes(background_genes, num_background_genes)
            
        elif background_strategy == 'random':
            background_genes = random_background_genes(background_candidates, num_background_genes)

        else:
            raise AssertionError('Background selection strategy {} not supported'.format(background_strategy))

    return create_label_dictionary(query_genes, background_genes)