
import numpy as np
import os
import random
from collections import Counter, defaultdict

#ADD DUPLICATE CHECKING

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
def create_label_vector(query_list, background_list):
    return list(query_list) + list(background_list), np.concatenate([np.ones(len(query_list)), np.zeros(len(background_list))])


#converts user-entered data to two gene lists: a query list, and a background list
def select_genes_for_chromatin_model(query_genes, gene_symbols, gene_refseqIDs, TAD_data, num_selected = 3000, user_background_genes = None, method = 'TAD'):
    '''
    query_genes: text-processed list of symbols or refseqIDs (or both) that user supplies
    gene_symbols: list of recognized gene_symbols
    gene_refseqIDs: list of corresponding refseqIDs
    TAD_data: tad and k4me3 cluster data for each gene. Determines distinct regulatory type
    num_selected: number of background genes to select
    user_background_genes: text-processed list of symbols or refseqIDs (or both) that user supplies as candidate background genes
    method: TAD or random, background selection strategy
    '''
    assert(method == 'TAD' or method == 'random'), 'Method may be either "TAD" or "random". User supplied {}'.format(method)

    #process and validate user-supplied query
    query_list = match_user_provided_genes(query_genes, gene_symbols, gene_refseqIDs)

    num_genes_supplied = len(query_list)
    assert(20 <= num_genes_supplied <= 200), 'User must provide list of between 20 and 200 genes. Provided {}'.format(int(num_genes_supplied))

    #if the user supplied a background candidates list, process that, otherwise take the background candidates to be all genes not in query
    if not user_background_genes is None:
        background_candidates = match_user_provided_genes(user_background_genes, gene_symbols, gene_refseqIDs)
    else:
        background_candidates = list(set(gene_symbols).difference(query_list))

    assert(len(background_candidates) >= num_selected), 'Number of background candidates must be greater than or equal number of genes to select as background.'

    #if no down-sampling is needed:
    if background_candidates == num_selected:
        return background_candidates

    if method == 'TAD':
        background_list = sampling_by_tad(background_candidates, gene_symbols, TAD_data, num_selected=num_selected)
        #this may have a number greater than num selected
        if len(background_list) > num_selected:
            background_list = np.random.choice(background_list, num_selected, replace = False)
    else:
        background_list = random_background_genes(background_candidates, num_selected)

    #get list of genes selected for query and background and get labels for chromatin model
    gene_names, label_vector = create_label_vector(query_list, background_list)

    return gene_names, label_vector