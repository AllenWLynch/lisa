
import numpy as np
import pandas as pd
import os
import random

#ADD DUPLICATE CHECKING 

#enforce sampling according to the TAD distribution of the genome.
def sampling_by_tad(background_gene_list, TAD_data, num_selected=3000):
    '''
    background_gene_list: list, genes from user-supplied list of background candidates, or genes that are background candidates because they are not in query
    TAD_data: TAD cluster and regulatory group for each gene
    num_selected: number of genes to sample for background
    '''
    assert(num_selected < len(TAD_data)), 'Number of genes to sample must be less than the number of background gene candidates'

    #intersect background gene list with TAD_list to eliminate reserved genes
    TAD_data = TAD_data[TAD_data.geneName.isin(background_gene_list)]
    
    #create bins for distinct regulatory states (combination of TAD and methylation states)
    TAD_data['bin'] = TAD_data.k4me3_order_cluster + TAD_data.tad_order_cluster

    #collect list of genes in each bin
    aggregated_by_TAD = TAD_data.groupby('bin')['geneName'].apply(list).rename(columns = {'list' : 'genes_in_tad'})

    #count number of genes in each bin
    aggregated_by_TAD['num_genes_in_TAD'] = TAD_data.genes_in_tad.str.len()

    #calculates the number of genes expected to be sampled from each TAD, with a bit of an enrichment (1.1x). Ensure all TADs have >= 1 expected genes
    aggregated_by_TAD['expected_samples'] = np.maximum(int(aggregated_by_TAD.num_genes_in_TAD / len(TAD_data) * num_selected * 1.1), 1.0)

    #samples the expected number of genes from the TAD (or the number of genes if expected is > actual)
    aggregated_by_TAD['sampled_genes'] = aggregated_by_TAD[aggregated_by_TAD.num_genes_in_TAD > 0].apply(
        lambda tad : np.random.choice(tad.genes_in_tad, np.min(tad.num_genes_in_TAD, tad.expected_samples), replace = False)
    )
    #flatten the list to yield background gene selections
    return set([gene_name for gene_list in aggregated_by_TAD.sampled_genes.values for gene_name in gene_list])


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
def process_user_genelist(user_genelist, gene_metadata)

    supplied_symbols = gene_metadata.geneName.isin(user_genelist)
    supplied_ids = gene_metadata.refseqID.isin(user_genelist)

    user_supplied_genes = np.logical_and(supplied_symbols, supplied_ids)

    return list(gene_metadata[user_supplied_genes].geneName.values)

#creates a label vector, indexed by gene symbol, with binary labels: {0 = background, 1 = query}
def create_label_vector(query_list, background_list):
    return pd.Series(
        index = query_list + background_list,
        values = np.concatenate(
            [np.ones(len(query_list)), np.zeros(len(background_list))]
        )
    )


#converts user-entered data to two gene lists: a query list, and a background list
def process_gene_lists(query_genes, gene_metadata, TAD_data, num_selected = 3000, user_background_genes = None, method = 'TAD'):
    '''
    query_genes: text-processed list of symbols or refseqIDs (or both) that user supplies
    gene_metdata: dataframe of chr, start, end, refseq, symbol
    num_selected: number of background genes to select
    user_background_genes: text-processed list of symbols or refseqIDs (or both) that user supplies as candidate background genes
    method: TAD or random, background selection strategy
    '''
    assert(method == 'TAD' or method == 'random'), 'Method may be either "TAD" or "random". User supplied {}'.format(method)

    #process and validate user-supplied query
    query_list = process_user_genelist(query_genes, gene_metadata)

    num_genes_supplied = len(query_list)
    assert(20 <= num_genes_supplied <= 200), 'User must provide list of between 20 and 200 genes. Provided {}'.format(int(user_supplied_genes.sum()))

    #if the user supplied a background candidates list, process that, otherwise take the background candidates to be all genes not in query
    if not user_background_genes is None:
        background_candidates = process_user_genelist(user_background_genes, gene_metadata)
    else:
        background_candidates = list(gene_metadata[ ~gene_metadata.geneName.isin(query_list)].geneName.values)

    assert(len(background_candidates) >= num_selected), 'Number of background candidates must be greater than or equal number of genes to select as background.'

    #if no down-sampling is needed:
    if background_candidates == num_selected:
        return background_candidates

    if method == 'TAD':
        background_list = sampling_by_tad(background_candidates, TAD_data, num_selected=num_selected)
        #this may have a number greater than num selected
        if len(background_candidates) > num_selected:
            background_list = np.random.choice(background_candidates, num_selected, replace = False)
    else:
        background_list = random_background_genes(background_candidates, num_selected)

    return create_label_vector(query_list, background_list)

    




















    


