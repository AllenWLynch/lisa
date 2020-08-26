
import numpy as np
import os
import random
from collections import Counter, defaultdict, OrderedDict
from utils import GeneSet, Gene

#enforce sampling according to the TAD distribution of the genome.
def sampling_by_tad(background_gene_list, num_selected=3000):
    '''
    background_gene_list: list, genes from user-supplied list of background candidates, or genes that are background candidates because they are not in query
    num_selected: number of genes to sample for background
    '''
    #intersect background gene list with TAD_list to eliminate reserved genes
    TAD_data = [(gene, gene.tad_domain) for gene in background_gene_list]

    #collect list of genes in each bin
    genes_in_TAD = defaultdict(list)
    for gene, tad_group in TAD_data:
        genes_in_TAD[tad_group].append(gene)
    
    #calculates number of genes in each TAD group
    num_genes_in_TAD = {tad_group : len(genes) for tad_group, genes in genes_in_TAD.items()}

    #calculates the number of genes expected to be sampled from each TAD, with a bit of an enrichment (1.1x). Ensure all TADs have >= 1 expected genes
    expected_samples = {
        tad_group : max(1, int(num_genes / len(TAD_data) * num_selected * 1.1))
        for tad_group, num_genes in num_genes_in_TAD.items()
    }

    #samples the expected number of genes from the TAD (or the number of genes if expected is > actual)
    selected_genes = GeneSet()
    for tad_group, num_expected in expected_samples.items():
        sampled_genes = np.random.choice( genes_in_TAD[tad_group], min(num_genes_in_TAD[tad_group] - 1, num_expected), replace = False)
        selected_genes.add_genes(sampled_genes)

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
        random_sample = GeneSet()
        random_sample.add_genes(np.random.choice(list(background_gene_list), num_selected, replace = False))

        return random_sample


#recieve user gene-of-interest list, process text, validate, and convert to gene symbols
def match_user_provided_genes(user_genelist, gene_set):
    """
    user_genelist: list of user-supplied text, supposedly containing some gene-like entries
    gene_set
    """
    selected_genes = GeneSet()
    for gene_candidate in user_genelist:
        try:
            selected_genes.add_gene( gene_set.get_gene_by_name(gene_candidate) )
        except KeyError:
            pass
    
    return selected_genes.get_distinct_genes_by_symbol()

#creates a label vector, indexed by gene symbol, with binary labels: {0 = background, 1 = query}
def create_label_dictionary(query_list, background_list):
    return {
        gene.get_location() : int(i < len(query_list))
        for i, gene in enumerate(list(query_list) + list(background_list))
    }, query_list.get_symbols(), background_list.get_symbols()
    
#converts user-entered data to two gene lists: a query list, and a background list
def select_genes(query_list, gene_set, num_background_genes = 3000, 
    background_strategy = 'regulatory', max_query_genes = 200, background_list = []):
    '''
    query_list: text-processed list of symbols or refseqIDs (or both) that user supplies
    num_background_genes: number of background genes to select
    background_strategy: regulatory, provided, or random background selection strategy
    user_background_genes: text-processed list of symbols or refseqIDs (or both) that user supplies as candidate background genes
    max_query_genes: config setting
    background_list: text-processed list of symbols or refseqIDs (or both)
    '''
                
    query_genes = match_user_provided_genes(query_list, gene_set)

    assert(20 <= len(query_genes) <= max_query_genes), 'User must provide list of 20 to {} unique genes. Provided {}'\
        .format(str(max_query_genes), str(len(query_genes)))

    if background_strategy == 'provided':

        background_genes = match_user_provided_genes(background_list, gene_set)

        assert( len(set(background_genes.get_symbols()).intersection(set(query_genes.get_symbols()))) == 0), 'Some genes in your query are also in your background genes'
        assert( len(background_genes) > len(query_genes) )

    else:

        assert( num_background_genes >= len(query_genes) and len(query_genes) <= 19000//2 ), "More query genes selected than background genes"

        background_candidates = gene_set.get_distinct_genes_by_symbol(excluding = query_genes.get_symbols())

        assert(len(background_candidates) >= num_background_genes), 'Number of background candidates must be greater than or equal number of genes to select as background.'

        #if no down-sampling is needed:
        if background_candidates == num_background_genes:
            background_genes = background_candidates

        elif background_strategy == 'regulatory':

            background_genes = sampling_by_tad(background_candidates, num_selected=num_background_genes)

            if len(background_genes) > num_background_genes:
                background_genes = random_background_genes(background_genes, num_background_genes)
            
        elif background_strategy == 'random':
            background_genes = random_background_genes(background_candidates, num_background_genes)

        else:
            raise AssertionError('Background selection strategy {} not supported'.format(background_strategy))

    #print(type(query_genes), len(query_genes), type(background_genes), len(background_genes))

    return create_label_dictionary(query_genes, background_genes)