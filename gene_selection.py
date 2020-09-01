
import numpy as np
import os
import random
from collections import Counter, defaultdict, OrderedDict
import random

class Gene:
    '''
    A gene has a unique genomic region, an accepted name, and aliases that also correspond to that gene or genomic region
    '''
    def __init__(self, chrom, tss_start, tss_end, names, tad_domain = None):
        self.chrom = chrom
        self.start = int(tss_start)
        self.end = int(tss_end)
        self.aliases = []
        self.tad_domain = tad_domain
        if isinstance(names, str):
            self.add_alias(names)
        else:
            self.add_aliases(names)
        
    def add_alias(self, alias, is_name = True):
                
        if not alias in self.aliases:
            self.aliases.append(alias)

    def get_name(self):
        return self.aliases[0]

    def add_aliases(self, aliases):
        assert( isinstance(aliases, list) )
        for alias in aliases:
            self.add_alias(alias)

    def get_location(self):
        return ':'.join([str(property_) for property_ in (self.chrom, self.start, self.end)])
    
    def __eq__(self, other):
        if isinstance(other, str):
            return other in self.aliases
        elif isinstance(other, Gene):
            return self.get_location() == other.get_location()
        else:
            return False

    def __repr__(self):
        return '\t'.join([str(x) for x in [self.get_location(), self.aliases[0], self.tad_domain, '|'.join(self.aliases)]])

    def __str__(self):
        return self.get_name()

    def get_RP_signature(self, bins, bin_index, delta = 10000, max_influence_distance = 100000):

        tss = self.start
        
        #find bin regions defining interesting RP region
        min_bin, split_bin, max_bin = np.digitize([tss - max_influence_distance, tss, tss + max_influence_distance], bins)
        #subset interesting region of chrom
        bin_indices = bin_index[min_bin: max_bin - 1]
        #split the bin holding the TSS into two bins
        bins = np.concatenate([bins[min_bin: split_bin], [tss], bins[split_bin: max_bin]])
        split_bin -= min_bin
        tss_bins = (split_bin - 1, split_bin)

        #get bin intervals to the left and right of TSS, then concatenate
        left_bins = np.abs(np.array(list(zip(bins[1:tss_bins[1] + 1], bins[:tss_bins[0] + 1]))) - tss)
        right_bins = np.abs(np.array(list(zip(bins[tss_bins[1]:-1], bins[tss_bins[1] + 1:]))) - tss)
        intervals = np.concatenate([left_bins, right_bins], axis = 0)

        #get integral of RP at boundary locations
        RP = intervals * (-np.log(1/3) / delta)
        RP = 2 * ( RP - np.log(np.exp(RP) + 1))
        #compute RP over area
        RP = np.subtract(RP[:,1], RP[:,0])
        #sum bins split by TSS
        summed_tss_rp = RP[list(tss_bins)].sum()
        RP[tss_bins[0]] = summed_tss_rp
        #remove split bin
        RP = np.delete(RP, tss_bins[1])
        #normalize so bin with TSS has RP of 1
        RP = RP / RP.max()

        return RP, bin_indices


class GeneSet:
    '''
    Enforces gene organization rules: 
        A genomic location may correspond to many names
        A name may correspond to many locations
        Primary organization should be by genomic location
    '''

    def __init__(self):
        self.genes_by_name = defaultdict(list)
        self.genes_by_chr = OrderedDict()

    def add_genes(self, new_genes):
        for new_gene in new_genes:
            self.add_gene(new_gene)

        return self

    def add_gene(self, new_gene):
        
        #finds if gene location is already inhabited
        if new_gene.get_location() in self.genes_by_chr:
            
            #get existing gene object
            existing_gene = self.genes_by_chr[new_gene.get_location()]

            #adds the new names to the existing object for this genomic location / gene
            existing_gene.add_aliases(new_gene.aliases)
            
            #adds pointers from these names to the existing gene
            for alias in new_gene.aliases:
                #if this alias is not registered
                if not alias in self.genes_by_name:
                    #add the location under the alias
                    self.genes_by_name[alias.upper()].append(existing_gene)

        else:
            
            #add this gene under its genomic location
            self.genes_by_chr[new_gene.get_location()] = new_gene

            #for its names, add this location
            for alias in new_gene.aliases:
                #uppercase the gene name so that capitalization is not a factor 
                self.genes_by_name[alias.upper()].append(new_gene)        

        return self

    def get_symbols(self):
        return [gene.get_name() for gene in self]

    def get_locations(self):
        return [gene.get_location() for gene in self]

    def get_distinct_genes_by_symbol(self, excluding = set()):

        names = self.get_symbols()

        distinct_names = set(names).difference(excluding)

        distinct_genes = GeneSet()
        for dinstinct_name in distinct_names:
            distinct_genes.add_gene(self.get_gene_by_name(dinstinct_name))

        return distinct_genes


    def get_gene_by_name(self, name):

        name = name.upper()
        if name not in self.genes_by_name:
            raise KeyError()
        else:
            return self.genes_by_name[name][0]


    def __str__(self):
        return '\t'.join(['location','gene_name','tad_domain','aliases']) + '\n' + '\n'.join([repr(gene) for gene in self.genes_by_chr.values()])

    def from_str(self, save_str):

        lines = save_str.split('\n')
        #skip header line
        for line in lines[1:]:
            location, name, tad, aliases = [x.strip() for x in line.split('\t')]
            new_gene = Gene(*location.split(':'), aliases.split('|'), tad_domain = tad)
            self.add_gene(new_gene)


    def get_genes_by_chrom(self, chromosome):
        return [
            gene for location_key, gene in self.genes_by_chr.items() if location_key.split(':')[0] == chromosome
        ]


    def __len__(self):
        return len(self.genes_by_chr)

    def __iter__(self):
        return iter(list(self.genes_by_chr.values()))


    def match_user_provided_genes(self, user_genelist):
        rejects = []
        selected_genes = GeneSet()
        for gene_candidate in user_genelist:
            try:
                selected_genes.add_gene( self.get_gene_by_name(gene_candidate) )
            except KeyError:
                rejects.append(gene_candidate)
        
        return selected_genes.get_distinct_genes_by_symbol()


    def random_sample(self, sample_num, seed = None):
        if not seed is None:
            np.random.seed(seed)
            random.seed(seed)

        assert(len(self) > sample_num), 'Background gene list provided must contain more than {} genes'.format(str(sample_num))
        
        #if same number of genes, skip sample
        if len(self) == sample_num:
            return self
        
        else:
            return GeneSet().add_genes(np.random.choice(list(self), sample_num, replace = False))

    #enforce sampling according to the TAD distribution of the genome.
    def sample_by_TAD(self, sample_num, seed = None):

        if not seed is None:
            np.random.seed(seed)
            random.seed(seed)

        #intersect background gene list with TAD_list to eliminate reserved genes
        TAD_data = [(gene, gene.tad_domain) for gene in self]

        #collect list of genes in each bin
        genes_in_TAD = defaultdict(list)
        for gene, tad_group in TAD_data:
            genes_in_TAD[tad_group].append(gene)
        
        #calculates number of genes in each TAD group
        num_genes_in_TAD = {tad_group : len(genes) for tad_group, genes in genes_in_TAD.items()}

        #calculates the number of genes expected to be sampled from each TAD, with a bit of an enrichment (1.1x). Ensure all TADs have >= 1 expected genes
        expected_samples = {
            tad_group : max(1, int(num_genes / len(TAD_data) * sample_num * 1.1))
            for tad_group, num_genes in num_genes_in_TAD.items()
        }

        #samples the expected number of genes from the TAD (or the number of genes if expected is > actual)
        selected_genes = GeneSet()
        for tad_group, num_expected in expected_samples.items():
            sampled_genes = np.random.choice( genes_in_TAD[tad_group], min(num_genes_in_TAD[tad_group] - 1, num_expected), replace = False)
            selected_genes.add_genes(sampled_genes)

        return selected_genes


#creates a label vector, indexed by gene symbol, with binary labels: {0 = background, 1 = query}
def create_label_dictionary(query_list, background_list):
    return { #this dict is used to create label vector for LISA algo
        gene.get_location() : int(i < len(query_list))
        for i, gene in enumerate(list(query_list) + list(background_list))
    }, dict( #return this dict for results (purely informational)
        query_symbols = query_list.get_symbols(),
        background_symbols = background_list.get_symbols(),
        query_locations = query_list.get_locations(),
        background_locations = background_list.get_locations(),
    )
    
#converts user-entered data to two gene lists: a query list, and a background list
def select_genes(query_list, gene_set, num_background_genes = 3000, 
    background_strategy = 'regulatory', max_query_genes = 200, background_list = [], seed = None):
    '''
    query_list: text-processed list of symbols or refseqIDs (or both) that user supplies
    num_background_genes: number of background genes to select
    background_strategy: regulatory, provided, or random background selection strategy
    user_background_genes: text-processed list of symbols or refseqIDs (or both) that user supplies as candidate background genes
    max_query_genes: config setting
    background_list: text-processed list of symbols or refseqIDs (or both)
    '''
                
    query_genes = gene_set.match_user_provided_genes(query_list)

    assert(20 <= len(query_genes) <= max_query_genes), 'User must provide list of 20 to {} unique genes. Provided {}'\
        .format(str(max_query_genes), str(len(query_genes)))

    if background_strategy == 'provided':

        background_genes = gene_set.match_user_provided_genes(background_list)

        assert( len(set(background_genes.get_symbols()).intersection(set(query_genes.get_symbols()))) == 0), 'Some genes in your query are also in your background genes'
        assert( len(background_genes) > len(query_genes) )

    else:

        assert( num_background_genes >= len(query_genes) and len(query_genes) <= 19000//2 ), "More query genes selected than background genes"

        background_candidates = gene_set.get_distinct_genes_by_symbol(excluding = query_genes.get_symbols())

        assert(len(background_candidates) >= num_background_genes), 'Number of background candidates must be greater than or equal number of genes to select as background.'

        #if no down-sampling is needed:
        if len(background_candidates) == num_background_genes:
            background_genes = background_candidates

        elif background_strategy == 'regulatory':

            background_genes = background_candidates.sample_by_TAD(num_background_genes, seed = seed)

            if len(background_genes) > num_background_genes:
                background_genes = background_genes.random_sample(num_background_genes, seed = seed)
            
        elif background_strategy == 'random':
            background_genes = background_candidates.random_sample(num_background_genes, seed = seed)

        else:
            raise AssertionError('Background selection strategy {} not supported'.format(background_strategy))

    #print(type(query_genes), len(query_genes), type(background_genes), len(background_genes))

    return create_label_dictionary(query_genes, background_genes)