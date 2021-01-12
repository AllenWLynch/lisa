
import numpy as np
import os
import random
from collections import Counter, defaultdict, OrderedDict
import random
from lisa.core.genome_tools import Region


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
            self.add_alias(names.upper())
        else:
            self.add_aliases([name.upper() for name in names])
        self.location = ':'.join([str(property_) for property_ in (self.chrom, self.start, self.end)])
        
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
        return self.location

    def get_tss_region(self):
        return Region(self.chrom, self.start, self.end, annotation = self)
    
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


class RefSeqGene(Gene):

    promoter_width_from_tss = 1500

    def __init__(self, name, chrom, strand,
        txStart, txEnd, exonStarts, exonEnds, symbol, tad_cluster, *, genome):

        if strand == '-':
            tss = (int(txEnd), int(txEnd) + 1)
        else:
            tss = (int(txStart), int(txStart) + 1)
            
        super().__init__(chrom, *tss, [symbol.upper(), name.upper()], tad_domain=tad_cluster)

        self.special_regions = dict()

        self.add_region(Region(chrom, *tss).slop(self.promoter_width_from_tss, genome))
        if exonStarts != '' and exonEnds != '':
            for exon_start, exon_end in zip(exonStarts.strip(',').split(','), exonEnds.strip(',').split(',')):
                self.add_region(Region(chrom, exon_start, exon_end))

        self.is_noncoding = name[:3] == "NR_"

    def get_exon_regions(self):
        return list(self.special_regions.values())

    def add_region(self, region):
        if not region.to_tuple() in self.special_regions:
            self.special_regions[region.to_tuple()] = region

    def add_regions(self, new_regions):
        for new_region in new_regions:
            self.add_region(new_region)


class GeneSet:
    '''
    Enforces gene organization rules: 
        A genomic location may correspond to many names
        A name may correspond to many locations
        Primary organization should be by genomic location
    '''

    @classmethod
    def from_file(cls, path):

        new_geneset = cls()

        with open(path, 'r') as f:
            for line in f.readlines()[1:]:
                location, name, tad, aliases = [x.strip() for x in line.split('\t')]
                new_gene = Gene(*location.split(':'), aliases.split('|'), tad_domain = tad)
                new_geneset.add_gene(new_gene)

        return new_geneset

    @classmethod
    def from_refseq(cls, path, genome):
        
        new_geneset = cls()
        with open(path, 'r') as f:
            for line in f.readlines():
                #print(line.strip().split('\t'))
                new_geneset.add_gene(RefSeqGene(*[x.strip() for x in line.strip().split('\t')], genome = genome))

        return new_geneset
            
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

            try:
                existing_gene.add_regions(new_gene.special_regions)
            except AttributeError:
                pass

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

    def get_distinct_genes_by_symbol(self, excluding = set(), exclude_nc_rna = True):

        names = self.get_symbols()

        distinct_names = set(names).difference(excluding)

        distinct_genes = GeneSet()
        for dinstinct_name in distinct_names:
            add_gene = self.get_gene_by_name(dinstinct_name)
            try:
                if not add_gene.is_noncoding:
                    distinct_genes.add_gene(add_gene)
            except AttributeError:
                distinct_genes.add_gene(add_gene)

        return distinct_genes


    def get_gene_by_name(self, name):

        name = name.upper()
        if name not in self.genes_by_name:
            raise KeyError()
        else:
            return self.genes_by_name[name][0]


    def __str__(self):
        return '\t'.join(['location','gene_name','tad_domain','aliases']) + '\n' + '\n'.join([repr(gene) for gene in self.genes_by_chr.values()])


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
        
        assert(len(self) > sample_num), 'Background gene list provided must contain more than {} genes'.format(str(sample_num))
        
        #if same number of genes, skip sample
        if len(self) == sample_num:
            return self
        
        else:
            return GeneSet().add_genes(np.random.choice(sorted(list(self), key = lambda x : x.location), sample_num, replace = False))

    #enforce sampling according to the TAD distribution of the genome.
    def sample_by_TAD(self, sample_num, seed = None):
        
        if not seed is None:
            np.random.seed(seed)

        #intersect background gene list with TAD_list to eliminate reserved genes, sorting to maintain order for seed repeatability
        TAD_data = [(gene, gene.tad_domain) for gene in sorted(self, key = lambda x : x.location)]

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
            sampled_genes = np.random.choice(sorted(genes_in_TAD[tad_group], key = lambda x : x.location), min(num_genes_in_TAD[tad_group] - 1, num_expected), replace = False)
            selected_genes.add_genes(sampled_genes)

        return selected_genes