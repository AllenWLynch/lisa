
from sys import stderr
from contextlib import contextmanager
from scipy import sparse
import numpy as np
from collections import defaultdict, OrderedDict

def ragged_array_to_sparse_matrix(indices, values, col_length):
    return sparse.hstack([
        sparse.csc_matrix(
            (val_column,
            (index_column, np.zeros(len(val_column)))),
            shape = (col_length, 1)
        )
        for index_column, val_column in zip(indices, values)
    ])

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
            
            #did not create new gene
            return False

        else:
            
            #add this gene under its genomic location
            self.genes_by_chr[new_gene.get_location()] = new_gene

            #for its names, add this location
            for alias in new_gene.aliases:
                #uppercase the gene name so that capitalization is not a factor 
                self.genes_by_name[alias.upper()].append(new_gene)        

            return True

    def get_symbols(self):
        return [gene.get_name() for gene in self]

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


class LoadingBar:
    
    def __init__(self, label, increments, length = 25, cold_start = False):
        self.increments = increments
        self.length = length
        self.label = label
        self.progress = 0
        self.cold_start = cold_start
        
    def __str__(self):
        if self.cold_start:
            self.cold_start = False
        else:
            self.increment()
        completed_steps = int(self.progress / self.increments * self.length)
        if completed_steps >= self.length:
            return '{}: [{}]'.format(self.label, "="*completed_steps) + '\n' if self.is_finished() else ''
        else:
            return '{}: [{}>{}]'.format(self.label, "="*completed_steps, " "*(self.length - completed_steps - 1))
    
    def increment(self):
        if not self.is_finished():
            self.progress += 1
        
    def is_finished(self):
        return self.progress >= self.increments


class Log:

    def __init__(self, target = stderr):
        self.target = target
        self.indents = 0

    @contextmanager
    def section(self, header):
        try:
            self.start_section(header)
            yield self
        finally:
            self.end_section()

    def start_section(self, section_header):
        self.append(section_header)
        self.indents += 1

    def end_section(self):
        self.indents -= 1 

    def append(self, text, end = '\n', update_line = False):
        linestart = '\r' if update_line else ''
        print(linestart + '\t'*self.indents + str(text), 
            end = '' if update_line else end, 
            file = self.target)