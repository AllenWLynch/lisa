
from numpy import inf

def parse_deseq_file(path, lfc_cutoff = 1.5, pval_cutoff = 0.05, sep = '\t'):

    up_genes, down_genes = [],[]

    with open(path, 'r') as f:
        
        skipped_header = False
        for i, line in enumerate(f):

            if i > 0 and not (line == '' or line == '\n'):
                fields = line.strip().split('\t')
                assert(len(fields) == 7), 'Correctly formatted DESeq2 file must have 7 fields.'

                try:
                    gene_name, lfc, pval = fields[0], float(fields[2]), float(fields[6])

                    if pval <= pval_cutoff and abs(lfc) >= lfc_cutoff:
                        if lfc > 0:
                            up_genes.append(gene_name)
                        elif lfc < 0:
                            down_genes.append(gene_name)

                except ValueError:
                    pass

    return up_genes, down_genes


def parse_macs_file(path):
    '''
**lisa.parse_macs_file** (path)

Parse regions and region scores from MACS2 .xls output file. Region scores are the number of reads at an associated position (pileup field)

Params:
    path (str):
        path to MACS2 .xls output file
    
Returns:
    region_fields (list):
        list of tuples of (chr, start, end) parsed from MACS file.
    region_scores (list):
        list of scores, read depth at each position

    '''

    with open(path, 'r') as bed:
        lines = bed.readlines()

    region_fields = []
    region_scores = []
    found_header = False
    for line in lines:
        if found_header:
            line = line.split('\t')
            assert(len(line) == 10), 'Correctly-formatted MACS2 .xls files have 10 fields.'
            region_fields.append(line[:3])
            region_scores.append(line[5])
        elif not line[0] == "#" and line[:3] == 'chr':
            line = line.split('\t')
            assert(len(line) == 10), 'Correctly-formatted MACS2 .xls files have 10 fields.'
            found_header = True

    return region_fields, region_scores


def parse_regions_file(path, header = False):
    '''
**lisa.parse_regions_file** (path, header = False)

Parse regions and region scores from bed-like file. Must have columns *chrom, start, end* with optional fourth *score* column.

Params:
    path (str):
        path to MACS2 .xls output file
    header (bool):
        If true, skip first line of bedfile
    
Returns:
    region_fields (list):
        list of tuples of (chr, start, end) parsed from MACS file.
    region_scores (list):
        list of scores, read depth at each position
    '''
    with open(path, 'r') as bed:
        lines = bed.readlines()

    assert(type(header) == bool)

    region_fields, region_scores = [],[]
    num_fields = len(lines[0].split('\t'))
    assert(num_fields in [3,4]), 'Bedfile must have only three columns (chr,start,end) or four (chr,start,end,score)'

    for i, line in enumerate(lines[(1 if header else 0) : ]):
        
        line = line.strip().split('\t')
        if len(line) != num_fields:
            raise AssertionError('Error at line #{}: expected {} fields, got {}'\
                .format(str(i), str(num_fields), str(len(line))))
        
        region_fields.append(line[:3])
        if num_fields == 4:
            region_scores.append(line[-1])
        elif num_fields == 3:
            region_scores.append(1)

    return region_fields, region_scores


def parse_bedfile(path, header = False):
    
    regions = []

    with open(path, 'r') as bed:
        for line in bed.readlines()[(1 if header else 0) : ]:
            fields = line.strip().split('\t')
            assert(len(fields) >= 3)
            chrom, start, end = fields[:3]
            regions.append((chrom, start, end))

    return regions