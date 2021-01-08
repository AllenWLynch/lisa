
from lisa.core.lisa_core import LISA_Core
from lisa.core.data_interface import PACKAGE_PATH, REQURED_DATASET_VERSION
from lisa.core.lisa_core import CONFIG_PATH as base_config_path
from lisa.core import gene_selection
import numpy as np
import multiprocessing
from scipy import sparse
from lisa.core import genome_tools
import os
import configparser
from collections.abc import Iterable
from lisa.lisa_user_data.assays import ISD_Assay
from lisa.core.data_interface import DataInterface
from collections import Counter

CONFIG_PATH = os.path.join(os.path.dirname(__file__),'config.ini')
_config = configparser.ConfigParser()
_config.read([base_config_path, CONFIG_PATH])

class ZeroScoreError(Exception):
    pass

class LISA(LISA_Core):
    '''
lisa.FromRegions
****************

Interface for performing LISA test for TF influence using user-provided regions. The regions may be accompanied with a positive weight or score that
notes the strength of that region for your metric of interest. Often, that metric is ATAC-seq or DNase read depth at that region, but you may provide any 
score as long as it is positive. LISA will rescale your metric automatically to suite the test. Regions should be more than 100 bp wide, but less than 1000 bp 
to ensure specificy for TF and motif hits within the regions. Peaks of ~300 bp are optimal since motifs hits can enrich ~75 to 100 bp away from peak summits. 

For optimal performance, your regions-of-interest should number > 20K and cover roughly the whole genome. If your regions are restricted to a certain chromosome,
You must manually provide background genes that are proximal to your regions.

This test also allows more flexibility to change LISA's function for mapping genomic regions' influence on nearby genes. By default, LISA uses 'Regulatory Potential' 
with a decay of 10000 bp, meaning the regions over a gene's TSS recieve maximum influence, and influence decays by half every 10K bp. This decay rate can be increased to 
allow long-range distal elements more weight, or reduced to prioritize promoter influence. The most powerful extension of this flexibility is the ability to specify a 
custom genes x regions matrix, where every region's influence is mapped to every gene. 

This interface outputs results in the same format as the ``FromGenes`` interface.

Example::

    # Read genelist file
    >>> genes_file = open('./genelist.txt', 'r')
    >>> genes = [x.strip() for x in genes_file.readlines()]
    >>> genes_file.close()
    # Instantiate lisa_regions object. You can pass your regions as a python list of lists, or as the filepath to a bedfile
    >>> lisa_regions = lisa.FromRegions('hg38', './regions.bed', isd_method = 'chipseq')
    # Run the LISA test on your genelist
    >>> results, metadata = lisa_regions.predict(genes, num_background_genes = 501)
    # Print results to stdout
    >>> print(results.to_tsv())
    # Get results as pandas DataFrame
    >>> results_df = pd.DataFrame(results.to_dict())

**For more, see `User Guide <docs/user_guide.rst>`_.**

    '''

    @classmethod
    def get_docs(cls):
        return '\n'.join(x.__doc__ for x in 
        [cls, cls.using_bedfile, cls.using_macs_output, cls.__init__, cls.predict])

    window_size = 100


    @classmethod
    def parse_macs_file(cls, path):

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

    @classmethod
    def parse_bedfile(cls, path, header):
        with open(path, 'r') as bed:
            lines = bed.readlines()

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

    @classmethod
    def using_macs_output(cls, species, xls_path, query_genes, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', 
            background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None):
        '''
        *classmethod*
        **lisa.FromRegions.using_macs_output(species, xls_path, query_genes, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None)**
        
        Use regions defined in MACS .xls file, and take the "pileup" field to be the region's score. 
        All arguments are the same as the "using_bedfile" method, except user must pass "xls_path" as path to MACS2 "{name}.xls" file.
        Header parameter has no effect.
        '''

        region_fields, region_scores = cls.parse_macs_file(xls_path)
        
        lisa = cls(species, region_fields, rp_map = rp_map, rp_decay=rp_decay, isd_method=isd_method, verbose=verbose, log=log)

        return lisa.predict(query_genes, region_scores = region_scores, background_list=background_list, 
            background_strategy=background_strategy, num_background_genes=num_background_genes, seed=seed)

    @classmethod
    def using_bedfile(cls, species, path, query_genes, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', 
            background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None):
        '''
    *classmethod*
    **lisa.FromRegions.using_bedfile(cls, species, path, query_genes, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None)**
    
    Run LISA FromRegions test using a bedfile.

    Params
    ------
    species : {'hg38', 'mm10'}
    path : str
        Path to tab-delineated bedfile with columns: chr start end [score]
        The score column is optional.
    query_genes : list
        Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
    rp_map : {"basic", "enhanced"}, scipy.sparse_matrix
        RP map type, currently supports "basic" and "enhanced". User may also pass their own RP map as scipy.sparse_matrix in the shape (genes x regions)
    rp_decay : float, int
        Decay rate of region influence on gene based on distance from TSS. Increase to prioritize distal regions, decrease to prioritize promoters. Default of 10000 bp is balanced.
    isd_method : {"chipseq", "motifs"} 
        Use ChIP-seq data or motifs to mark TF binding locations.
    background_list : list
        User-specified list of background genes to compare with query_list. Must contain more genes than query list and entire list will be used. If provided, "background_strategy" must be set to "provided".
    background_strategy : {"regulatory","random","provided"}
        Regulatory will sample background genes from a stratified sample of TADs and regulatory states, random will randomly sample from all non-query genes.
    num_background_genes : int
        Number of genes to use as comparison to query genes. More background genes make test slower, but more stable.
    seed : int
        Seed for gene selection and regression model initialization.
    header : bool
        Skip first line of bedfile
    verbose: int
        Number of levels of log messages to print to stderr

    Returns
    -------
    results
        lisa.core.utils.LISA_Results with each key representing a table column, sorted by "summary_p_value" field. The dictionary can be passed directly to a the pandas constructor: ``results_df = pd.DataFrame(results.to_dict())``.
    metadata 
        Dictionary with test metadata. Includes query genes provided and background genes that were selected.
        '''
    
        assert(type(header) == bool)

        region_fields, region_scores = cls.parse_bedfile(path, header)

        lisa = cls(species, region_fields, rp_map = rp_map, rp_decay=rp_decay, isd_method=isd_method, verbose=verbose, log=log)

        return lisa.predict(query_genes, region_scores = region_scores, background_list=background_list, 
            background_strategy=background_strategy, num_background_genes=num_background_genes, seed=seed)


    def __init__(self, species, regions, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', verbose = 4, log = None):
        '''
    **lisa.FromRegions(species, regions, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', verbose = 4, log = None)**
    
    Initialize the LISA test using user-defined regions.

    Params
    ------
    species : {'hg38', 'mm10'} 
    regions : list of lists/tuples with format [('chr', start, end), ... ]
        User-defined regions. 
    rp_map : {"basic", "enhanced"}, scipy.sparse_matrix
        RP map type, currently supports "basic" and "enhanced". User may also pass their own RP map as scipy.sparse_matrix in the shape (genes x regions)
    rp_decay : float, int 
        Decay rate of region influence on gene based on distance from TSS. Increase to prioritize distal regions, decrease to prioritize promoters. Default of 10000 bp is balanced.
    isd_method : {"chipseq", "motifs"} 
        Use ChIP-seq data or motifs to mark TF binding locations.
    verbose: int
        Number of levels of log messages to print to stderr
    
    Returns
    -------
    lisa.lisa_user_data.lisa.LISA object
        '''

        super().__init__(species, _config, 100, isd_method= isd_method, verbose=verbose, log = log)

        if isinstance(rp_map, str):
            rp_map_styles = self._config.get('lisa_params','rp_map_styles').split(',')
            assert(rp_map in rp_map_styles), 'RP map must be numpy/scipy.sparse array, or be one of provided maps: {}'.format(','.join(rp_map_styles))
        else:
            assert( isinstance(rp_map, np.ndarry) or isinstance(rp_map, scipy.sparse)), 'RP map must be either numpy ndarry or scipy.sparse matrix'
        self.rp_map = rp_map

        #self.genome = genome_tools.Genome.from_file(self._config.get('paths','genomes').format(package_path = PACKAGE_PATH, species = self.species), window_size=100)

        assert(isinstance(regions, (list, tuple))), '"regions" parameter must be list of region tuples in format [ (chr,start,end [,score]), (chr,start,end [,score]) ... ] or name of bed file.'
        
        self.log.append('Validation user-provided regions ...')

        self.num_regions_supplied = len(regions)

        regions = self._check_region_specification(regions)

        self.region_set = genome_tools.RegionSet(regions, self.data_interface.genome)
        self.region_score_map = np.array([r.annotation for r in self.region_set.regions])

        assert(isinstance(rp_decay, (int, float)) and rp_decay > 0), 'RP decay parameter must be positive int/float'
        self.rp_decay = rp_decay

        assert(len(regions) >= 1000 and len(regions) < 1000000), 'User must provide atleast 1000 reigons, and less than 1 million.'

    def _check_region_specification(self, regions):

        invalid_chroms = Counter()
        valid_regions = []
        for i, region in enumerate(regions):
            assert(isinstance(region, (tuple, list)) and len(region) == 3), 'Error at region #{}: Each region passed must be in format (string \"chr\",int start, int end'\
                .format(str(i))
            
            try:
                new_region = genome_tools.Region(*region, annotation = i)
                self.data_interface.genome.check_region(new_region)
                valid_regions.append(new_region)
            except ValueError:
                raise AssertionError('Error at region #{}: Could not coerce positions into integers'.format(str(i)))
            except genome_tools.NotInGenomeError as err:
                invalid_chroms[region[0]]+=1
                #raise AssertionError('Error at region #{}: '.format(str(i+1)) + str(err) + '\nOnly main chromosomes (chr[1-22,X,Y] for hg38, and chr[1-19,X,Y] for mm10) are permissible for LISA test.')
            except genome_tools.BadRegionError as err:
                raise AssertionError('Error at region #{}: '.format(str(i+1)) + str(err))
        
        if len(invalid_chroms) > 0 :
            self.log.append('WARNING: {} regions encounted from unknown chromsomes: {}'.format(
                str(sum(invalid_chroms.values())), str(','.join(invalid_chroms.keys()))
            ))

        return valid_regions
    
    def _load_factor_binding_data(self):

        self.factor_binding, self.factor_dataset_ids, self.factor_metadata = self.data_interface.get_binding_data(self.isd_method)

        m2m_region_map = np.array(self.region_set.map_genomic_windows(regions_to_bins=False)).astype(int)

        self.factor_binding = self.data_interface.project_sparse_matrix(self.factor_binding, m2m_region_map, 
            len(self.region_set), binarize=True)

        return self.factor_binding, self.factor_dataset_ids, self.factor_metadata

    def _load_rp_map(self):

        if isinstance(self.rp_map, str):
            if self.rp_map == 'basic':
                self.rp_map = self.data_interface._make_basic_rp_map(self.data_interface.gene_loc_set, self.region_set, self.rp_decay)
            elif self.rp_map == 'enhanced':
                self.rp_map = self.data_interface._make_enhanced_rp_map(self.data_interface.gene_loc_set, self.region_set, self.rp_decay)
            else:
                NotImplementedError()
        else:
            desired_shape = (len(self.rp_map_locs), len(self.regions))
            assert(self.rp_map.shape == desired_shape), 'RP_map must be of shape (n_genes, n_regions): {}'.format(str(desired_shape))

        return self.rp_map

    def _initialize_assays(self, **assay_kwargs):

        self.add_assay(
            ISD_Assay(**assay_kwargs, technology = 'Regions')
        )

    def predict(self, query_genes, region_scores = None, background_list = [], background_strategy = 'all', num_background_genes = 3000, seed = 2556):
        '''
        **predict(query_genes, region_scores = None, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556)**
        
        Predict TF influence given a set of genes.
        
        Params
        ------
        query_genes : list
            Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
        region_scores : list or np.ndarray of shape (len(regions), )
            Region scores/weights. Must be same length as regions. If not passed, all regions will be given score of 1.
        background_list : list
            User-specified list of background genes to compare with query_list. Must contain more genes than query list and entire list will be used. If provided, ```background_strategy``` must be set to "provided".
        background_strategy : {"regulatory","random","provided"}
            Regulatory will sample background genes from a stratified sample of TADs and regulatory states, random will randomly sample from all non-query genes.
        num_background_genes : int
            Number of genes to use as comparison to query genes. More background genes make test slower, but more stable.
        seed : int
            Seed for gene selection and regression model initialization.

        Returns
        -------
        results
            lisa.core.utils.LISA_Results with each key representing a table column, sorted by "summary_p_value" field. The results can be passed directly to a the pandas constructor by calling the "to_dict()" command: ``results_df = pd.DataFrame(results.to_dict())``.
        metadata 
            Dictionary with test metadata. Includes query genes provided and background genes that were selected.
        '''

        if region_scores is None:
            region_scores = np.ones(self.num_regions_supplied)

        else:
            assert(isinstance(region_scores, (np.ndarray, list))), 'Passed scores must be of type list or numpy.ndarray'
            assert(len(region_scores) == self.num_regions_supplied), 'Must provided a score for each region passed. Proved {} regions, and {} scores'\
                .format(str(len(self.region_set)), str(self.num_regions_supplied))
            
            try:
                region_scores = np.array(region_scores).astype(np.float64)
            except ValueError as err:
                raise AssertionError('Region score could not be cast to float: ' + str(err))

            assert(len(region_scores.shape) == 1), 'Region scores must be 1-D list or array.'
            
            assert(np.all(region_scores >= 0)), 'All scores must be non-negative'

            region_scores = region_scores[self.region_score_map]

        return super().predict(query_genes, background_list=background_list, background_strategy=background_strategy, 
            num_background_genes= num_background_genes, seed=seed, region_scores = region_scores)