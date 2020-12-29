
from lisa.core.lisa_core import LISA_Core, PACKAGE_PATH
from lisa.core.lisa_core import CONFIG_PATH as base_config_path
import numpy as np
import multiprocessing
from scipy import sparse
from lisa.lisa_user_data import genome_tools
import os
import configparser
from collections.abc import Iterable
from lisa.lisa_user_data.assays import ISD_Assay

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

*Example:*

.. code:: python

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

**For more, see `User Guide <docs/user_guide.rst>`_.**

    '''

    @classmethod
    def document(cls):
        return cls('mm10', [('chr1',100,200)]).get_docs()

    def __init__(self, species, regions, region_scores = None, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', **kwargs):
        '''
**lisa.FromRegions(self, species, regions, region_scores = None, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq')**
    Initialize the LISA test using user-defined regions.

    Params
    ------
    species : {'hg38', 'mm10'} 
    regions : list or lists/tuples with format [('chr', start, end[, score]), ... ]
        User-defined regions. The score column is optional and if not provided, all regions will be given same weight. This parameter may also be the filename of a bed file with the same format.
    region_scores : list or np.ndarray of shape (len(regions), ) (optional) 
        Region scores/weights. Must be same length as regions. User may not provide regions with a score column and this parameter at the same time.
    rp_map : str, list, scipy.sparse_matrix, np.ndarray 
        RP map type, currently only supports "basic". User may also pass their own RP map of scipy.sparse_matrix or np.ndarry type in the shape (genes x regions)
    rp_decay : float, int 
        Decay rate of region influence on gene based on distance from TSS. Increase to prioritize distal regions, decrease to prioritize promoters. Default of 10000 bp is balanced.
    isd_method : {"chipseq", "motifs"} 
        Use ChIP-seq data or motifs to mark TF binding locations.
    
    Returns
    -------
    lisa object
        '''

        super().__init__(species, _config, isd_method= isd_method, **kwargs)

        if isinstance(rp_map, str):
            rp_map_styles = self._config.get('lisa_params','rp_map_styles').split(',')
            assert(rp_map in rp_map_styles), 'RP map must be numpy/scipy.sparse array, or be one of provided maps: {}'.format(','.join(rp_map_styles))
        else:
            assert( isinstance(rp_map, np.ndarry) or isinstance(rp_map, scipy.sparse)), 'RP map must be either numpy ndarry or scipy.sparse matrix'
        self.rp_map = rp_map

        self.genome = genome_tools.Genome.from_file(self._config.get('paths','genomes').format(package_path = PACKAGE_PATH, species = self.species), window_size=100)

        if isinstance(regions, str):
            assert(os.path.isfile(regions)), 'File not found: {}. If providing regions as string, must point to existing tab-delineated bedfile'.format(regions)
            with open(regions, 'r') as f:
                regions = [x.strip().split('\t') for x  in f.readlines()]

        else:
            assert(isinstance(regions, (list, tuple))), '"regions" parameter must be list of region tuples in format [ (chr,start,end [,score]), (chr,start,end [,score]) ... ] or name of bed file.'

        regions = self._check_region_specification(regions, region_scores)

        self.region_set = genome_tools.RegionSet(regions, self.genome)

        assert(isinstance(rp_decay, (int, float)) and rp_decay > 0), 'RP decay parameter must be positive int/float'
        self.rp_decay = rp_decay

    @staticmethod
    def _check_region_specification(lines, region_scores):

        def _check_score(score):
            try:
                score = float(score)
                assert(score > 0), 'Region scores must be non-negative'
                if score == 0:
                    raise ZeroScoreError()
                return score
            except ValueError:
                raise AssertionError('Cannot cast score: {} to float.'.format(str(score)))

        regions = []
        if not region_scores is None:
            assert(len(lines) == len(region_scores)), 'Number of regions and number of region scores must be identical.'
            assert(isinstance(region_scores, (np.ndarray, list))), 'Region scores must be passed as list or numpy array'
            region_scores = np.array(region_scores)
            assert(len(region_scores.shape) == 1), 'Region scores must be 1-D array'
            assert(np.all(region_scores >= 0)), 'Region scores must be non-negative'

        num_fields = len(lines[0])
        assert(num_fields in [3,4] if region_scores is None else 3), 'Bedfile must have only three columns (chr, start, end), or four columns (chr, start, end, score)\nUser may not provide a bedfile with scores and fill the region_scores parameter at the same time.'

        for line_no, fields in enumerate(lines):
            assert(num_fields == len(fields)), 'ERROR, region #{}: expected {} fields, got {}'.format(str(line_no), str(num_fields), str(len(fields)))
            
            if num_fields == 3:
                score = 1 if region_scores is None else region_scores[line_no]
            elif num_fields == 4:
                score = fields[-1]

            try:
                regions.append(genome_tools.Region(*fields[:3], annotation = _check_score(score)))
            except ZeroScoreError:
                pass
            except ValueError:
                raise AssertionError('Error at region #{}: Could not coerce positions into integers'.format(str(line_no)))

        assert(len(regions) > 0), 'No regions with positive weights were passed'
        
        return regions
            
    @staticmethod
    def _make_basic_rp_map(gene_loc_set, region_set, decay):

        distance_matrix = gene_loc_set.map_intersects(region_set, lambda x,y : x.get_genomic_distance(y), slop_distance=50000)

        distance_matrix.data = np.power(2, -distance_matrix.data/decay)

        return distance_matrix.tocsr()

    def _make_enhanced_rp_map(gene_loc_set, region_set, decay):

        #make regions x exons map
        exons = []

        for loc_key in self.gene_locs:
            exons.extend(self.all_genes.genes_by_chr[loc_key].get_exon_regions())
        
        exons = genome_tools.RegionSet(exons, self.genome)
        region_exon_map = region_set.map_intersects(exons, slop_distance=0) #REGIONS X EXONS

        #make exons x genes map
        gene_loc_dict = dict(zip(self.gene_locs, enumerate(len(self.gene_locs))))

        rows, cols = [],[]
        for row, exon in enumerate(exons.regions):
            col = gene_loc_dict[exon.annotation]
            rows.append(row)
            cols.append(col)

        exon_gene_map = sparse.csc_matrix((np.ones_like(rows), (rows, cols)), shape = (len(exons), len(self.gene_locs)))

        region_exon_map = region_exon_map.dot(exon_gene_map).astype(np.bool)

        not_exon_promoter = 1 - region_exon_map.sum(axis = 1).todense().astype(np.bool)

        basic_rp_map = self._make_basic_rp_map(gene_loc_set, region_set, decay)

        enhanced_rp_map = basic_rp_map.multiply(not_exon_promoter) + region_exon_map

        return enhanced_rp_map

        
    def _load_factor_binding_data(self):
        self.factor_binding, self.factor_dataset_ids = super()._load_factor_binding_data('100')

        m2m_region_map = np.array(self.region_set.map_genomic_windows()).astype(int)

        index_converted = self.factor_binding.tocsc()[m2m_region_map[:,0], ].tocoo()

        self.factor_binding = sparse.coo_matrix(
            (index_converted.data, (m2m_region_map[index_converted.row, 1], index_converted.col)), 
            shape = (len(self.region_set), len(self.factor_dataset_ids))
        ).tocsr()

        self.factor_binding.data = np.ones_like(self.factor_binding.data)

        return self.factor_binding

    def _load_rp_map(self):

        if isinstance(self.rp_map, str):

            gene_loc_set = genome_tools.RegionSet([genome_tools.Region(*x.strip().split(':')) for x in self.rp_map_locs], self.genome)

            self.rp_map_locs = np.array([":".join([r.chromosome, str(r.start), str(r.end)]) for r in gene_loc_set.regions])

            if self.rp_map == 'basic':
                self.rp_map = self._make_basic_rp_map(gene_loc_set, self.region_set, self.rp_decay)
            elif self.rp_map == 'enhanced':
                self.rp_map = self._make_enhanced_rp_map(gene_loc_set, self.region_set, self.rp_decay)
            else:
                NotImplementedError()
        
        else:
            desired_shape = (len(self.rp_map_locs), len(self.regions))
            assert(self.rp_map.shape == desired_shape), 'RP_map must be of shape (n_genes, n_regions): {}'.format(str(desired_shape))

        return self.rp_map

    def _initialize_assays(self, **assay_kwargs):

        region_scores = np.array([r.annotation for r in self.region_set.regions])

        self.add_assay(
            ISD_Assay(region_scores, **assay_kwargs)
        )

    def predict(self, query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556):
        '''
    **self.predict(self, query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556)**
        Predict TF influence given a set of genes.
        
        Params
        ------
        query_list : list
            Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
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
            Dictionary with each key representing a table column, sorted by "summary_p_value" field. The dictionary can be passed directly to a the pandas constructor: ``results_df = pd.DataFrame(results.todict())``.
        metadata 
            Dictionary with test metadata. Includes query genes provided and background genes that were selected.
        '''

        return super().predict(query_list, background_list=background_list, background_strategy=background_strategy, 
            num_background_genes= num_background_genes, seed=seed)