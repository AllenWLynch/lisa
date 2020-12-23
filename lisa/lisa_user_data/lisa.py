
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
"""
LISA_Core implements the main methods for results formatting and data loading that make
up a standard LISA application. Extensions of the core method may instantiate different assays
depending on the data available to asses TF influence in different ways.

If only expression data is available, the base LISA interface can be instantiated, which 
will use public data to estimate TF influence.
"""
CONFIG_PATH = os.path.join(os.path.dirname(__file__),'config.ini')
_config = configparser.ConfigParser()
_config.read([base_config_path, CONFIG_PATH])

class LISA(LISA_Core):

    def __init__(self, species, regions, region_scores = None, config = _config, rp_map = 'basic', rp_decay = 10000, **kwargs):
        super().__init__(species, config, **kwargs)

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
                lines = [x.strip().split('\t') for x  in f.readlines()]

            regions = self.check_region_specification(lines, region_scores)

        elif not isinstance(regions, (list, tuple)):

            raise AssertionError('"regions" parameter must be list of region tuples in format [ (chr,start,end [,score]), (chr,start,end [,score]) ... ] or name of bed file.')

            regions = self.check_region_specification(regions, region_scores)

        regions = self.convert_to_region_objs(regions)

        self.region_set = genome_tools.RegionSet(regions, self.genome)

        assert(isinstance(rp_decay, (int, float)) and rp_decay > 0), 'RP decay parameter must be positive int/float'
        self.rp_decay = rp_decay

    @staticmethod
    def check_score(score):
        try:
            score = float(score)
            assert(score >= 0), 'Region scores must be non-negative'
            return score
        except ValueError:
            raise AssertionError('Cannot cast score: {} to float.'.format(str(score)))
    
    @staticmethod
    def check_region_specification(lines, region_scores):

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

            regions.append(fields[:3], self.check_score(score))

    @staticmethod
    def convert_to_region_objs(region_data):
        region_objs = []
        for i, region in enumerate(region_data):
            assert(isinstance(region, (list, tuple))), 'Error at region #{}: Region must be a list or tuple containing (chr, start ,end)'.format(str(i))
            assert(len(region) == 3), 'Error at region #{}: Regions must be of format (chr, start, end)'.format(str(i))
            try:
                region_objs.append(genome_tools.Region(*region[:3], annotation=region[-1]))
            except ValueError:
                raise AssertionError('Error at region #{}: Could not coerce positions into integers'.format(str(i)))
        return region_objs

    @staticmethod
    def _make_basic_rp_map(distance_matrix, decay):
        distance_matrix.data = np.power(2, -distance_matrix.data/decay)
        return distance_matrix
        
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

            inter_region_distances = gene_loc_set.distance_intersect(self.region_set, lambda x,y : x.get_genomic_distance(y), max_distance=50000)

            if self.rp_map == 'basic':
                self.rp_map = self._make_basic_rp_map(inter_region_distances, self.rp_decay).tocsr()
            else:
                NotImplementedError()

            self.rp_map_locs = np.array([":".join([r.chromosome, str(r.start), str(r.end)]) for r in gene_loc_set.regions])
        
        else:
            desired_shape = (len(self.rp_map_locs), len(self.regions))
            assert(self.rp_map.shape == desired_shape), 'RP_map must be of shape (n_genes, n_regions): {}'.format(str(desired_shape))

        return self.rp_map

    def _initialize_assays(self, **assay_kwargs):

        region_scores = np.array([r.annotation for r in self.region_set.regions])

        self.add_assay(
            ISD_Assay(region_scores, **assay_kwargs)
        )