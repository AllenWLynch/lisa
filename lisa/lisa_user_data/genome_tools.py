
import numpy as np
from collections import Counter
from scipy import sparse

class BadRegionError(Exception):
    pass

class Region:

    def __init__(self, chromosome, start, end, annotation = None):
        self.chromosome = str(chromosome)
        self.start = int(start)
        self.end = int(end)
        self.center = int((self.end - self.start)/ 2) + self.start
        self.annotation = annotation

        assert(self.end > self.start), 'End position must be greater than start position'

    @staticmethod
    def _overlap_distance(min1, max1, min2, max2):
        return max(0, min(max1, max2) - max(min1, min2))

    def __len__(self):
        return self.end - self.start

    def overlaps(self, region, min_overlap_proportion = 0):

        assert(isinstance(region, Region))

        if self.chromosome == region.chromosome:
            overlap_dist = self._overlap_distance(self.start, self.end, region.start, region.end)
            return overlap_dist > 0 and overlap_dist >= len(self) * min_overlap_proportion
        else:
            return False

    def get_center(self):
        return self.center

    def get_genomic_distance(self, region):
        if self.chromosome == region.chromosome:
            return np.abs(self.get_center() - region.get_center())
        else:
            return np.inf

    def slop(self, d, genome):
        return Region(self.chromosome, max(0, self.center - d), min(self.center + d, genome.get_chromlen(self.chromosome)))

    def __str__(self):
        return '({},{},{})'.format(self.chromosome, self.start, self.end)

    def __eq__(self, region):
        return self.chromosome == region.chromosome and self.start == region.start and self.end == region.end

class Genome:

    @classmethod
    def from_file(cls, path, window_size = 100):
        with open(path, 'r') as f:
            chromosomes, lengths = list(zip(*[l.strip().split('\t') for l in f.readlines()]))

        return cls(chromosomes, lengths, window_size= window_size)

    def __init__(self, chromosomes, lengths, window_size = 100):
        
        self.window_size = int(window_size)
        self.chromosomes, self.lengths = chromosomes, lengths

        self.lengths = list(map(int, self.lengths))

        self.chrom_idx = dict(zip(self.chromosomes, range(len(self.chromosomes))))

        self.indptr = [0]
        for l in self.lengths:
            self.indptr.append(self.indptr[-1] + self.get_num_windows(l, self.window_size))

    def get_chromlen(self, chromosome):

        assert(isinstance(chromosome, str))
        try:
            return self.lengths[self.chrom_idx[chromosome]]
        except KeyError:
            raise AssertionError(chromosome + ' not in genome.')

    @staticmethod
    def get_num_windows(chrom_len, window_size):
        return int(chrom_len / window_size) + int(chrom_len % window_size > 0)

    def check_region(self, region):
        try:
            assert(region.chromosome in self.chrom_idx), 'Chr {} from region {} not in genome'.format(str(region.chromosome), str(region))
            assert(region.start >= 0),'Region may not start at negative index: {}'.format(str(region))
            assert(region.end <= self.get_chromlen(region.chromosome)), 'Region {} extends past end of chromosome'.format(str(region), region.chromosome)
        except AssertionError as err:
            raise BadRegionError(str(err))

    def get_window_from_position(self, chromosome, position):

        assert(isinstance(chromosome, str))
        self.check_region(Region(chromosome, position, position + 1))

        chr_start_idx = int(position / self.window_size)
        window_idx = self.indptr[self.chrom_idx[chromosome]] + chr_start_idx

        return Region(chromosome, chr_start_idx * self.window_size, min( (chr_start_idx + 1) * self.window_size, self.get_chromlen(chromosome) )), window_idx

    def get_next_window(self, region):
        assert(isinstance(region, Region))
        return self.get_window_from_position(region.chromosome, region.end)

    def num_windows_in_genome(self):
        return self.indptr[-1]

    def get_region(self, window_idx):
        assert(isinstance(window_idx, int) and window_idx >= 0)
        assert(window_idx < self.indptr[-1]), 'Invalid window idx, genome does not contain that many windows'
        for i, num_windows in enumerate(self.indptr[1:]):
            if num_windows > window_idx:
                window_idx -= self.indptr[i]
                start_pos = window_idx * self.window_size
                return self.get_window_from_position(self.chromosomes[i], start_pos)

        raise Exception('An error occured while finding window')
        

class Endpoint:

    def __init__(self,point_val, endpoint_side, *, segment_obj, list_idx, list_id):
        assert(list_id in ['i','j'])
        assert(endpoint_side in ['L','R'])
        self.point_val = point_val
        self.list_id = list_id
        self.endpoint_side = endpoint_side
        self.segment_obj = segment_obj
        self.list_idx = list_idx

def get_endpoints(regions, slop_distance, genome, list_id):
    endpoints = []
    for i, region in enumerate(regions):
        slop = region.slop(slop_distance, genome)
        endpoint_args = dict(
            segment_obj = region,
            list_idx = i,
            list_id = list_id,
        )
        endpoints.append(Endpoint(slop.start, 'L', **endpoint_args))
        endpoints.append(Endpoint(slop.end,'R', **endpoint_args))
    return endpoints

def get_pairs(i_endpoints, j_endpoints, distance_function):
    interactions = []
    all_endpoints = sorted([*i_endpoints, *j_endpoints], key = lambda x : (x.point_val, x.endpoint_side))
    active_i_segments = {}
    active_j_segments = {}
    for endpoint in all_endpoints:
        if endpoint.endpoint_side == 'L':
            if endpoint.list_id == 'i':
                interactions.extend([
                    (endpoint.list_idx, segment.list_idx, distance_function(endpoint.segment_obj, segment.segment_obj))
                    for segment in active_j_segments.values()
                ])
                active_i_segments[endpoint.list_idx] = endpoint
            else:
                interactions.extend([
                    (segment.list_idx, endpoint.list_idx, distance_function(endpoint.segment_obj, segment.segment_obj))
                    for segment in active_i_segments.values()
                ])
                active_j_segments[endpoint.list_idx] = endpoint
        else:
            if endpoint.list_id == 'i':
                del active_i_segments[endpoint.list_idx]
            else:
                del active_j_segments[endpoint.list_idx]
    
    return interactions

class RegionSet:
    #stores set of genomic regions stored for quick chrom access
    def __init__(self, regions, genome):

        for region in regions:
            genome.check_region(region)

        self.genome = genome
        self.regions = sorted(regions, key = lambda r : (r.chromosome, r.start))

        self.chrom_counts = Counter([r.chromosome for r in self.regions])
        self.chrom_order = sorted(self.chrom_counts.keys())

        indptr = [0]
        self.chrom_indptr = {}
        for chrom in self.chrom_order:
            indptr.append(indptr[-1] + self.chrom_counts[chrom])
            self.chrom_indptr[chrom] = (indptr[-2], indptr[-1])

    def get_regions_by_chrom(self, chrom):
        start_idx, end_idx = self.chrom_indptr[chrom]
        return self.regions[start_idx : end_idx]

    def map_intersects(self, regions, distance_function = lambda x, y : 1, slop_distance = 1e5):

        assert(isinstance(max_distance, (float, int)))
        assert(isinstance(regions, RegionSet))
        
        matrix_rows = []
        for chrom1 in self.chrom_order:
            sparse_subsections = []
            for chrom2 in regions.chrom_order:
                submatrix_shape = (self.chrom_counts[chrom1], regions.chrom_counts[chrom2])
                if chrom1 == chrom2:
                    i_endpoints = get_endpoints(self.get_regions_by_chrom(chrom1), slop_distance, self.genome, 'i')
                    j_endpoints = get_endpoints(regions.get_regions_by_chrom(chrom1), slop_distance, self.genome, 'j')
                    interaction_pairs = get_pairs(i_endpoints, j_endpoints, distance_function)
                    i,j,data = list(zip(*interaction_pairs))
                    sparse_subsections.append(sparse.csr_matrix((data, (i,j)), shape = submatrix_shape))
                else:
                    sparse_subsections.append(sparse.csr_matrix(submatrix_shape))

            matrix_rows.append(sparse.hstack(sparse_subsections))

        return sparse.vstack(matrix_rows)

    def __str__(self):
        return '\n'.join([str(r) for r in self.regions])

    def __len__(self):
        return len(self.regions)

    def map_genomic_windows(self, min_window_overlap_proportion = 0.4):

        m2m_map = []

        for i, region in enumerate(self.regions):

            min_start_loc = max(region.start - (1 - min_window_overlap_proportion) * self.genome.window_size, 0)
            max_start_loc = max(region.end - min_window_overlap_proportion * self.genome.window_size, 0)

            window, j = self.genome.get_window_from_position(region.chromosome, min_start_loc)

            while window.start <= max_start_loc:

                if window.overlaps(region, min_overlap_proportion =  min_window_overlap_proportion):
                    m2m_map.append((j, i))
                
                try:
                    window, j = self.genome.get_next_window(window)
                except BadRegionError:
                    break

        return m2m_map

