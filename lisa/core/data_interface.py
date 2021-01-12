import h5py as h5
import os
import configparser
import numpy as np
from scipy import sparse

from lisa.core.utils import indices_list_to_sparse_array
from lisa.core import gene_selection
from lisa.core import genome_tools

from lisa._version import __file__ as version_file
from lisa._version import __version__
from lisa.core.utils import Log
from urllib import request, error

PACKAGE_PATH = os.path.dirname(version_file)
CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'h5_config.ini')
INSTALL_PATH = os.path.join(PACKAGE_PATH, 'data')
REQURED_DATASET_VERSION = '2.0'

h5_config = configparser.ConfigParser()
h5_config.read(CONFIG_PATH)

COMPRESSION = 'gzip'

class DatasetNotFoundError(KeyError):
    pass

class DataInterface:

    _config = h5_config
    data_path = os.path.join(PACKAGE_PATH, 'data')

    @classmethod
    def get_window_bedfile_str(cls, species, window_size):

        genome = cls.load_genome(species, window_size)

        window_strs = []
        for i, window in enumerate(genome.list_windows()):
            window_strs.append(str(window) + '\t' + str(i))

        return '\n'.join(window_strs)

    @classmethod
    def get_metadata_headers(cls, technology):
        return cls._config.get('metadata', technology + '_headers').split(',')
    
    @classmethod
    def get_dataset_url(cls, species, window_size):
        return h5_config.get('lisa_params','h5_path')\
            .format(path = cls._config.get('cistrome','data_url'), species = species, 
            version = REQURED_DATASET_VERSION, window = str(window_size))

    @classmethod
    def get_dataset_path(cls, species, window_size):
        return cls._config.get('lisa_params','h5_path').format(path = cls.data_path, species = species,
                version = REQURED_DATASET_VERSION, window = str(window_size))

    #___ DATASET DOWNLOADING ____
    @classmethod
    def fetch_from_cistrome(cls, species, window_size):

        dataset_url = cls.get_dataset_url(species, window_size)

        if not os.path.isdir(cls.data_path):
            os.mkdir(cls.data_path)

        filename, _ = request.urlretrieve(dataset_url)
        os.rename(filename, cls.get_dataset_path(species, window_size))

    @classmethod
    def load_genome(cls, species, window_size):
        return genome_tools.Genome.from_file(cls._config.get('genome','genome')\
            .format(package_path = PACKAGE_PATH, species = species), 
                window_size= window_size)

    def __init__(self, species, window_size = 1000, 
        download_if_not_exists = True, make_new = False, log = None, 
        path = None, load_genes = True):

        self.species = species
        self.window_size = int(window_size)

        if log is None:
            self.log = Log()
        else:
            self.log = log

        if path is None:
            self.path = self.get_dataset_path(self.species, self.window_size)
        else:
            self.path = path

        if make_new:
            h5.File(self.path, 'w').close()
        elif not os.path.isfile(self.path):
            if download_if_not_exists and path is None:
                self.download_data()
            else:
                h5.File(self.path, 'w').close()

        #___ LOAD GENE DATA FROM PACKAGE _____
        self.genome = self.load_genome(self.species, self.window_size)

        if load_genes:
            self.load_genes()

        
    def load_genes(self):
        self.log.append('Loading gene info ...')
        self.genes = gene_selection.GeneSet.from_refseq(self._config.get('genome','genes')\
            .format(package_path = PACKAGE_PATH, species = self.species), self.genome)

        self.gene_loc_set = genome_tools.RegionSet([gene.get_tss_region() for gene in self.genes], self.genome)

        self.rp_map_locs = np.array([r.annotation.get_location() for r in self.gene_loc_set.regions])


    def get_install_path(self):
        return self.data_path

    def get_windows(self):
        return '\n'.join(
            str(r) for r in self.genome.list_windows()
        )

    # ____ RP MAP DATA _____

    @staticmethod
    def _make_basic_rp_map(gene_loc_set, region_set, decay):

        distance_matrix = gene_loc_set.map_intersects(region_set, lambda x,y : x.get_genomic_distance(y), slop_distance=50000)

        distance_matrix.data = np.power(2, -distance_matrix.data/decay)

        return distance_matrix.tocsr()

    def _make_enhanced_rp_map(self, gene_loc_set, region_set, decay):

        #make regions x exons map and exons x genes map
        try:
            indptr, indices, exons = [0],[],[]
            for locus in gene_loc_set.regions:
                new_exons = locus.annotation.get_exon_regions()
                exons.extend(new_exons)
                indices.extend(range(indptr[-1], indptr[-1] + len(new_exons)))
                indptr.append(indptr[-1] + len(new_exons))

            exon_gene_map = sparse.csc_matrix((np.ones(len(exons)), indices, indptr), shape = (len(exons), len(gene_loc_set.regions)))
            
            exons = genome_tools.RegionSet(exons, self.genome)
            region_exon_map = region_set.map_intersects(exons, distance_function = lambda x,y : x.overlaps(y, min_overlap_proportion=0.4),slop_distance=0) #REGIONS X EXONS

            region_exon_map = region_exon_map.dot(exon_gene_map).astype(np.bool)

            not_exon_promoter = 1 - region_exon_map.sum(axis = 1).astype(np.bool)

            basic_rp_map = self._make_basic_rp_map(gene_loc_set, region_set, decay)

            enhanced_rp_map = basic_rp_map.transpose().multiply(not_exon_promoter) + region_exon_map

            return enhanced_rp_map.transpose()

        except Exception as err:
            print(repr(err))
            return region_exon_map, exon_gene_map


    def build_binned_rp_map(self, style, rp_decay):

        region_set = genome_tools.RegionSet(list(self.genome.list_windows()), self.genome)

        if style == 'basic':
            return self._make_basic_rp_map(self.gene_loc_set, region_set, rp_decay)
        elif style == 'enhanced':
            return self._make_enhanced_rp_map(self.gene_loc_set, region_set, rp_decay)
        else:
            NotImplementedError()

    @staticmethod
    def set_attributes(dataset, attr_dict):
        for key, value in attr_dict.items():
            dataset.attrs[key] = value

    def get_rp_map_shape(self):
        return (len(self.genes), len(self.genome))

    def add_rp_map(self, style, rp_map):

        assert(rp_map.shape == self.get_rp_map_shape()), \
            'RP map must be of shape (num genes, num bins): ' + str(self.get_rp_map_shape())

        rp_map_path = self._config.get('rp_map','rp_map').format(style = style)

        rp_map = rp_map.tocsr()

        with h5.File(self.path, 'a') as data:

            if rp_map_path in data:
                del data[rp_map_path]

            group = data.create_group(rp_map_path)

            group.create_dataset('indptr', data = rp_map.indptr, dtype = np.int32, compression=COMPRESSION)
            group.create_dataset('indices', data = rp_map.indices, dtype = np.int32, compression=COMPRESSION)
            group.create_dataset('data', data = rp_map.data, dtype = np.float32, compression=COMPRESSION)

            self.set_attributes(group,dict(shape = rp_map.shape))

    def get_rp_maps(self):
        
        try:
            with h5.File(self.path, 'a') as data:
                return list(data['rp_maps'].keys())
        except KeyError:
            return []

    def get_rp_map(self, style):

        rp_map_path = self._config.get('rp_map','rp_map').format(style = style)

        with h5.File(self.path, 'r') as data:

            try:
                group = data[rp_map_path]

                rp_map = sparse.csr_matrix(
                    (group['data'][...], group['indices'][...], group['indptr'][...]), shape = group.attrs['shape']
                )
            except KeyError:
                raise DatasetNotFoundError(rp_map_path)
        
        return rp_map

    #___ BIN PROJECTION _____

    def check_bin_map_unique(self, bin_map):
        return len(np.unique(bin_map)) == len(bin_map)


    def project_indices(self, indices, bin_map):

        input_hits = sparse.csc_matrix(
            (np.ones_like(indices), indices, [0, len(indices)]),
        )

        input_hits = self.project_sparse_matrix(input_hits, bin_map, None)

        return input_hits.tocoo().row

    @staticmethod
    def project_array(arr, bin_map, num_bins):
        #assert(check_bin_map_unique(bin_map[:,0]) and check_bin_map_unique(bin_map[:,1])), 'To project array, bin_map must have all one-to-one mappings'
        new_arr = np.zeros(num_bins)

        new_arr[bin_map[:,1]] = arr[bin_map[:,0]]

        return new_arr

    @staticmethod
    def project_sparse_matrix(input_hits, bin_map, num_bins, binarize = False):

        index_converted = input_hits.tocsc()[bin_map[:,0], :].tocoo()

        input_hits = sparse.coo_matrix(
            (index_converted.data, (bin_map[index_converted.row, 1], index_converted.col)), 
            shape = (num_bins, input_hits.shape[1]) if not num_bins is None else None 
        ).tocsr()

        if binarize:
            input_hits.data = np.ones_like(input_hits.data)

        return input_hits


    #___ BINDING FACTOR DATA _____
    def get_factor_hit_path(self, technology, dataset_id):
        return self._config.get('factor_binding','hits').format(technology = technology, dataset_id = dataset_id)

    def get_metadata(self, attributes, technology, dataset_id):
        return {dataset_id : {key : attributes[key] for key in self.get_metadata_headers(technology)}}

    def transpose_metadata(self, metadata, technology):

        headers = self.get_metadata_headers(technology)
        sample_ids = list(metadata.keys())

        return {'sample_id' : sample_ids, **{key : [metadata[sample][key] for sample in sample_ids] for key in headers}}

    def add_binding_data(self, technology, dataset_id, hit_bins, **metadata):

        hits_path = self.get_factor_hit_path(technology, dataset_id)

        with h5.File(self.path, 'a') as data:
            if hits_path in data:
                del data[hits_path]

            hits = data.create_dataset(hits_path, data = np.array(hit_bins), dtype = np.int32, compression=COMPRESSION)
            self.set_attributes(hits, metadata)

    def get_binding_dataset(self, technology, dataset_id):

        metadata_headers = self.get_metadata_headers(technology)

        with h5.File(self.path, 'r') as data:
            
            factor_dataset_path = self.get_factor_hit_path(technology, dataset_id)
            
            try:
                hit_bins = np.array(data[factor_dataset_path][...])

                attributes = data[factor_dataset_path].attrs
            except KeyError:
                raise DatasetNotFoundError(factor_dataset_path)

            metadata = self.get_metadata(attributes, technology, dataset_id)
        
        return hit_bins, metadata

    def get_binding_data(self, technology):
        
        with h5.File(self.path, 'r') as data:

            dataset_ids = list(data[self._config.get('factor_binding','root').format(technology=technology)].keys())

            indices = []
            metadata = dict()
            for dataset_id in dataset_ids:
                hit_bins, sample_meta = self.get_binding_dataset(technology, dataset_id)

                metadata.update(sample_meta)
                indices.append(hit_bins)

        hits_matrix = indices_list_to_sparse_array(indices, len(self.genome))

        return hits_matrix.transpose(), np.array(dataset_ids), self.transpose_metadata(metadata, technology)

    def remove_binding_dataset(self, technology, dataset_id):

        factor_dataset_path = self.get_factor_hit_path(technology, dataset_id)

        with h5.File(self.path, 'a') as data:
            del data[factor_dataset_path]

    def list_binding_datasets(self, technology):

        try:
            with h5.File(self.path, 'r') as data:

                dataset_ids = list(data[self._config.get('factor_binding','root').format(technology=technology)].keys())

            return dataset_ids

        except KeyError:
            return []
        

    #____ PROFILE DATA _____
    def add_profile_data(self, technology, dataset_id, profile, rp_maps, rp_map_styles, 
        norm_depth = 1e5, **metadata):
        
        assert(len(rp_maps) == len(rp_map_styles))

        profile_path = self._config.get('profiles','profile').format(technology = technology, dataset_id = dataset_id)

        profile = np.array(profile)
        if len(profile.shape) == 1:
            profile = profile[:,np.newaxis]
        assert(len(profile.shape) == 2)
        assert(profile.shape[0] == self.genome.num_windows_in_genome())
        
        if not norm_depth is None:
            profile = profile/profile.sum() * norm_depth

        with h5.File(self.path, 'a') as data:

            if profile_path in data:
                del data[profile_path]
            
            hits = data.create_dataset(profile_path, data = profile, dtype = np.float16, compression=COMPRESSION)
            self.set_attributes(hits, metadata)

            for rp_map, style in zip(rp_maps, rp_map_styles):

                rp_matrix_path = self._config.get('profiles','rp_matrix_col').format(technology = technology, style = style, dataset_id = dataset_id)

                if rp_matrix_path in data:
                    del data[rp_matrix_path]

                rp_matrix_col = data.create_dataset(rp_matrix_path, data = rp_map.dot(profile), dtype = np.float32, compression=COMPRESSION)
                self.set_attributes(rp_matrix_col, metadata)

    def remove_profile(self, technology, dataset_id):

        profile_path = self._config.get('profiles','profile').format(technology = technology, dataset_id = dataset_id)

        with h5.File(self.path, 'a') as data:
            del data[profile_path]
            
            for style in self.get_rp_maps():
                rp_matrix_col_path = self._config.get('profiles','rp_matrix_col').format(technology = technology, style = style, dataset_id = dataset_id)
                del data[rp_matrix_col_path]


    def get_profile(self, technology, dataset_id):

        profile_path = self._config.get('profiles','profile').format(technology = technology, dataset_id = dataset_id)

        with h5.File(self.path, 'r') as data:

            try:
                profile = np.array(data[profile_path][...])

                attributes = data[profile_path].attrs
            except KeyError:
                raise DatasetNotFoundError(profile_path)

            metadata = self.get_metadata(attributes, technology, dataset_id)

        return profile, metadata

    def list_profiles(self, technology):

        profiles_dir = self._config.get('profiles','root').format(technology = technology)

        try:
            with h5.File(self.path, 'r') as data:
                
                dataset_ids = list(data[profiles_dir].keys())

            return dataset_ids
        except KeyError:
            return []

    def get_rp_matrix(self, technology, style):

        with h5.File(self.path, 'r') as data:
            
            rp_matrix_dir = self._config.get('profiles','rp_matrix').format(technology = technology, style = style)

            dataset_ids = list(data[rp_matrix_dir].keys())

            slices = []
            for _id in dataset_ids:
                slices.append(
                    np.array(data[rp_matrix_dir][_id][...])
                )

        return np.concatenate(slices, axis = 1), np.array(dataset_ids)

    def download_data(self):

        with self.log.section('Grabbing {} data (~15 minutes):'.format(self.species)):
            
            self.log.append('Downloading from database ...')

            try:
                self.fetch_from_cistrome(self.species, self.window_size)
            except error.URLError as err:
                raise AssertionError('ERROR: Cannot connect to cistrome.org for data (usually due to security settings on some servers)!\nView github pages for manual dataset install instructions.')

            self.log.append('Done')