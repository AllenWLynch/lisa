from lisa.core.lisa_core import LISA_Core
from lisa.lisa_public_data.genes_test import FromGenes
from lisa.core.data_interface import PACKAGE_PATH, REQURED_DATASET_VERSION, INSTALL_PATH
from lisa.core.lisa_core import CONFIG_PATH as base_config_path
import numpy as np
from scipy import sparse
from lisa.core import genome_tools
import os
import configparser
from collections.abc import Iterable
from lisa.lisa_user_data.assays import ISD_Assay
from lisa.core.data_interface import DataInterface
from collections import Counter
from lisa.core.utils import Log, LoadingBar
import subprocess
import tempfile

CONFIG_PATH = os.path.join(os.path.dirname(__file__),'config.ini')
_config = configparser.ConfigParser()
_config.read([base_config_path, CONFIG_PATH])


class FromCoverage(LISA_Core):

    window_size = FromGenes.window_size

    @classmethod
    def using_bigwig(cls, species):
        passed

    @classmethod
    def _get_genome_bin_path(cls, species):
        return _config.get('bam_test_params', 'genome_bins').format(data_path = INSTALL_PATH, species = species, window_size = str(cls.window_size))

    @classmethod
    def _write_genome_bins(cls, species):
        
        bedstr = DataInterface.get_window_bedfile_str(species, cls.window_size)

        with open(cls._get_genome_bin_path(species, window_size), 'w') as bed:
            bed.write(bedstr)

    @classmethod
    def convert_bigwig(cls, bigwig, species, bigwig_cmd_path):

        log = Log()

        genome = DataInterface.load_genome(species, cls.window_size)
        coverage_array = np.zeros(len(genome))

        log.append('Converting BigWig file to coverage array ...')

        if not os.path.exists(cls._get_genome_bin_path(species)):
            log.append('Writing bins ...')
            cls._write_genome_bins(species)

        try:

            temp = tempfile.NamedTemporaryFile('w', delete=False)
            temp.close()

            process = subprocess.run([bigwig_cmd_path, bigwig, cls._get_genome_bin_path(species), temp.name], capture_output=True)

            if process.returncode == 0:

                with open(temp.name, 'r') as cmd_output:
                    for line in cmd_output:
                        fields = line.strip().split('\t')
                        coverage_array[int(fields[0])] = fields[4]
                    
                return coverage_array
            
            else:
                raise AssertionError(process.stderr.decode('utf-8'))
        finally:
            os.remove(temp.name)


    def __init__(self, species, coverage_array, rp_map = 'enhanced_10K', isd_method = 'chipseq', verbose = 4, log = None):
        '''
*class*
**lisa.FromRegions** (species, regions, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', verbose = 4, log = None)**

    Initialize the LISA test using user-defined regions.

    Parameters:
        species: {'hg38', 'mm10'}

        coverage_array: (1D or Nx1 np.ndarray):
            Array of scores over 1kb bins.
        sd_method {"chipseq", "motifs"}:
            Use ChIP-seq data or motifs to mark TF binding locations.
        rp_map {"basic_10K", "enhanced_10K"}:
            Choice of RP map, which maps the regulatory influence of a region to a gene. The "basic_10K" model is based simply off distance, with the "enhanced_10K" model masks out the promoter and exon regions of other nearby genes.
        verbose (int):
            Number of levels of log messages to print to stderr
    
    Returns:
        lisa object
        '''

        super().__init__(species, _config, self.window_size, isd_method= isd_method, verbose=verbose, log = log)
        
        assert(isinstance(coverage_array, np.ndarray)), 'Coverage array must be of type numpy.ndarray'
        assert(len(coverage_array.shape) == 1 or (1 in coverage_array.shape and len(coverage_array.shape) == 2)), 'Coverage array must be 1D array or column/row vector'

        coverage_array = coverage_array.reshape(-1)
        assert(len(coverage_array) == len(self.data_interface.genome)), 'Coverage array must be of same length as genome: {} bins'.format(str(len(self.data_interface.genome)))
        self.coverage_array = coverage_array

        rp_map_styles = self._config.get('bam_test_params','rp_map_styles').split(',')
        assert(rp_map in rp_map_styles), 'RP map must be numpy/scipy.sparse array, or be one of provided maps: {}'.format(','.join(rp_map_styles))
        self.rp_map_style = rp_map

    def _load_factor_binding_data(self):
        return self.data_interface.get_binding_data(self.isd_method)
        
    def _load_rp_map(self):
        return self.data_interface.get_rp_map(self.rp_map_style)

    def _initialize_assays(self, **assay_kwargs):
        self.add_assay(
            ISD_Assay(**assay_kwargs, technology = 'BamTest')
        )

    def predict(self, query_list, background_list = [], background_strategy = 'all', num_background_genes = 3000, 
        seed = 2556):
        '''
    *method*
    **.predict** (self, query_list, background_list = [], background_strategy = 'all', num_background_genes = 3000, seed = 2556)
    
        Predict TF influence given a set of genes.
        
        Params:
            query_list (list):
                Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
            background_list (list):
                User-specified list of background genes to compare with query_list. Must contain more genes than query list and entire list will be used. If provided, ```background_strategy``` must be set to "provided".
            background_strategy {"regulatory","random","provided"}:
                Regulatory will sample background genes from a stratified sample of TADs and regulatory states, random will randomly sample from all non-query genes.
            num_background_genes (int):
                Number of genes to use as comparison to query genes. More background genes make test slower, but more stable.
            seed (int):
                Seed for gene selection and regression model initialization.

        Returns:
            results (lisa.core.utils.LISA_Results):
                Can be passed directly to a the pandas constructor: ``results_df = pd.DataFrame(results.to_dict())``.
            metadata (dict):
                Dictionary with test metadata. Includes query genes provided and background genes that were selected. This metadata dict also contains information on the accessibility datasets that were selected to represent the chromatin landscape around you genes-of-interest, for example, the tissue and cell line from which the profiles were derived.
        
        '''
        return super().predict(query_list, region_scores = self.coverage_array, background_list=background_list, background_strategy=background_strategy, 
            num_background_genes= num_background_genes, seed=seed)