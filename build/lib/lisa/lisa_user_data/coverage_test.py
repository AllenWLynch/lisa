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
import pyBigWig as bw

CONFIG_PATH = os.path.join(os.path.dirname(__file__),'config.ini')
_config = configparser.ConfigParser()
_config.read([base_config_path, CONFIG_PATH])


class FromCoverage(LISA_Core):
    '''
lisa.FromCoverage
****************

Inputs:
    * Genes of interest
    * BigWig file, coverage over genome

Outputs:
    * Predicted TF influence

Use LISA to infer TF influence on your geneset using your own coverage data. This test is better suited than the "regions" test when your measure produces wide peaks/areas of influence.
An example of this is H3K27ac data, which correlates with gene expression similarly to accessibility, but produces wide peaks that may span many distinct TF binding locations.

This interface outputs results in the same format as the ``FromGenes`` interface.

Example::

    with open('./genelist.txt', 'r') as genes_file:
        genes = [x.strip() for x in genes_file.readlines()]

    results, metadata = lisa.FromRegions.using_bigwig('hg38', genes, './sample.bigwig')

    results_df = pd.DataFrame(results.to_dict())

For more, see `User Guide <user_guide.md>`_.
    '''

    window_size = FromGenes.window_size

    @classmethod
    def get_docs(cls):
        return '\n'.join(x.__doc__ for x in 
        [cls, cls.using_bigwig, cls.__init__, cls.predict])

    @classmethod 
    def using_bigwig(cls, species, query_genes, bigwig_path, rp_map = 'enhanced_10K', isd_method = 'chipseq', background_list = [], 
        background_strategy = 'all', num_background_genes = 3000, seed = 2556, verbose = 4, log = None):
        '''
*classmethod*
**lisa.FromCoverage.using_bigwig** (species, query_genes, bigwig_path, rp_map = 'basic', rp_decay = 10000, isd_method = 'chipseq', background_list = [], background_strategy = 'all', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None)

    Run LISA FromCoverage test using a bigwig coverage file.

    Parameters:
        species: {'hg38', 'mm10'}

        query_genes (list): 
            Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
        bigwig_path (str): 
            Path to bigwig file

    Returns:
        results (lisa.core.utils.LISA_Results): 
            With each key representing a table column, sorted by "summary_p_value" field. The dictionary can be passed directly to a the pandas constructor: ``results_df = pd.DataFrame(results.to_dict())``.
        metadata (dict): 
            Test metadata. Includes query genes provided and background genes that were selected.
        '''

        if log is None:
            log = Log()
        
        coverage_array = cls.convert_bigwig(bigwig_path, species, log = log)

        return cls(species, coverage_array, rp_map = rp_map, isd_method=isd_method, verbose=verbose, log=log)\
            .predict(query_genes, background_list=background_list, background_strategy=background_strategy, num_background_genes=num_background_genes,
            seed=seed)

    @classmethod
    def convert_bigwig(cls, bigwig, species, log = None):

        if log is None:
            log = Log()

        genome = DataInterface.load_genome(species, cls.window_size)
        coverage_array = np.zeros(len(genome))

        log.append('Converting BigWig file to coverage array ...')

        bar = LoadingBar('Progress', len(genome) // 1000 + 1, cold_start=True)

        try:
            coverage_bw = bw.open(bigwig)

            log.append(bar, update_line=True)

            for i, window in enumerate(genome.list_windows()):
                
                if window.chromosome in coverage_bw.chroms():
                    mean_coverage = coverage_bw.stats(*window.to_tuple())[0]
                    coverage_array[i] = mean_coverage

                if i%1000 == 0:
                    log.append(bar, update_line = True)

            return np.nan_to_num(coverage_array)

        finally:
            coverage_bw.close()

    
    def __init__(self, species, coverage_array, rp_map = 'enhanced_10K', isd_method = 'chipseq', verbose = 4, log = None):
        '''
*class*
**lisa.FromCoverage** (species, regions, rp_map = 'enhanced_10K', rp_decay = 10000, isd_method = 'chipseq', verbose = 4, log = None)

    Initialize the LISA test using user-defined regions.

    Parameters:
        species: {'hg38', 'mm10'}

        coverage_array: (1D or Nx1 np.ndarray):
            Array of scores over 1kb bins.
        isd_method {"chipseq", "motifs"}:
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
                Test metadata. Includes query genes provided and background genes that were selected, as well as reg-scores for top 100 factors on selected genes.
        '''
        return super().predict(query_list, region_scores = self.coverage_array, background_list=background_list, background_strategy=background_strategy, 
            num_background_genes= num_background_genes, seed=seed)