
*******************************************
LISA: Landscape In-Silico deletion Analysis
*******************************************

LISA is a statistical test that finds the most influential Transcription Factors related to a set of genes. We leverage integrative modeling of a comprehensive dataset 
of 100s of chromatin accessiblity samples and 1000s of ChIP experiments to make predictions. Particularly, LISA models how much the *cis*-regulatory elements around 
a gene are influenced by deleting elements associated with a TF (a process we call *insilico* deletion). For more information, see 
`<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>`_.

.. contents:: Interfaces


lisa.FromGenes
**************

Inputs:
    * Genes of interest

Outputs:
    * Predicted TF influence

Interface for performing LISA test for TF influence using public chromatin accessibility data. Given just a set of genes, LISA will identify a subset of a large database
of public H3K27ac and DNase profiles that represents a model for multiple chromatin states around those genes. LISA then assesses the influence of TF binding 
on your genes-of-interest vs a sampling of background genes through through those representative datasets, and aggregates the effects to produce a final p-value.

This model is useful for integrating accessibility and binding data when you have strictly a list of associated genes (from scRNA-seq, for example). If you have 
genes-of-interest, as well as regions-of-interest, you may use the more specific test provided by ``lisa.FromRegions``.

Example::

    with open('./genelist.txt', 'r') as genes_file:
        genes = [x.strip() for x in genes_file.readlines()]

    results, metadata = lisa.FromGenes('hg38').predict(genes)

    results_df = pd.DataFrame(results.to_dict())

For more, see `user guide <user_guide.md>`_.

    

*class*
**lisa.FromGenes** (species, rp_map = 'enhanced_10K', assays = ['Direct','H3K27ac','DNase'], isd_method = 'chipseq', verbose = True, log = None)

    Initialize the LISA test using public data.

    Params:
        species {'hg38', 'mm10'}:
        assays (list of {"Direct","H3K27ac","DNase"}):
            default is all tests
        isd_method {"chipseq", "motifs"}:
            Use ChIP-seq data or motifs to mark TF binding locations.
        rp_map {"basic_10K", "enhanced_10K"}:
            Choice of RP map, which maps the regulatory influence of a region to a gene. The "basic_10K" model is based simply off distance, with the "enhanced_10K" model masks out the promoter and exon regions of other nearby genes.
        verbose (int):
            Number of levels of log messages to print to stderr

    Returns:
        lisa object
        

    *method*
    **.predict** (self, query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556)
    
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
        
        


lisa.FromRegions
****************

Inputs:
    * Genes of interest
    * Regions (from bulk MACS2 peakcalling or peak-count matrix)
    * Region scores (optional)

Outputs:
    * Predicted TF influence

Interface for performing LISA test for TF influence using user-provided regions and genes-of-interest. The regions may be accompanied with a positive weight or score that
notes the strength of that region for your metric of interest. Often, that metric is ATAC-seq or DNase read depth at that region, but you may provide any 
score as long as it is positive. Regions should be more than 100 bp wide, but less than 1000 bp 
to ensure specificy for TF and motif hits within the regions. 

For optimal performance, your regions-of-interest should number > 20K and cover roughly the whole genome. If your regions are restricted to a certain chromosome,
You must manually provide background genes that are proximal to your regions.

This test also allows more flexibility to change LISA's function for mapping genomic regions' influence on nearby genes. By default, LISA uses 'Regulatory Potential' 
with a decay of 10000 bp, meaning the regions over a gene's TSS recieve maximum influence, and influence decays by half every 10K bp. This decay rate can be increased to 
allow long-range distal elements more weight, or reduced to prioritize promoter influence.

This interface outputs results in the same format as the ``FromGenes`` interface.

Example::

    with open('./genelist.txt', 'r') as genes_file:
        genes = [x.strip() for x in genes_file.readlines()]

    #Using Bedfile
    results, metadata = lisa.FromRegions.using_bedfile('hg38', genes, './path_to_bedfile.bed')

    #Using MACS output
    results, metadata = lisa.FromRegions.using_macs_output('hg38', './path_to_macs.xls', genes)

    #Using Estimator
    regions, scores = lisa.parse_bedfile('./path_to_bedfile.bed', header = False)
    results, metadata = lisa.FromRegions('hg38', regions).predict(genes)

    results_df = pd.DataFrame(results.to_dict())

For more, see `User Guide <user_guide.md>`_.

    

*classmethod*
**lisa.FromRegions.using_bedfile** (species, query_genes, bed_path, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None)**

    Run LISA FromRegions test using a bedfile.

    Parameters:
        species: {'hg38', 'mm10'}

        query_genes (list): 
            Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
        bed_path (str): 
            Path to tab-delineated bedfile with columns: chr start end [score]. The score column is optional.

    Returns:
        results (lisa.core.utils.LISA_Results): 
            With each key representing a table column, sorted by "summary_p_value" field. The dictionary can be passed directly to a the pandas constructor: ``results_df = pd.DataFrame(results.to_dict())``.
        metadata (dict): 
            Test metadata. Includes query genes provided and background genes that were selected, as well as reg-scores for top 100 factors on selected genes.
        

*classmethod*
**lisa.FromRegions.using_macs_output** (species, query_genes, xls_path, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556, header = False, verbose = 4, log = None)

    Use regions defined in MACS .xls file, and take the "pileup" field to be the region's score. 
    All arguments are the same as the "using_bedfile" method, except user must pass "xls_path" as path to MACS2 "{name}.xls" file.
    Header parameter has no effect.
        

*class*
**lisa.FromRegions** (species, regions, rp_map = 'enhanced', rp_decay = 10000, isd_method = 'chipseq', verbose = 4, log = None)**

    Initialize the LISA test using user-defined regions.

    Parameters:
        species: {'hg38', 'mm10'}

        regions (list of lists/tuples with format [('chr', start, end), ... ]):
            User-defined regions. 
        rp_map ({"basic", "enhanced"}, scipy.sparse_matrix):
            RP map type, currently supports "basic" and "enhanced". User may also pass their own RP map as scipy.sparse_matrix in the shape (genes x regions)
        rp_decay (float, int):
            Decay rate of region influence on gene based on distance from TSS. Increase to prioritize distal regions, decrease to prioritize promoters. Default of 10000 bp is balanced.
        isd_method {"chipseq", "motifs"}:
            Use ChIP-seq data or motifs to mark TF binding locations.
        verbose (int):
            Number of levels of log messages to print to stderr
    
    Returns:
        lisa object
        

    *method*
    **.predict** (query_genes, region_scores = None, background_list = [], background_strategy = 'all', num_background_genes = 3000, seed = 2556)**
    
        Predict TF influence given a set of genes.
        
        Params:
            query_genes (list):
                Genes-of-interest, in either Symbol of RefSeqID format. Must provide between 20 to 500 genes.
            region_scores (list or np.ndarray of shape (len(regions), ):
                Region scores/weights. Must be same length as regions. If not passed, all regions will be given score of 1.
            background_list (list):
                User-specified list of background genes to compare with query_list. Must contain more genes than query list and entire list will be used. If provided, ```background_strategy``` must be set to "provided".
            background_strategy {"regulatory","random","provided"}:
                Regulatory will sample background genes from a stratified sample of TADs and regulatory states, random will randomly sample from all non-query genes.
            num_background_genes (int):
                Number of genes to use as comparison to query genes. More background genes make test slower, but more stable.
            seed (int):
                Seed for gene selection and regression model initialization.

        Returns
            results:
                lisa.core.utils.LISA_Results with each key representing a table column, sorted by "summary_p_value" field. The results can be passed directly to a the pandas constructor by calling the "to_dict()" command: ``results_df = pd.DataFrame(results.to_dict())``.
            metadata: 
                Test metadata. Includes query genes provided and background genes that were selected, as well as reg-scores for top 100 factors on selected genes.
        

    *method*
    **.get_rp_map** ()

        Return RP map calculated for your regions, along with the gene and regions metadata.

        Returns:
            rp_map (scipy.sparse_matrix):
                The calculated RP map, (genes x regions)
            gene_metdata (list of tuples):
                List of genes, with columns chr,start,end,symbol
            region_metadata (list of tuples):
                List of supplied regions, sorted. Columns are chr,start,end

        Formatting to Anndata::

            lisa_test = lisa.FromRegions('hg38', regions)
            rp_map, gene_metadata, region_metadata = lisa_test.get_rp_map()

            rp_map_anndata = anndata.AnnData(X = rp_map, 
                obs = pd.DataFrame(gene_metadata, columns = ['chr','start','end','symbol']),
                var = pd.DataFrame(region_metadata, columns = ['chr','start','end'])
            )
        

    *method*
    **.get_binding_matrix** ()

        Returns a binary matrix of factor hits within your regions. Rows are regions, columns are either ChIP-seq samples or motifs. 

        Returns:
            factor_binding_matrix (scipy.sparse_matrix):
                Sparse binary matrix of shape (regions x factor). Is "1" if factor is predicted to bind at region, "0" if not.
            region_metadata (list of tuples):
                List of supplied regions, sorted. Columns are chr,start,end
            factor_metadata (dict of lists):
                Metadata for each factor binding profile.

        Formatting to Anndata::

            lisa_test = lisa.FromRegions('hg38', regions)
            factor_binding, regions, factors = lisa_test.get_binding_matrix()

            binding_anndata = anndata.AnnData(X = factor_binding,
                obs = pd.DataFrame(regions, columns = ['chr','start','end'])
                var = pd.DataFrame(factors)
            )

        

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
        
