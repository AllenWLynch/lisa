
********
Lisa CLI
********

Installing LISA using pip or conda adds the "lisa" command to your path. LISA's functionality is divided into three main subcommands:

* `lisa oneshot`_ : one genelist
* `lisa multi`_ : multiple genelists
* `lisa regions`_ : one genelist and a list of regions
* `lisa coverage`_ : one genelist and bigwig file of coverage over genome
* `lisa deseq`_ : results from DE analysis using DESeq2

Which are used depending on the evidence you have on hand. 

See the `User Guide <user_guide.rst>`_ for more usage information.
See the `Python API <python_api.rst>`_ for more in-depth description of tests and parameters.

lisa oneshot
------------

You have:

* one genelist

Use LISA to infer influential TFs from one gene list, with background epigenetic landscape modeled using public data. 
If you have multiple lists, this option will be slower than using "multi" due to data-loading time. 

Example::

    $ lisa oneshot hg38 ./genelist.txt -b 501 --seed=2556 --save_metadata > results.tsv

::

    usage: lisa oneshot [-h] [-o OUTPUT_PREFIX]
                                [--background_strategy {regulatory,random,provided,all}]
                                [--background_list BACKGROUND_LIST | -b NUM_BACKGROUND_GENES]
                                [-v VERBOSE]
                                [-a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]]
                                [--rp_map_style {enhanced_10K,basic_10K}]
                                [--seed SEED] [--use_motifs] [--save_metadata]
                                {hg38,mm10} query_list

    positional arguments:
      {hg38,mm10}           Find TFs associated with human (hg38) or mouse (mm10) genes
      query_list            user-supplied gene lists. One gene per line in either symbol or refseqID format

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output file prefix. If left empty, will write results to stdout. (default: None)
      --background_strategy {regulatory,random,provided,all}
                            Background genes selection strategy. LISA samples background genes to compare to user's genes-of-interest from a diverse
                                    regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
                                    (default: regulatory)
      --background_list BACKGROUND_LIST
                            user-supplied list of backgroung genes. Used when --background_strategy flag is set to "provided" (default: None)
      -b NUM_BACKGROUND_GENES, --num_background_genes NUM_BACKGROUND_GENES
                            Number of sampled background genes to compare to user-supplied genes (default: 3000)
      -v VERBOSE, --verbose VERBOSE
      -a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...], --assays {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]
                            Which set of insilico-deletion assays to run. (default: ['Direct', 'H3K27ac', 'DNase'])
      --rp_map_style {enhanced_10K,basic_10K}
                            Which style of rp_map to assess influence of regions on genes. "basic" is stricly distance-based, while "enhanced" masks the exon and promoter regions of nearby genes. (default: enhanced_10K)
      --seed SEED           Random seed for gene selection. Allows for reproducing exact results. (default: 2556)
      --use_motifs          Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database). (default: chipseq)
      --save_metadata       Save json-formatted metadata from processing each gene list. (default: False)


lisa multi
----------

You have:

* multiple genelists

Use LISA to infer influential TFs from multiple lists. This function processes each genelist independently in the same manner as the "oneshot" command, but reduces data loading time. Useful when performing 
the test on up and down-regulated genes from multiple RNA-seq clusters.

Example::

    $ lisa multi hg38 ./genelists/*.txt -b 501 -o ./results/

::

    usage: lisa multi [-h] -o OUTPUT_PREFIX [-v VERBOSE]
                              [-b NUM_BACKGROUND_GENES] [--random_background]
                              [-a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]]
                              [--rp_map_style {enhanced_10K,basic_10K}]
                              [--seed SEED] [--use_motifs] [--save_metadata]
                              {hg38,mm10} query_lists [query_lists ...]

    positional arguments:
      {hg38,mm10}           Find TFs associated with human (hg38) or mouse (mm10) genes
      query_lists           user-supplied gene lists. One gene per line in either symbol or refseqID format

    required arguments:
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output file prefix. (default: None)

    optional arguments:
      -h, --help            show this help message and exit
      -v VERBOSE, --verbose VERBOSE
      -b NUM_BACKGROUND_GENES, --num_background_genes NUM_BACKGROUND_GENES
                            Number of sampled background genes to compare to user-supplied genes. These genes are selection from other gene lists. (default: 3000)
      --random_background   Use random background selection rather than "regulatory" selection. (default: regulatory)
      -a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...], --assays {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]
                            Which set of insilico-deletion assays to run. (default: ['Direct', 'H3K27ac', 'DNase'])
      --rp_map_style {enhanced_10K,basic_10K}
                            Which style of rp_map to assess influence of regions on genes. "basic" is stricly distance-based, while "enhanced" masks the exon and promoter regions of nearby genes. (default: enhanced_10K)
      --seed SEED           Random seed for gene selection. Allows for reproducing exact results. (default: 2556)
      --use_motifs          Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database). (default: chipseq)
      --save_metadata       Save json-formatted metadata from processing each gene list. (default: False)


lisa regions
------------

You have:

* one genelist
* regions (250 - 1000 bp wide) of interest related to that list
* optional: a positive score/weight associated with each region (you may pass zero-weight regions, but they do not affect the test and will be filtered out)

Use LISA to infer TF influence on your geneset, but provide your regions-of-interest rather than building a background epigenetic model using public data. When providing 
your own regions, LISA uses higher resolution, more precise binding data to increase the power of the test. Your regions should be between ~250 and 1000 bp in width, and the 
associated score should be positive. Scores are often read-depth at those regions, but can be any metic you think may influence gene regulation.

Example::

    $ lisa regions -r ./regions.bed -q ./genelist.txt --save_metadata > results.tsv
    $ lisa regions -r ./macs_peaks.xls -q ./genelist.txt --macs_xls > results.tsv

::

    usage: lisa regions -q QUERY_GENES -r REGIONS [--header] [--macs_xls]
                                [--rp_map_style {enhanced,basic}]
                                [--rp_decay RP_DECAY] [-o OUTPUT_PREFIX]
                                [--background_strategy {regulatory,random,provided,all}]
                                [--background_list BACKGROUND_LIST | -b NUM_BACKGROUND_GENES]
                                [-v VERBOSE] [--seed SEED] [--use_motifs]
                                [--save_metadata] [-h]
                                {hg38,mm10}

    positional arguments:
      {hg38,mm10}           Find TFs associated with human (hg38) or mouse (mm10) genes

    required arguments:
      -q QUERY_GENES, --query_genes QUERY_GENES
                            user-supplied gene list. One gene per line in either symbol or refseqID format (default: None)
      -r REGIONS, --regions REGIONS
                            Tad-delineated bed file with columns: chr, start, end[, score]. The score column is optional. If not provided, LISA will assign each region a uniform weight. (default: None)

    optional arguments:
      --header              Bed file has header row as first row. The header row may contain  (default: False)
      --macs_xls            If provided, regions file is a MACS2 .xls output file, and the "pileup" field is taken to be the region score. (default: False)
      --rp_map_style {enhanced,basic}
      --rp_decay RP_DECAY   Distance in base-pairs in which the influence of a region on a gene decays by half. Increase for more weight on distal elements, decrease for more weight on promoter elements. (default: 10000)
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output file prefix. If left empty, will write results to stdout. (default: None)
      --background_strategy {regulatory,random,provided,all}
                            Background genes selection strategy. LISA samples background genes to compare to user's genes-of-interest from a diverse
                                    regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
                                    (default: all)
      --background_list BACKGROUND_LIST
                            user-supplied list of backgroung genes. Used when --background_strategy flag is set to "provided" (default: None)
      -b NUM_BACKGROUND_GENES, --num_background_genes NUM_BACKGROUND_GENES
                            Number of sampled background genes to compare to user-supplied genes (default: 3000)
      -v VERBOSE, --verbose VERBOSE
      --seed SEED           Random seed for gene selection. Allows for reproducing exact results. (default: 2556)
      --use_motifs          Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database). (default: chipseq)
      --save_metadata       Save json-formatted metadata from processing each gene list. (default: False)
      -h, --help


lisa coverage
------------

You have:

* one genelist
* bigwig of coverage over the genome

Use LISA to infer TF influence on your geneset using your own coverage data. This test is better suited than the "regions" test when your measure produces wide peaks/areas of influence.
An example of this is H3K27ac data, which correlates with gene expression similarly to accessibility, but produces wide peaks that may span many distinct TF binding locations.

Example::

    $ lisa coverage -bw ./sample.bigwig -q ./genelist.txt --save_metadata > results.tsv

::

    usage: lisa coverage -q QUERY_GENES -bw BIGWIG_PATH
                                [--rp_map_style {enhanced_10K,basic_10K}]
                                [-o OUTPUT_PREFIX]
                                [--background_strategy {regulatory,random,provided,all}]
                                [--background_list BACKGROUND_LIST | -b NUM_BACKGROUND_GENES]
                                [-v VERBOSE] [--seed SEED] [--use_motifs]
                                [--save_metadata] [-h]
                                {hg38,mm10}

    positional arguments:
      {hg38,mm10}           Find TFs associated with human (hg38) or mouse (mm10) genes

    required arguments:
      -q QUERY_GENES, --query_genes QUERY_GENES
                            user-supplied gene list. One gene per line in either symbol or refseqID format (default: None)
      -bw BIGWIG_PATH, --bigwig_path BIGWIG_PATH
                            Bigwig file describing coverage over the genome. (default: None)

    optional arguments:
      --rp_map_style {enhanced_10K,basic_10K}
                            Which style of rp_map to assess influence of regions on genes. "basic" is stricly distance-based, while "enhanced" masks the exon and promoter regions of nearby genes. (default: enhanced_10K)
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output file prefix. If left empty, will write results to stdout. (default: None)
      --background_strategy {regulatory,random,provided,all}
                            Background genes selection strategy. LISA samples background genes to compare to user's genes-of-interest from a diverse
                                    regulatory background (regulatory - recommended), randomly from all genes (random), or uses a user-provided list (provided).
                                    (default: all)
      --background_list BACKGROUND_LIST
                            user-supplied list of backgroung genes. Used when --background_strategy flag is set to "provided" (default: None)
      -b NUM_BACKGROUND_GENES, --num_background_genes NUM_BACKGROUND_GENES
                            Number of sampled background genes to compare to user-supplied genes (default: 3000)
      -v VERBOSE, --verbose VERBOSE
      --seed SEED           Random seed for gene selection. Allows for reproducing exact results. (default: 2556)
      --use_motifs          Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database). (default: chipseq)
      --save_metadata       Save json-formatted metadata from processing each gene list. (default: False)
      -h, --help

lisa deseq
----------

You have:

* RNA-seq differential expression results from DESeq2

Use LISA to infer influential TFs given differentially expressed genes found using DESeq2. Will seperate up-regulated and down-regulated genes into their own LISA tests.

Example::

    $ lisa deseq hg38 ./deseq_results.tsv -o deseq/ -b 501 --seed=2556 --save_metadata

::

    usage: lisa deseq [-h] [-lfc LFC_CUTOFF] [-p PVAL_CUTOFF] [--sep SEP]
                -o OUTPUT_PREFIX [-v VERBOSE]
                [-b NUM_BACKGROUND_GENES] [--random_background]
                [-a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]]
                [--rp_map_style {basic_10K,enhanced_10K}]
                [--seed SEED] [--use_motifs] [--save_metadata]
                {hg38,mm10} deseq_file

    required arguments:
      -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Output file prefix.

    positional arguments:
      {hg38,mm10}           Find TFs associated with human (hg38) or mouse (mm10) genes
      deseq_file            DEseq differential expression output file. Will be parsed for differentially up and down-regulated genes.

    optional arguments:
      -h, --help            show this help message and exit
      -lfc LFC_CUTOFF, --lfc_cutoff LFC_CUTOFF
                            Log2 fold-change cutoff. For up-regulated genes, must have LFC > cutoff. For down-regulated genes, less than -1 * cutoff. 
                            Default of 1 means genes must be up or down-regulated by a factor of 2 to be included in query. (default: 1)
      -p PVAL_CUTOFF, --pval_cutoff PVAL_CUTOFF
                            Adjusted p-value cutoff. Gene must have pval below cutoff to be a query gene. (default: 0.1)
      --sep SEP             Field separator for DESeq output file. (default: tab)
      -v VERBOSE, --verbose VERBOSE
      -b NUM_BACKGROUND_GENES, --num_background_genes NUM_BACKGROUND_GENES
                            Number of sampled background genes to compare to user-supplied genes. These genes are selection from other gene lists. (default: 3000)
      --random_background   Use random background selection rather than "regulatory" selection. (default: regulatory)
      -a {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...], --assays {Direct,H3K27ac,DNase} [{Direct,H3K27ac,DNase} ...]
                            Which set of insilico-deletion assays to run. (default: ['Direct', 'H3K27ac', 'DNase'])
      --rp_map_style {basic_10K,enhanced_10K}
                            Which style of rp_map to assess influence of regions on genes. "basic" is stricly distance-based, while "enhanced" masks the exon and promoter regions of nearby genes. (default: basic_10K)
      --seed SEED           Random seed for gene selection. Allows for reproducing exact results. (default: 2556)
      --use_motifs          Use motif hits instead of ChIP-seq peaks to represent TF binding (only recommended if TF-of-interest is not represented in ChIP-seq database). (default: chipseq)
      --save_metadata       Save json-formatted metadata from processing each gene list. (default: False)
