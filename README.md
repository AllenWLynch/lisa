
# LISA: Landscape In-Silico deletion Analysis

## About

LISA is a statistical test for the influence of Transcription Factors on a set of genes which leverages integrative modeling of chromatin accessiblity and factor binding to make predictions that go beyond simple co-expression analysis. Particularly, LISA models the effects of deleting the influence of a TF on the cis regulatory elements of your genes-of-interest. For more information, see <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6</a>. This implementation extends the original, running faster, reducing dependencies, and adding useful CLI functions for pipeline integration.

## Requirements

* Mac or Linux OS
* Python 3.6+
* 15 GB of available storage space

## Installation

*LISA will install data into the virutal environment's ```site_packages``` directory, so ensure the env's location can store ~15GB.*

### PyPI

It is recommended to install lisa to a virtual environment:

```bash
$ python3 -m venv .venvs/lisa_env
$ source .venvs/lisa_env/bin/activate
```
Install LISA to this virtual env using this command:

```bash
(lisa_env) $ pip install lisa2
```

Or, if you want the newest version from github (not a sanctioned release):

```bash
(lisa_env) $ pip install git+git://github.com/liulab-dfci/lisa2.git#egg=lisa2
```

### Conda

First, create a virtual environment:

```bash
(base) $ conda create --name lisa_env
(base) $ conda activate lisa_env
```

Then install from Conda:

```bash
(lisa_env) $ conda install -c allenwlynch lisa2
```

### Dataset Installation Issues

If you successfully install lisa but the program fails while downloading data, follow these [manual dataset installation instructions.](docs/troubleshooting.md)

## Usage

Installing LISA adds a command to your path:

```bash
(lisa_env) $ lisa 
Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6 X. Shirley Liu Lab, 2020

positional arguments:
  {oneshot,multi,one-vs-rest,download,unpack,backcompat,run-tests}
                        commands
    oneshot             Use LISA to infer genes from one gene list. If you have multiple lists, this option will be slower than using "multi" due to data-loading
                        time.
    multi               Process multiple genelists. This reduces data-loading time if using the same parameters for all lists.
    one-vs-rest         Compare gene lists in a one-vs-rest fashion. Useful downstream of cluster analysis.
    download            Download data from CistromeDB. Use if data recieved is incomplete or malformed.
    unpack              Helper command for manually installing Lisa's data
    backcompat          Interface with LISA using the version 1 command line constructs

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

LISA's functionality is divided into three subcommands. If you have one set of genes-of-interest, run the "oneshot" command. If you have many gene lists from different experiments, run "multi", which treats each list independently but reduces data-loading time compared to running "oneshot" multiple times.

### oneshot usage:

To try out LISA, download a sample gene list from Cistrome.org. This genelist contains 149 genes that are differentially expressed when the transciption factor Sox2 is knocked out. Let's see if LISA can recover the source the change in expression.

```bash
(lisa_env) $ wget http://cistrome.org/~alynch/data/lisa_data/test_genelists/sox2_down.txt
```

Now, run the "lisa oneshot" command, which will be fastest since it only loads data to memory relevant to one genelist. The first time you run this command, it will download LISA's required data from the Cistrome server, which may take about 15 minutes. Once the data is downloaded, LISA will execute in ~30 seconds.
**To make this command run much faster, fill the "-c/--cores" parameter, up to 10**

```bash
(lisa_env) $ lisa oneshot hg38 sox2_down.txt -b 300 -c 1 --seed=2556 --save_metadata> results.tsv
```

The example above shows common a usage pattern of the "oneshot" command using 300 genes as a comparitive background and with a seed supplied so that results are repeatable. You may want to use more background genes (up to 3000) to ensure a more stable prediction. 

The user must also specify the genes' species of origin, in this case human, hg38.
This command prints a table of TFs sorted by regulatory effect on the genes-of-interest, seen here saved to ```results.tsv```. For detailed descriptions of your options, run "lisa oneshot --help".

Quickly inspecting the results, LISA found the effects of SOX2 regulation on this genelist to be highly-ranked! (7th of >8000 ChIP samples). The other TF found to regulate this genelist, NANOG, functions in concert with SOX2 to establish cell identity, so this is a strong prediction as well.  

```bash
$ head -n10 results.tsv | cut -f1,3
Rank 	 factor 
1 	 NANOG 
2 	 NANOG 
3 	 NANOG 
4 	 NANOG 
5 	 NANOG 
6 	 NANOG 
7 	 SOX2 
8 	 NANOG 
1 	 NANOG
```

### multi usage:

To try this command, download a folder of different genelists from Cistrome server, and unpack them. These genelists contain differentially-expressed genes resulting from the knockout and activation of four transcription factors.

```bash
(lisa_env) $ wget http://cistrome.org/~alynch/data/lisa_data/genelists.tar.gz
(lisa_env) $ tar -xvf genelists.tar.gz
```

Now run ```lisa multi```, pointed at the directory of genelists. You many also provide a list of files, but each genelist must have a unique filename, as this filename is used to save the results.

```bash
(lisa_env) $ mkdir results
(lisa_env) $ lisa multi hg38 test_genelists/*.txt -o results/ -c 10 -b 500 --seed=2556
```

The command above independently processes all genes lists in the "genelists" folder, and saves the results tables to the "results" folder. The top factors influencing each gene list are then printed to stdout as a summary table.

## Output File Formats

#### Tabular Results - ChIP-seq mode

| Column Name | Description |
|--|--|
| Rank | Ranking of factor influence with respect to "summary_p_value" |
| sample_id | Cistrome.org sample ID |
| factor | Factor gene symbol |
| cell_line | cell line of ChIP-seq sample | 
| cell_type | cell type |
| tissue | tissue type |
| DNase_p_value | TF influence assessed through DNase accessibility |
| ChIP-seq_p_value | TF influence through direct binding locality enrichment |
| H3K27ac_p_value | TF influence through H3K27ac accessibility |
| DNase_factor_accessibility_z_score | Z-normalized score of accessibility of chromatin around the assessed factor's TSS. This may indicate if a particular factor scores highly for influence, but is not expressed in the accessibility samples used to assess that influence.
| H3K27ac_factor_accessibility_z_score | same as above |
| summary_p_value | Cauchy combined p-value aggregating results of all tests conducted. This is the best indicator of overall factor influence on your genes of interest. |


#### Unstructured results (use --save_metadata option to keep)

Contains various calculations used to conduct LISA's tests. Select keys shown below:
| Key | Value |
|--|--|
| query_symbols | User-provided symbols representing genes-of-interest |
| background_symbols | Background genes provided by user or selected by LISA |
| DNase/H3K27ac -> chromatin_model -> coefs | Weights assigned to the contribution of each accessibility dataset to the LISA test |
| DNase/H3K27ac -> selected_datasets | Selected accessibility datasets' metadata and attributes |

The metadata and weights of accessibility datasets used in the LISA test may be important for performing your analysis of the results, and can show which tissues are highly accessibility around your genes of interest.

## Python Module

Once installed, LISA may also be used as a python module. Documentation can be viewed [here](docs/lisa_base.md).

## Changelog

### [2.1.0] - 2020-12-01

* Bigfixes in output of "lisa multi" test
* Refactored classes for future extension to user-supplied fragment files and peaks
* Added integration testing
* Added factor accessibility introspection to results printout
* Made RP maps substitutable for future tests
* Made assays modular so users can specify which statistical tests they are interested in

### [2.0.6] - 2020-11-22

* Support for Lisa version 1 API for integration with LISA website
* Bugfixes in motif mode results
* Slight speedups in parallelization of insilico-delition computing

## Support

If you have questions, requests, or issues, please email alynch@ds.dfci.harvard.edu.