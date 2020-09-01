
# LISA: Landscape In-Silico deletion Analysis

## About

LISA predicts which TFs regulate a set of genes using integrative modeling of chromatin accessiblity and ChIP-seq binding. Particularly, LISA models the effects of deleting the influence of a TF on the expression of the genes-of-interest. For more information, see <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6</a>. This implementation extends the original, running 30x faster, reducing dependencies, and adding useful CLI functions for pipeline integration. 


## Installation

LISA is available for install from PyPI. Simply set up a virtual environment, then run:
```
pip install LISA
```


## Usage

```bash
>>>LISA
___________________________________________________________________________________________________________________________
Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6 X.
Shirley Liu Lab, 2020 
___________________________________________________________________________________________________________________________

positional arguments:
  {oneshot,multi,one-vs-rest}
                        commands
    oneshot             Use LISA to infer genes from one gene list. If you
                        have multiple lists, this option will be slower than
                        using "multi" due to data-loading time.
    multi               Process multiple genelists. This reduces data-
                        loading time if using the same parameters for all
                        lists.
    one-vs-rest         Compare gene lists in a one-vs-rest fashion. Useful
                        downstream of cluster analysis.

optional arguments:
  -h, --help            show this help message and exit
```

LISA's functionality is divided into three subcommands. If you have one set of genes-of-interest, run the "oneshot" command. If you have many gene lists from different experiments, run "multi", which treats each list independently but reduces data-loading time compared to running "oneshot" multiple times.

Lastly, if you are performing cluster-based differential expression analysis, a common workflow involves identifying lists of cluster-specific differentially-expressed genes. These genelists can be used by LISA in a one-vs-rest fashion to identify TFs that regulate each cluster. 

### oneshot usage:

```bash
>>> LISA oneshot hg38 test_genelist.txt -c 10 > results.tsv
```
The example above shows common usage of the "oneshot" command running with 10 cores, which will find TFs that regulate the genes contained in "test_genelist.txt". The user must also specify the genes' species of origin, in this case human, hg38. 
This command prints a table of TFs sorted by regulatory effect on the genes-of-interest, seen here saved to "results.tsv".

### multi and one-vs-rest usage:

```bash
>>>LISA multi hg38 genelists/*.txt -o results/ -c 10
```
| Sample |  Top Regulatory Factors (p < 0.05) |
|--|--|
| nfe2l2_down.txt | NR3C1, RELA, CEBPB, NFE2L2, FOSL2, NKX2-1, HES2, JUN, FOXA1 |
| nfe2l2_up.txt  | HOXB13, NR3C1, SMARCA4, FOXA2, FOXA1, GATA3 |
| notch1_down.txt | ETS1, MYB, LMO1, TAL1, T, RUNX1, GATA3, TLX1 |
| notch1_up.txt  | TAL1, RBPJ, NOTCH1, MYB, SND1, GATA3, BRD4, TCF12, TCF7L1, RUNX1, ETS1 |
| pou5f1_down.txt  | NANOG, SMAD3, CTNNB1, SOX2, POU5F1 |
| pou5f1_up.txt |  SMAD3, FOXA2, CTNNB1, EOMES, FOXH1, EZH2, FOXA1, SMAD2/3, LHX2, SMAD2, SMAD4, TEAD1, NANOG, YAP1 |
| sox2_down.txt |  NANOG, SOX2, SMAD3, POU5F1, CTNNB1, NIPBL, OTX2 |
| sox2_up.txt  |   EOMES, FOXA1, CTNNB1, NANOG, EZH2, SMAD3, SMAD2, FOXH1, SMAD1, TCF4, JARID2, OTX2 |
| stat3_down.txt | MAZ, YAP1, HOXB13, FOXA2, FOS, CEBPB, IRF1, PTPN11, FOXA1, FOSL2, SNAI2, KDM2B, TRIM28, NANOG |
| stat3_up.txt  |  FOSL2, FOS, STAT3, JUND, JUN, GATA3, FOSL1, TCF12, CEBPB, TEAD1, YAP1 |


The command above independently processes all genes lists in the "genelists" folder, and saves the results tables to the "results" folder. The top factors influencing each gene list are then summarized to stdout as shown in the table above. The CLI for the "multi" and "one-vs-rest" tools are identical, but the manner in which background genes are selected differs.


## Python Module

Once installed, LISA may also be used as a python module:

```python
from LISA import LISA

lisa = LISA('hg38', cores = 10)

with open('test_genelist.txt', 'r') as f:
    query_list = f.readlines()

results, metadata = lisa.predict(query_list)

print(results.subset(range(20)).to_tsv())
```

The LISA module is implemented in the style of an sklearn estimator. First, the user instantiates a LISA model object, then uses that object to predict TFs from any number of genelists. "lisa.predict" returns a ```results``` object with handy data manipulation methods such as  ```sortby``` and ```subset```. To convert the results to a python dictionary, run ```results.todict()```, and to convert the results table to a pandas dataframe, ```pd.DataFrame(results.todict())```. Note that pandas is not required to use LISA.