
# LISA: Landscape In-Silico deletion Analysis

## About

LISA predicts which TFs regulate a set of genes using integrative modeling of chromatin accessiblity and ChIP-seq binding. Particularly, LISA models the effects of deleting the influence of a TF on the expression of the genes-of-interest. For more information, see <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6</a>. This implementation extends the original, running 30x faster, reducing dependencies, and adding useful CLI functions for pipeline integration. 

## Requirements

* Mac or Linux OS
* Python 3.6+
* 15 GB of available storage space

## Installation

LISA is available for install from github. It is recommended to install lisa to a virtual environment:

```bash
>>> python3 -m venv .venvs/lisa_env
>>> source .venvs/lisa_env/bin/activate
```
Install LISA to this virtual env using this command:

```bash
(lisa_env) >>> pip install --upgrade git+git://github.com/AllenWLynch/lisa.git#egg=lisa
```

## Usage

Installing LISA via pip adds a command to your path:

```bash
(lisa_env) >>> lisa 
Lisa: inferring transcriptional regulators through integrative modeling of
public chromatin accessibility and ChIP-seq data
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6 X.
Shirley Liu Lab, 2020

positional arguments:
  {oneshot,multi,one-vs-rest,download}
                        commands
    oneshot             Use LISA to infer genes from one gene list. If you
                        have multiple lists, this option will be slower than
                        using "multi" due to data-loading time.
    multi               Process multiple genelists. This reduces data-loading
                        time if using the same parameters for all lists.
    one-vs-rest         Compare gene lists in a one-vs-rest fashion. Useful
                        downstream of cluster analysis.
    download            Download data from CistromeDB. Use if data recieved is
                        incomplete or malformed.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
```

LISA's functionality is divided into three subcommands. If you have one set of genes-of-interest, run the "oneshot" command. If you have many gene lists from different experiments, run "multi", which treats each list independently but reduces data-loading time compared to running "oneshot" multiple times.

Lastly, if you are performing cluster-based differential expression analysis, a common workflow involves identifying lists of cluster-specific differentially-expressed genes. These genelists can be used by LISA in a one-vs-rest fashion to identify TFs that regulate each cluster. 

### oneshot usage:

To try out LISA, download a sample gene list from Cistrome.org. This genelist contains 149 genes that are differentially expressed when the transciption factor Sox2 is knocked out. Let's see if LISA can recover the source the change in expression.

```bash
(lisa_env) >>> wget http://cistrome.org/~alynch/data/lisa_data/test_genelists/sox2_down.txt
--  http://cistrome.org/~alynch/data/lisa_data/test_genelists/sox2_down.txt
Resolving cistrome.org (cistrome.org)... 155.52.218.90
Connecting to cistrome.org (cistrome.org)|155.52.218.90|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 1030 (1.0K) [text/plain]
Saving to: ‘sox2_down.txt’

sox2_down.txt                          100%[===========================================================================>]   1.01K  --.-KB/s    in 0s      
2020-09-04 09:59:08 (75.6 MB/s) - ‘sox2_down.txt’ saved [1030/1030]
```

Now, run the "lisa oneshot" command, which will be fastest since it only loads data to memory relevant to one genelist. The first time you run this command, it will download LISA's required data from the Cistrome server, which may take about 15 minutes. Once the data is downloaded, LISA will execute in ~30 seconds.

```bash
(lisa_env) >>> lisa oneshot hg38 sox2_down.txt -c 10 -b 3000 --seed=2556 > results.tsv
Data not found, must download from CistromeDB ...
Grabbing hg38 data (~15 minutes):
	Downloading from database ...
	Extracting data ...
        Done!                                   
        
Using 10 cores ...
Matching genes and selecting background ...
	Selected 133 query genes and 500 background genes.
Loading data into memory (only on first prediction):
	Loading binding data ...
	Loading regulatory potential map ...
	Loading ChIP-seq RP matrix ...
	Loading DNase RP matrix ...
	Loading H3K27ac RP matrix ...
	Done!
Calculating ChIP-seq peak-RP p-values ...
Modeling DNase purturbations:
	Selecting discriminative datasets and training chromatin model ...
	Calculating in-silico deletions:
		Reading DNase data: [====================]
		Performing in-silico knockouts ...
		Calculating Δ regulatory score ...
		Calculating p-values ...
	Done!
Modeling H3K27ac purturbations:
	Selecting discriminative datasets and training chromatin model ...
	Calculating in-silico deletions:
		Reading H3K27ac data: [====================]
		Performing in-silico knockouts ...
		Calculating Δ regulatory score ...
		Calculating p-values ...
	Done!
Mixing effects using Cauchy combination ...
Formatting output ...
Done!
```

The example above shows common a usage pattern of the "oneshot" command, running with 10 cores, using 3000 genes as a comparitive background, and with a seed supplied so that results are repeatable.  The user must also specify the genes' species of origin, in this case human, hg38. 
This command prints a table of TFs sorted by regulatory effect on the genes-of-interest, seen here saved to "results.tsv".

```bash
(lisa_env) >>> cat results.tsv | cut -f1,3,7-8 | head -n10
Rank	factor	combined_p_value	    combined_p_value_adjusted
1	    NANOG	  0.0	                  0.0
2	    NANOG	  0.0	                  0.0
3	    NANOG	  5.551115123125783e-17	3.836930773104541e-13
4	    NANOG	  5.551115123125783e-17	3.836930773104541e-13
5	    NANOG	  5.551115123125783e-17	3.836930773104541e-13
6	    SOX2	  5.773159728050814e-15	3.9904080040287226e-11
7	    SMAD3	  9.325873406851315e-15	6.446043698815629e-11
8	    NANOG	  9.126033262418787e-14	6.307914190983865e-10
9	    SMAD3	  9.614531393253856e-14	6.645564099017065e-10
```

LISA found the effects of SOX2 regulation on this genelist to be statistically significant (p = 3.9904080040287226e-11 << 0.01)! The other TF found to regulate this genelist, NANOG, functions in concert with SOX2 to establish cell identity, so this is a strong prediction as well.  

### multi and one-vs-rest usage:

To try this command, download a folder of different genelists from Cistrome server, and unpack them. These genelists contain differentially-expressed genes resulting from the knockout and activation of four transcription factors.

```bash
(lisa_env) >>> wget http://cistrome.org/~alynch/data/lisa_data/genelists.tar.gz
(lisa_env) >>> tar -xvf genelists.tar.gz
```

Now run "lisa multi", pointed at the directory of genelists. You many also provide a list of files, but each genelist must have a unique filename, as this filename is used to save the results.

```bash
(lisa_env) >>> mkdir results
(lisa_env) >>>lisa multi hg38 test_genelists/*.txt -o results/ -c 10 -b 500 --seed=2556
```

The command above independently processes all genes lists in the "genelists" folder, and saves the results tables to the "results" folder. The top factors influencing each gene list are then summarized to stdout as shown in the table above. 

The CLI for "lisa multi" and "lisa one-vs-rest" commands is similar, but the difference is how these commands select background genes to compare with your genes-of-interest. In "lisa multi", background genes are chosen randomly or sampled from different regulatory states (default), but with "lisa one-vs-rest", Each list is compared against all genes in the other lists. This may provide a more robust analysis if each genelist were derived from diffentially-expressed genes in an upstream clustering analysis. 

## Python Module

Once installed, LISA may also be used as a python module:

```python
from lisa import LISA
import pandas as pd

lisa = LISA('hg38', cores = 10)

with open('sox2_down.txt', 'r') as f:
    query_list = f.readlines()

results, metadata = lisa.predict(query_list)

results_df = pd.DataFrame(results.todict())
```

The LISA module is implemented in the style of an sklearn estimator. First, the user instantiates a LISA model object, then uses that object to predict TFs from any number of genelists. "lisa.predict" returns a ```results``` object with handy data manipulation methods such as  ```sortby``` and ```subset```. Converting from the results object to a Pandas dataframe for analysis is shown above. Note that pandas is not required to use LISA.

## Support

If you have questions, requests, or issues, please email alynch@ds.dfci.harvard.edu.
