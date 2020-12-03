

# Python API

## Usage

The LISA module is implemented in the style of an sklearn estimator. First, the user instantiates a LISA model object, then uses that object to predict TFs from any number of genelists. "lisa.predict" returns a results object with handy data manipulation methods such as sortby and subset. Converting from the results object to a Pandas dataframe for analysis is shown above. Note that pandas is not required to use LISA.

```python
from lisa import LISA

lisa = LISA('hg38', cores = 10)

with open('sox2_down.txt', 'r') as f:
    query_list = f.readlines()

results, metadata = lisa.predict(query_list)
```

To convert the "results" object to a python dictionary, use:

```python
results_dict = results.to_dict()
```

This "results_dict" can also be converted directly to a Pandas DataFrame, as shown below:

```python
import pandas as pd
results_df = pd.DataFrame(results_dict)
```
## API


### lisa.LISA(species, assays = ['Direct','H3K27ac','DNase'], cores = 1, isd_method = 'chipseq', num_datasets_selected = 10, verbose = True, log = None)

### Parameters
<b>species : {'hg38', 'mm10'}</b><br>
Specify from which genes to choose for analysis<br><br>
<b>assays : {'Direct', 'H3K27ac', 'Direct'}, default=['Direct','H3K27ac','DNase']</b><br>
Which assays to use to compute TF influence of your genes of interest. User may pass list of any assays from the options above.<br><br>
<b>cores : int, default=1</b><br>
Number of cores. Recommended to set equal to "num_datasets_selected". In the default case, use at most 10 cores.<br><br>
<b>isd_method : {'chipseq', 'motifs'}, default='chipseq'</b><br>
Use ChIP-seq or Motifs data to assess TF binding positions on the genome.<br><br>
<b>num_datasets_selected : int, default=10</b><br>
Identify N accessibility datasets that show differential accessibility around your genes-of-interest.<br><br>
<b>verbose : bool, int, default=True</b><br>
Verbosity of LISA test status updates to stderr. To reduce, set verbose to <i>int</i> decribing number of levels to print. <i>True</i> sets no limit on verbosity, <i>False</i> and <i>0</i> permit no printing.<br><br>
<b>log : Log object, default=None</b><br><br>

### Methods

### predict(query_list, background_list = [], background_strategy = 'regulatory', num_background_genes = 3000, seed = 2556)<br><br>

<b>Parameters</b><br>
<b>query_list : list</b><br>
List of gene names (symbols or refseq IDs) that make up your genes of interest. Must contain more than 20 genes and less than 500.<br><br>
<b>background_strategy : {'regulatory', 'random', 'provided'}, default='regulatory'</b><br>
How to choose background genes to compare TF influence compared to your genes-of-interest. Regulatory option performs a stratified sampling of TADs and basal accessibility. If "provided", user must specify background genes manually.<br><br>
<b>background_list : list, default=[]</b><br>
If user chooses "provided" background strategy, must provide list of gene names here.<br><br>
<b>num_background_genes : int, default=3000</b><br>
Number of background genes to compare to genes of interest. More background genes increase stability of test, but require longer computation time.<br><br>
<b>seed : int, default=2556</b><br>
Seed controls random sampling of genes and initilization weights of models.<br><br>
<b>Returns</b><br>
<b>results : <i>class</i> LISA_Results</b><br>
A "LISA_Results" object which stores tabular results. For description of fields, see main github page. To convert LISA_Results to pandas DataFrame, use:<br>
<b>df = pandas.DataFrame(results.todict())</b><br><br>
<b>metadata : dict</b><br>
Dictionary of metadata used LISA models, for description of fields, see main github page.<br>