*******************************************
LISA: Landscape In-Silico deletion Analysis
*******************************************

.. contents:: Table of Contents

About
-----

LISA is a statistical test for the influence of Transcription Factors on a set of genes which leverages integrative modeling of chromatin accessiblity and factor binding to make predictions that go beyond simple co-expression analysis. Particularly, LISA models the effects of deleting the influence of a TF on the cis regulatory elements of your genes-of-interest. For more information, see <a href=https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6</a>. This implementation extends the original, running faster, reducing dependencies, and adding useful CLI functions for pipeline integration.

Requirements
------------

* Mac or Linux OS
* Python 3.6+
* 15 GB of available storage space

Installation
------------

*LISA will install data into the virutal environment's ```site_packages``` directory, so ensure the env's location can store ~15GB.*

PyPI
~~~~

It is recommended to install lisa to a virtual environment:

.. code-block:: bash

  $ python3 -m venv .venvs/lisa_env
  $ source .venvs/lisa_env/bin/activate
  
Install LISA to this virtual env using this command:

.. code-block:: bash

  (lisa_env) $ pip install lisa2

Conda
~~~~~

First, create a virtual environment:

.. code-block:: bash

  (base) $ conda create --name lisa_env
  (base) $ conda activate lisa_env

Then install from Conda:

.. code-block:: bash

  (lisa_env) $ conda install -c liulab-dfci lisa2

Dataset Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you successfully install lisa but the program fails while downloading data, follow these `manual dataset installation instructions <docs/troubleshooting.md>`_.

Usage
-----

Command Line Interface
~~~~~~~~~~~~~~~~~~~~~~

LISA's cli offers convenient methods for the most common use cases. See the `API <docs/cli.rst>`_, or try:

.. code-block::

  (lisa_env) $ lisa {command} --help

for parameter descriptions. See the `User Guide <docs/user_guide.rst>`_ for best practices.

Python Interface
~~~~~~~~~~~~~~~~

The python module allows more control over the LISA test and more convenient data analysis. See the `Python API <docs/python_api.rst>`_ and the `User Guide <docs/user_guide.rst>`_.

Output File Formats
-------------------

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

Changelog
---------

**[2.1.0] - 2020-12-01**

* Bigfixes in output of "lisa multi" test
* Refactored classes for future extension to user-supplied fragment files and peaks
* Added integration testing
* Added factor accessibility introspection to results printout
* Made RP maps substitutable for future tests
* Made assays modular so users can specify which statistical tests they are interested in

**[2.0.6] - 2020-11-22**

* Support for Lisa version 1 API for integration with LISA website
* Bugfixes in motif mode results
* Slight speedups in parallelization of insilico-delition computing

Support
-------

If you have questions, requests, or issues, please email alynch@ds.dfci.harvard.edu.
