*******************************************
LISA: Landscape In-Silico deletion Analysis
*******************************************

.. contents:: Table of Contents

About
-----

LISA is a statistical test for the influence of Transcription Factors on a set of genes. We leverage integrative modeling of public chromatin accessiblity and factor binding to make predictions that go beyond simple co-expression analysis. 
The minimum you need to run LISA is a list of genes-of-interest, but you can also supply your own epigenetic background. For more information, see `Qin et al., 2020 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>`_. 
This implementation extends the original, running faster, reducing dependencies, and adding useful CLI functions for pipeline integration.

The Model
-----

The key components of the LISA test are the:
  1. profile, a distribution of accessibility over regions in the genome, supplied by user or predicted from public data
  2. hits, the regions where a TF is predicted to bind (through ChIP-seq or motif)
  3. region-gene map, maps the influence of a region to nearby genes.

First, LISA constructs a null model of gene influence, which assumes each accessible region is occupied by its associated factors, and that all factor-bound regions exert influence on nearby genes. 
LISA then tests for the influence of a factor on a gene by calculating what proportion of that gene's influence could be attributed to that factor binding nearby regions.
When you provide genes-of-interest, LISA finds factors that preferentially affects these genes over a sampling of background genes.

.. image:: docs/model_diagram.png
  :width: 300

See the `User Guide <docs/user_guide.md>`_ to see it in action. 

Requirements
------------

* Mac or Linux OS
* Python 3.6+
* 15 GB of available storage space

Installation
------------

**LISA will install data into the virutal environment's "site_packages" directory, so ensure the env's location can store ~10GB.**

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

for parameter descriptions. See the `User Guide <docs/user_guide.md>`_ for best practices.

Python Interface
~~~~~~~~~~~~~~~~

The python module allows more control over the LISA test and more convenient data analysis. See the `Python API <docs/python_api.rst>`_ and the `User Guide <docs/user_guide.md>`_.

Changelog
---------

**[2.2.0] - 2021-01-10**

Added
~~~~~

* Added "FromRegions" test, and moved all older functionalities to "FromGenes". New feature allows user to run LISA test with their own regions-of-interest
* Added "query_reg_score" and "background_reg_score" matrices to output metadata of "FromRegions" test, which allows user to see which genes are likely regulated by each factor.
* New backend interface for faster file transfers
* Added ability to append more data to backend for future updates, including ATAC-seq epigenetic backgrounds
* Added more documentation and user guide
* Appended new ATAC data and reprocessed motifs using JASPAR database

Removed
~~~~~~~

* Removed "cores" option from multi and oneshot tests, and removed mutliprocessing from package. 
* Removed "one-vs-rest" test because proved to provide unstable results

**[2.1.0] - 2020-12-01**

* Bugfixes in output of "lisa multi" test
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
