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

*LISA will install data into the virutal environment's ```site_packages``` directory, so ensure the env's location can store ~10GB.*

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
