

from setuptools import setup, find_packages

setup(
    name = 'lisa2',
    description = "Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data. X. Shirley Liu Lab, 2020",
    version = '2.2.5',
    url = 'https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6',
    author = 'Allen Lynch',
    author_email = 'alynch@ds.dfci.harvard.edu',
    packages = find_packages(),
    zip_safe = True,
    scripts = ['bin/lisa'],
    install_requires = [
        'numpy>=1.19,<2',
        'scipy>=1.5,<2',
        'h5py>=2.10.0,<3',
        'scikit-learn>=0.23.2,<1',
        'pyBigWig>=0.3.17,<1'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    python_requires='>=3.6,<4',
    long_description_content_type='text/markdown',
    long_description = '''
# About

LISA is a statistical test for the influence of Transcription Factors on a set of genes. We leverage integrative modeling of public chromatin accessiblity and factor binding to 
make predictions that go beyond simple co-expression analysis. The minimum you need to run LISA is a list of genes-of-interest, 
but you can also supply your own epigenetic background. For more information, see Qin et al., 2020. 
This implementation extends the original, running faster, reducing dependencies, and adding useful CLI functions for pipeline integration.

# Documentation

Please see Lisa's [github repo](https://github.com/liulab-dfci/lisa2) for source code and tutorials.

# Installation

### Requirements

* Mac or Linux OS
* Python 3.6+
* 15 GB of available storage space

### Installation

LISA will install data into the virutal environment's "site_packages" directory, so ensure the env's location can store ~10GB.

### PyPI

It is recommended to install lisa to a virtual environment:

    $ python3 -m venv .venvs/lisa_env
    $ source .venvs/lisa_env/bin/activate

Install LISA to this virtual env using this command:

    (lisa_env) $ pip install lisa2

### Conda

First, create a virtual environment:

    (base) $ conda create --name lisa_env
    (base) $ conda activate lisa_env

Then install from Conda:

    (lisa_env) $ conda install -c liulab-dfci lisa2

### Dependencies

* numpy
* scipy
* h5py
* pyBigWig
* scikit-learn

    ''',
)