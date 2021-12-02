'''
*******************************************
LISA: Landscape In-Silico deletion Analysis
*******************************************

LISA is a statistical test that finds the most influential Transcription Factors related to a set of genes. We leverage integrative modeling of a comprehensive dataset 
of 100s of chromatin accessiblity samples and 1000s of ChIP experiments to make predictions. Particularly, LISA models how much the *cis*-regulatory elements around 
a gene are influenced by deleting elements associated with a TF (a process we call *insilico* deletion). For more information, see 
`<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6>`_.

.. contents:: Interfaces

'''

from .lisa_user_data.regions_test import FromRegions
from .lisa_public_data.genes_test import FromGenes
from .lisa_user_data.coverage_test import FromCoverage
from .core.io import parse_regions_file, parse_macs_file
from ._version import __version__