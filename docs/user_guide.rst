Lisa User Guide
===============

This notebook demonstrates the use of LISA’s CLI and Python interfaces
on real single-cell RNA-seq and ATAC-seq data. We’ll start our
investigation into mouse skin cell differentiation using the FromGenes
function, then incorporate our own accessibility data with the
FromRegions function.

.. code:: ipython3

    from lisa import FromGenes, FromRegions

