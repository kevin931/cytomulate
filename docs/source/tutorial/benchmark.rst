######################
Benchmark Analyses
######################

Since Cytomulate is open source, we also would like to share our pipelines for obtaining the benchmark
results in `our paper <https://doi.org/10.1101/2022.06.14.496200>`_. Here, we write a brief tutorial
on how to use the codes and how to benefit from it!

--------------------------------

************************
Downloading Source Codes
************************

We released our source codes as a GitHub release, which is linked `here <https://github.com/kevin931/cytomulate/releases/tag/benchmark.rev.1>`_.
All you will have to do is the following:

1. Download the ``benchmark.zip`` under the "Assets" tab.
2. Decompress the folder with the software available for your OS and access the contents.

While we provided documentation in the form of comments in the files themselves, here we
would like to point out a few more things that may be helpful in your Cytomulate journey.

Dependencies
---------------

To run the codes, you will need the following softwares installed on your system:

- Python with Cytomulate (You can download the source codes on the same page or use Cytomulate v0.2.0 release if you prefer.)
- R installation with all the packages listed in the ``library()`` calls.

Datasets
----------

You will also have to download the necessary datasets used in our paper. All the accession methods and
their availablility is in the XXX section of our paper!


--------------------------

**************************
Codes and Functionalities
**************************

In this section, we explain each part of the code (sorted in directories) and what they do in accordance
with our paper's analyses.

.. note::

    The ``FileIO`` class and the ``KLdivergence`` function are included in multiple Python scripts for
    convenience purposes only. In reality, it is okay to write a separate module to house these! In fact,
    the ``FileIO`` class is now part of PyCytoData, which makes life easier.


Directory: batch
------------------

This directory contains codes used to benchmark batch correction methods as shown in the
**Comparing Batch Normalization Methods using Cytomulate** section. The ``data_gen.py``
generates datasets with multiple batches using the complex simulation functionalities
of Cytomulate; then, ``batch_correction.R`` performs batch correction. Finally,
``benchmark.py`` computes the benchmarks as shown in paper.


Directory: clustering
----------------------

This directory contains codes used to benchmark clustering methods as shown in the 
**Validating Clustering Performance using Cytomulate** The overall structure is similar
to that of batch correction codes with the exception of ``clustering.R`` which performs
clustering rather than batch correction and ``benchmark.py`` which uses a different
metric.

Directory: masking
--------------------

These codes are used to randomly mask cells in the Levine_32dim dataset to assess the performance
of all of the methods. The ``data_gen.py`` generates the masked cell types and datasets, whereas
the ``compute_kl.py`` computes the KL benchmark each method as presented in Fig. 5 of our paper.
These results are presented in the **Cytomulate is robust against cell-type resolution** secion.

Directory: metric_computation
------------------------------

This directory contains codes to compute the main metrics used in our paper. The three
main metrics on mean, covariance, and KL divergence are computed in python. The pMSE metric
from ``synthpop`` is computed in R.

Directory: processing_time
---------------------------

The R script in this directory benchmarks the processing times of Cytomulate competitors in R.
Since Cytomulate is the only Python method in our paper, it is not included in the R script.
Rather, if the benchmarking of Cytomulate is desired, you can use the ``/usr/bin/time`` in linux
for timing Cytomulate's CLI or use the various timing modules in python.

The results of processing time and Cytomulate's efficiency are included in the **Cytomulate is efficient**
section of our paper.

Directory: simulations
------------------------

This directory contains codes to generate the simulation results for each method. Each script
is named according to the method. Note that only Cytomulate is in Python, while all others
are using R. These results are used throughout our paper.
