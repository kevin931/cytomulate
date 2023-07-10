######################
Installation Guide
######################

You can easily install ``cytomulate`` one with command! Follow the installation guide and
you will be able to cytomutate within minutes!

---------

***********
Conda
***********

You can easily install ``cytomulate`` from ``conda``:

.. code-block::

    conda install -c normalizingflow cytomulate


---------

***********
PyPI
***********

If you prefer ``PyPI``, the installation is easy as well:

.. code-block:: 

    pip install cytomulate

---------

*************
Dependencies
*************

These are the core dependencies that you will need for cytomulate. They should
be automatically installed with ``pip`` or ``conda``, but if there is an issue,
you can elect to install on your own.

* scikit-learn
* numpy
* scipy
* networkx
* matplotlib
* tqdm

--------------

********************************
Optional Dependency: PyCytoData
********************************

If you would like compatibility with ``PyCytoData``, you can install an optional
dependency in the following way:

.. code-block::

    pip install PyCytoData

or via ``conda``:

.. code-block::

    conda install pycytodata -c kevin931 -c bioconda

Now, you can have the option to output your simulation results in a ``PyCytoData`` object.

.. note::

    ``PyCytoData`` requires ``Python>=3.7``, which is stricter than ``cytomulate``.
    If you are still running an older version, please consider upgrading.

.. image:: ../../assets/pycytodata.jpg
   :width: 600
   :alt: PyCytoData Alliance
