#############################################################
Cytomulate: Accurate and Efficient Simulation of CyTOF data
#############################################################

Cytomulate is a package to simulation realistic data for Mass Cytometry or Cytometry by Time-of-Flight (CyTOF).
We support both model-based through **Creation Mode** and real-data-based simulation through **Emulation Mode**.
Cytomulate serves as solutions to benchmarking, method validation, prototyping, and more. You can easily generate
realistic and accurately CyTOF simulations within seconds.

Using our unified framework compatible with `PyCytoData <https://pycytodata.readthedocs.io/en/latest/>`_,
you can stay within our ecosystem: starting from simulation to precessing, you can then select the best dimension
reduction method for your dataset and proceed onto other downstream analyses. All you need is a few simple commands!
For more details, **read our tutorials and documentations linked below!** Or try this example:

.. code-block:: python

   >>> from cytomulate import CreationCytofData

   >>> cytof_data = CreationCytofData()
   >>> cytof_data.initialize_cell_types()
   >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)

**When in doubt, Cytomulate it!**

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 1
   :caption: Tutorial

   tutorial/emulation
   tutorial/creation
   tutorial/cli
   tutorial/complex
   tutorial/visualization
   tutorial/pycytodata
   tutorial/benchmark

.. toctree::
   :maxdepth: 1
   :caption: Technical Details


.. toctree::
   :maxdepth: 1
   :caption: Development

   change/contribution
   change/build
   change/development
   change/index
   license

.. toctree::
   :maxdepth: 1
   :caption: Full API Reference

   documentation/index

.. toctree::
   :maxdepth: 1
   :caption: Resources
   :hidden:

   references
   Our Paper <https://doi.org/10.1101/2022.04.26.489549>
   Dr. Xinlei (Shery) Wang <https://people.smu.edu/swang/>
   Dr. Tao Wang <https://qbrc.swmed.edu/labs/wanglab/aboutpi.php>
   DBAI <https://dbai.biohpc.swmed.edu/>
   GitHub <https://github.com/kevin931/cytomulate/>


***********************
Resources
***********************

For our paper, read our preprint `here <https://doi.org/10.1101/2022.04.26.489549>`_.

For more resources on our labs, collaborators, and related projects, please visit the following:

   * `Dr. Xinlei (Shery) Wang faculty page <https://people.smu.edu/swang/>`_
   * `Dr. Tao Wang Lab <https://qbrc.swmed.edu/labs/wanglab/aboutpi.php>`_
   * `Database for Actionable Immunology (DBAI) for more computational immunology-related tools <https://dbai.biohpc.swmed.edu/>`_

For project development and newest updates, consider visiting the following:

   * Our `Github <https://github.com/kevin931/cytomulate/>`_