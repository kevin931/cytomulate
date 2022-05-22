####################
Quickstart Guide
####################

We get it: you want a TLDR for how to use this package. Here you go. We will walk you from installation
to cytomulating entire datasets.


*****************************************************
Creation Mode: Probabilistic Model-Based Simulation
*****************************************************

By default, the Creation Mode simulates 1 batch of 10 cell types with 20 markers.
:

.. code:: python

    from cytomulate.creation.cytof_data import CreationCytofData
    # Instantiate a CreationCytofData object
    cytof_data = CreationCytofData()
    # Then, we generate models for each cell type
    cytof_data.initialize_cell_types()
    # Finally, we can generate random samples
    expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)


*****************************************************
Emulation Mode: Real Data-Based Simulation
*****************************************************

Sample code:

.. code:: python

    import numpy as np
    from cytomulate.emulation.cytof_data import EmulationCytofData