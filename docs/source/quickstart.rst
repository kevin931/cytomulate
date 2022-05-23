####################
Quickstart Guide
####################

We get it: you want a TLDR for how to use this package. Here you go. We will walk you from installation
to cytomulating entire datasets.


*****************************************************
Creation Mode: Probabilistic Model-Based Simulation
*****************************************************

By default, the Creation Mode simulates 1 batch of 10 cell types with 20 markers.
No cell differentiation will be simulated.


.. code:: python

    import numpy as np
    from cytomulate.creation.cytof_data import CreationCytofData
    np.random.seed(42)
    # Instantiate a CreationCytofData object
    cytof_data = CreationCytofData()
    # Then, we generate models for each cell type
    cytof_data.initialize_cell_types()
    # Finally, we can generate random samples
    expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)


*****************************************************
Emulation Mode: Real Data-Based Simulation
*****************************************************

To use the Emulation Mode, we need some real cell expressions as well as the corresponding cell type labels.

For the sake of demonstration, we will simulate some cell expressions and cell type labels.
No cell differentiation will be simulated.

.. code:: python

    # For simulating from multivariate Gaussian
    import numpy as np
    from cytomulate.emulation.cytof_data import EmulationCytofData
    np.random.seed(42)
    # We simulate two cell types with different menas
    X1 = np.random.multivariate_normal(mean=np.zeros(20), cov=np.eye(20), size=500)
    X2 = np.random.multivariate_normal(mean=np.ones(20), cov=np.eye(20), size=500)
    expression_matrix = np.vstack((X1, X2))
    # Truncate the expressions to make them CyTOF-like
    expression_matrix = np.clip(expression_matrix, 0, None)
    cell_labels = np.repeat([0,1], 500)
    # Instantiate an EmulationCytofData object
    cytof_data = EmulationCytofData()
    cytof_data.initialize_cell_types(expression_matrix=expression_matrix,
                                     labels=cell_labels)
    expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)



