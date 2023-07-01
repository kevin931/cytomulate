####################
Quickstart Guide
####################

``Cytomulate`` has an incredibly simple API for our two simulation modes:

- Creation Mode
- Emulation

In this tutorial, we will walk you through the basics of each mode along
with some easy tips and tricks that may be beneficial to you as you start
cytomulating away in your everyday life.

Before we get started, let's set a seed so that we can get reproducible
results:

.. code:: python

    import numpy as np
    np.random.seed(42)

Now, let's start our journey!

--------------------------------

*****************************************************
Creation Mode: Probabilistic Model-Based Simulation
*****************************************************

Creation mode, as its name implies, allows Cytomulate to formulate its
own recipe to create datasets. In other words, we don't mimic any existing
datasets. Rather, you can specify all the parameters and get exactly that
in return (think of it as customizing your frozen yogurt instead of having
a set menu item). 

.. code:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)

By default, the Creation Mode simulates 1 batch of 10 cell types with 20 markers
and your specified number of cells. Now, let's look at the results we get:

.. code:: python

    >>> expression_matrices
    {0: array([[0.0174032 , 0.56522606, 0.30617146, ..., 1.45473617, 0.31720944,
                0.56082489],
               [1.51403051, 0.58787876, 1.47135531, ..., 1.46453964, 1.47289462,
                0.56293968],
               [1.50258188, 0.30403334, 0.55254835, ..., 0.55964366, 0.11598993,
                0.11538903],
               ...,
               [0.02217756, 0.58209763, 0.29106511, ..., 1.43931079, 0.31428325,
                0.55701247],
               [0.52981196, 1.46558442, 0.56604271, ..., 1.46975045, 0.1134796 ,
                0.10429606],
               [1.47349014, 0.62232074, 1.46480141, ..., 1.45502123, 1.4689097 ,
                0.56239449]])}
    >>> labels
    {0: array([2, 0, 7, 4, 2, 6, 1, 4, 0, 8, 5, 2, 0, 2, 7, 1, 2, 5, 0, 2, 4, 0,
               0, 2, 0, 1, 7, 2, 2, 2, 3, 2, 1, 2, 5, 2, 7, 7, 2, 3, 5, 0, 9, 2,
               8, 2, 2, 2, 4, 2, 0, 2, 0, 4, 7, 7, 0, 5, 4, 2, 2, 4, 4, 7, 7, 3,
               0, 9, 6, 0, 7, 2, 9, 2, 7, 5, 0, 1, 0, 8, 1, 2, 7, 2, 4, 0, 2, 0,
               6, 2, 2, 4, 0, 5, 2, 0, 0, 2, 5, 0])}

As some of your keen-eyed readers may notice, the outputs are dictionaries:
this is because Cytomulate can accomodate multiple samples as indexed by
the dictionary keys. Of course, you can procceed to extract the expression
matrix and then work with it in downstream analyses.

PyCytoData Output
------------------------

For those of you who are familiar with ``PyCytoData`` or want a cleaner interface
to work with, Cytomulate can output a ``PyCytoData`` object.

.. note::
    
    ``PyCytoData`` is required to be installed for this to work. Since it is an
    optional dependency, read our `Installation Guide <https://cytomulate.readthedocs.io/en/dev/installation.html>`_
    for further details. Once installed, ``PyCytoData`` is fully compatible with ``Cytomulate``


To do this, simply use the following method instead:

.. code:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)
    
Now, let's look at the object and access the expression matrix and labels: 

.. code:: python

    >>> type(dataset)
    >>> dataset.expression_matrix
    array([[0.03390264, 0.03323944, 0.79319831, ..., 0.00431289, 1.56704157,
            0.11522665],
           [0.00213084, 0.32081423, 0.04375508, ..., 0.03236736, 1.59080603,
            0.03952084],
           [0.03625106, 0.03741527, 0.80485429, ..., 0.00405477, 1.5738564 ,
            0.07236162],
           ...,
           [0.78996499, 0.80564232, 0.03399493, ..., 0.06597879, 0.03527863,
            0.31189172],
           [0.        , 0.32236194, 0.05363561, ..., 0.02794309, 1.58739998,
            0.0293298 ],
           [0.03699177, 0.04021649, 0.80394265, ..., 0.00309843, 1.57274021,
            0.07097411]])
    >>> dataset.cell_types
    array([1, 7, 1, 9, 4, 7, 7, 6, 9, 1, 6, 1, 6, 4, 3, 4, 1, 1, 1, 1, 9, 9,
           4, 6, 0, 4, 1, 7, 1, 4, 4, 4, 3, 1, 1, 3, 7, 3, 3, 1, 1, 5, 4, 3,
           1, 1, 4, 6, 1, 1, 1, 1, 1, 9, 6, 6, 1, 3, 1, 1, 4, 3, 1, 4, 1, 4,
           7, 1, 7, 1, 6, 6, 3, 9, 6, 1, 6, 3, 6, 9, 4, 1, 6, 6, 1, 9, 6, 6,
           4, 1, 4, 6, 4, 4, 4, 6, 3, 3, 7, 1])
    >>> dataset.sample_index
    array(['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0'], dtype='<U1')

As you can see, ``PyCytoData`` manages uses simple array instead of dictionaries
because it has capabilities to manage batches samples. Of course, the details
of this package is out of the scope of this project, but to find out more
about ``PyCytoData``, you can read the detailed documentation written by its
lovely devs `here <https://pycytodata.readthedocs.io/en/latest/>`_.

As always, you can use the results for downstream analyses as usual.


-------------------------------------

*****************************************************
Emulation Mode: Real Data-Based Simulation
*****************************************************

If you already have the expression matrices and datasets of your dream but you
still want to experience the glory of Cytomulate, let us introduce Emulation
Mode. In this mode, Cytomulate uses an existing dataset as a basis for
generating new expressions. The key advantage of this mode is that it can
quickly replicate existing data without the need of resampling.

To use this mode, we require prior information on cell types, which will
ensure the best approximation. To do this, let's use ``PyCytoData`` again!
First, let's load our existing datasets:

.. code:: python

    >>> from PyCytoData import DataLoader
    >>> data = DataLoader.load_dataset(dataset="Levine13")
    >>> data.expression_matrix
    array([[ 1.24334908e+02,  6.28371582e+01, -6.17444396e-01, ...,
             1.06896072e+02,  6.39934635e+00,  7.14621687e+00],
           [ 1.22633148e+02,  5.52684593e+01, -3.17519844e-01, ...,
             1.27218781e+02, -3.17452759e-01,  1.12626851e+00],
           [ 3.30561943e+01,  1.73848724e+01, -7.71313131e-01, ...,
             3.32087189e+02, -2.46072114e-01,  8.84189606e-01],
           ...,
           [ 3.49014664e+01,  4.32544184e+00,  8.33491230e+00, ...,
             2.79086884e+02,  1.60285759e+01,  3.90819855e+01],
           [ 1.70956116e+01,  9.30270076e-01, -1.08385071e-01, ...,
             3.84983948e+02,  4.54559469e+00,  9.67729034e+01],
           [ 1.04753265e+01, -7.23805502e-02, -5.91436803e-01, ...,
             5.08439331e+02,  2.38833976e+00,  1.06308832e+01]])
    >>> data.cell_types
    array(['Plasmacytoid DC', 'Plasmacytoid DC', 'Plasmacytoid DC', ...,
       'MEP', 'MEP', 'MEP'], dtype='<U17')

For those of you who are familiar the ``Levine13`` dataset, this will be right
at home! For others, this is a well-known benchmark dataset.

Now, to start cytomulating, the overall interface is very similar but with a
different class:

.. code:: python

    >>> from cytomulate import EmulationCytofData

    >>> cytof_data = EmulationCytofData()
    >>> cytof_data.initialize_cell_types(expression_matrix=data.expression_matrix,
                                         labels=data.cell_types)
    >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)

Now, let's look at our outputs:

.. code:: python

    >>> expression_matrices
    {0: array([[9.89849695e+01, 9.49059169e+00, 8.13185395e-01, ...,
                3.45029702e+01, 3.18472044e-01, 5.29247895e+02],
               [2.15393928e+02, 2.98800278e+01, 1.59270843e+00, ...,
                1.67047977e+02, 5.34828676e+00, 1.03305364e+02],
               [3.82536701e+02, 2.91190531e+02, 1.05922645e+02, ...,
                0.00000000e+00, 0.00000000e+00, 4.10342314e+01],
               ...,
               [2.21592856e+02, 1.08856275e+00, 6.48076690e-01, ...,
                0.00000000e+00, 0.00000000e+00, 3.48087446e+02],
               [1.72226786e+01, 0.00000000e+00, 7.60774586e+00, ...,
                0.00000000e+00, 4.61570855e+00, 1.45485442e+02],
               [0.00000000e+00, 1.21741897e+01, 2.83614833e+00, ...,
                5.37901263e+00, 6.80561586e+00, 3.41060042e+01]])}
    >>> labels
    {0: array(['Mature CD4+ T', 'NotGated', 'NotGated', 'Mature CD38lo B',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Naive CD4+ T',
               'NotGated', 'NotGated', 'CD11bhi Monocyte', 'CD11bmid Monocyte',
               'NotGated', 'Mature CD4+ T', 'Naive CD8+ T', 'CD11b- Monocyte',
               'NotGated', 'NotGated', 'Mature CD4+ T', 'NotGated',
               'Mature CD4+ T', 'Megakaryocyte', 'NotGated', 'NotGated',
               'NotGated', 'NotGated', 'Mature CD4+ T', 'NotGated', 'NotGated',
               'NotGated', 'NotGated', 'Megakaryocyte', 'NotGated',
               'Mature CD8+ T', 'NotGated', 'Mature CD8+ T', 'Mature CD4+ T',
               'NotGated', 'NotGated', 'Naive CD4+ T', 'NotGated',
               'CD11bhi Monocyte', 'NotGated', 'NotGated', 'NotGated', 'NotGated',
               'Megakaryocyte', 'NotGated', 'NK', 'NotGated', 'CD11bhi Monocyte',
               'Naive CD8+ T', 'Naive CD8+ T', 'NotGated', 'NotGated',
               'Mature CD4+ T', 'Naive CD8+ T', 'NotGated', 'NotGated',
               'Mature CD8+ T', 'NotGated', 'Mature CD38lo B', 'NotGated', 'NK',
               'NotGated', 'Mature CD8+ T', 'NotGated', 'NotGated',
               'Mature CD8+ T', 'CD11bhi Monocyte', 'NotGated', 'NotGated',
               'Mature CD8+ T', 'NotGated', 'HSC', 'Erythroblast', 'NotGated',
               'Mature CD8+ T', 'NotGated', 'NotGated', 'NotGated', 'NotGated',
               'Erythroblast', 'Mature CD8+ T', 'Mature CD4+ T', 'Megakaryocyte',
               'Mature CD8+ T', 'NotGated', 'NotGated', 'NotGated',
               'Megakaryocyte', 'NotGated', 'NotGated', 'Naive CD4+ T',
               'NotGated', 'NotGated', 'Mature CD4+ T', 'NotGated',
               'Erythroblast'], dtype='<U17')}


PyCytoData Output
------------------------

If you have fallen in love with ``PyCytoData``, good news: the emulation mode is compatible with
``PyCytoData`` output as well! The procedure is exactly the same as the Creation Mode:

.. code:: python

    >>> from cytomulate import EmulationCytofData

    >>> cytof_data = EmulationCytofData()
    >>> cytof_data.initialize_cell_types(expression_matrix=data.expression_matrix,
                                         labels=data.cell_types)
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)


It's as simple as this! The rest is the same as the Creation Mode!


**Congratulations!!** You've officially made it through the Quickstart Guide! You're
on track to become a Cytomulate expert! Now, you can read more about settings and complex
simulation situations `in this tutorial <https://cytomulate.readthedocs.io/en/dev/tutorial/complex.html>`_.
