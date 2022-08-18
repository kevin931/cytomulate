####################
Complex Simulations
####################

In the `Quickstart Guide <https://cytomulate.readthedocs.io/en/dev/quickstart.html>`_, we've covered the
basics of simulation with each mode. In the basic configuration, we've changed the number of cells and
cell types. There are a few scenarios that are more complex than the basic scenario:

- Batch Effect
- Temporal Effect
- Cell Trajectories

This tutorial will walk you through how you can simulate each of the effects and what to expect. The interface
remains largely the same as the default interface with a few changes!

------------------------------

***************
Batch Effect
***************

As documented by literature and highlighted in `our study <https://doi.org/10.1101/2022.06.14.496200>`_,
batch effect exists when we obtain two similar samples. All else being equal, there can still
be differences between the two. Typically, the batch effect is a nuisance effect because researchers
usually want to analyze their samples together. The question you may ask is the following: why do
we want to generate such a nuisance effect? It turns out that we've thought about this as well!
The answer is rather simple and cute: For those studies and packages that perform batch normalization,
it is critical that they can have access to data with batch effects! In this case, Cytomulate will
be their dream package! 

We will leave the details to the paper if you are so inclined. To get started, let's generate
some batches:

.. code-block:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData(n_batches = 2)
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_overall_batch_effects(variance=0.1)
    >>> cytof_data.generate_local_batch_effects(variance=0.1)
    >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)

Generating two batches is as simple as this. For a more dramatic batch effect, use a larger variance, but
don't go overboard. Let's look at the data we have:

.. code-block:: python

    >>> expression_matrices
    {0: array([[0.99954052, 0.05126085, 0.04319021, ..., 0.05639346, 0.03011954,
                0.01908699],
               [0.99437961, 0.0360193 , 0.01607398, ..., 0.05581756, 0.04001113,
                0.02366643],
               [0.99559529, 0.03337728, 1.26513989, ..., 0.0286802 , 0.05832995,
                0.03264525],
               ...,
               [0.99919022, 0.02938261, 0.01392464, ..., 0.04645918, 0.05905717,
                0.02834046],
               [0.58031363, 0.03062228, 0.02477842, ..., 1.00096896, 0.04413759,
                0.01627813],
               [0.04929953, 0.04749511, 0.55356669, ..., 0.01676425, 0.05299929,
                0.01542087]]),
     1: array([[0.06279869, 0.57782953, 0.57845272, ..., 0.01318464, 1.29094045,
                0.04493943],
               [0.0378855 , 0.05237571, 0.56344615, ..., 0.01938484, 0.03924808,
                0.01164543],
               [0.04778554, 0.56998005, 0.58171915, ..., 0.03695538, 1.28998698,
                0.04139365],
               ...,
               [0.08989788, 0.04364689, 0.55921678, ..., 0.0162798 , 0.04618635,
                0.01405968],
               [0.05427921, 0.56828206, 0.58283406, ..., 0.02975604, 1.267897  ,
                0.04349445],
               [0.03539387, 0.04002488, 0.54994117, ..., 0.022102  , 0.03705235,
                0.01502765]])}

As you can see, we have two different samples in this case, which are stored in a dictionary. To
access the sample, you can simply use the dictionary key:

.. code-block:: python

    >>> expression_matrices[0]
    array([[0.99954052, 0.05126085, 0.04319021, ..., 0.05639346, 0.03011954,
            0.01908699],
           [0.99437961, 0.0360193 , 0.01607398, ..., 0.05581756, 0.04001113,
               0.02366643],
           [0.99559529, 0.03337728, 1.26513989, ..., 0.0286802 , 0.05832995,
               0.03264525],
           ...,
           [0.99919022, 0.02938261, 0.01392464, ..., 0.04645918, 0.05905717,
               0.02834046],
           [0.58031363, 0.03062228, 0.02477842, ..., 1.00096896, 0.04413759,
               0.01627813],
           [0.04929953, 0.04749511, 0.55356669, ..., 0.01676425, 0.05299929,
               0.01542087]])

The labels are stored in the same way. The procedures and outputs will remain the same for
Emulation Mode.

Global vs. Local Batch Effects
-----------------------------------

As you may have noticed from above, we have two additional steps: generating both overall
and local batch effects. They are so named because of the following reasons:

1. Batches can be globally different.
2. Batches may differ according to specific channels. In other, channels and global batch effect can have an interaction.

We refer to the former as overall batch effect, whereas the latter is named local batch effect.
Of course, you don't have to call both functions if you prefer only one of the effects. However,
we do recommend both for the most accurate results. The same logic applies for more than 2 batches.


``PyCytoData`` Object
------------------------

As you suspect, we can do this using ``PyCytoData``! Yay! To do this, simply have it output
to ``PyCytoData``:

.. code-block:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData(n_batches = 2)
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_overall_batch_effects(variance=0.1)
    >>> cytof_data.generate_local_batch_effects(variance=0.1)
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)

Now, let's look at the samples:

.. code-block:: python

    >>> dataset.n_samples
    2
    >>> dataset.sample_index
    array(['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',
           '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',
           '1', '1', '1', '1', '1'], dtype='<U1')

``PyCytoData`` stores the data differently by combining both samples in one object.
You can of course subset it accordingly.


*****************
Temporal Effect
*****************

Well, you may think: "Fine, we have batch effect. But are there really more effects?"
The answer is of course yes! There is also the temporal effect. For those who have preprocessed
CyTOF datasets, you know that one common step is to perform bead normalization, which is to
correct the temporal effect. So, we can generate temporal effects for those fine folks who
need it as well. To do so, it is very easy: 


.. code-block:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_temporal_effects(variance=0.1)
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)

As usual, you can change the variance to control how much of a temporal effect there will be.
By default, this uses the **Brownian Bridge** method. This overall code structure is also valid
for Emulation mode.

Polynomial Temporal Effect
-----------------------------

If you don't like Brownian Bridge, you can use a different method to geterate temporal effect.
In this case, you can specify the coefficients of a polynomial equation:

.. code-block:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_temporal_effects(variance=0.1, coefficients=[1,-1,0.5])
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)

This will fit the following polynomial:

.. math::
    
    1 - x + 0.5 \times x^2

This interface is the same as Numpy's polynomial.

.. note::

    The ``coefficients`` only specifies the shape of the polynomial, but not the exact
    polynomial that will be fitted.

You will still need to specify the variance.


Spline Temporal Effect
-------------------------

The third option is using spline. In this case, you will need to specify the interval
points:

.. code-block:: python

    >>> from cytomulate import CreationCytofData
    >>> import numpy as np

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_temporal_effects(x={0:np.linspace(0, 1, 10),
                                                1:np.linspace(0, 1, 10)},
                                             y={0:np.zeros(10),
                                                1:np.zeros(10)})

Here, you don't need to specify the variance.

----------------------------------

**************************
Cellular Trajectory
**************************

Flashback to high school: your biology teacher repeats how cells differentiate and replicate.
This is very much relavant in this case because different cell types can be related: as cells
mature, they can change. Thus, for related cell types, we can infer its trajectory by considering 
the relationship between them. 

Cytomulate, by default, takes this relationship into consideration. In the Emulation Mode, this is
true by definition because we are emulating what already exists in a given dataset. In the Creation
Mode, we use trees to mimic the relationships between cell types. However, the difference between
the default and trajectory is that cells do not "differentiate" in the default settings. In other
words, cells don't really "move" between the nodes on the tree. In order to have continuous movement
to mimic a differentiation path, this is what it's about! Now, before we bore you with text, let's
start cytomulating:

.. code-block:: python

    >>> from cytomulate import CreationCytofData

    >>> cytof_data = CreationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_cell_graph()
    >>> expression_matrices, labels, pseudo_time, children_cell_type = cytof_data.sample(n_samples = 100)

Notice that we did save two extra outputs. Let's look at the outputs:

.. code-block:: python

    >>> pseudo_time
    {0: array([[9.74536513e-01, 9.96081942e-01, 9.99992100e-01, ...,
                4.48842711e-06, 1.37760352e-02, 7.73697266e-01],
               [1.05592573e-03, 7.38913540e-01, 5.83574946e-03, ...,
                1.81810218e-01, 2.45327695e-04, 7.06481155e-02],
               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, ...,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
               ...,
               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, ...,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00],
               [9.55611504e-02, 2.67576663e-01, 9.92857547e-01, ...,
                9.44378391e-01, 4.19920559e-01, 9.01649202e-01],
               [0.00000000e+00, 0.00000000e+00, 0.00000000e+00, ...,
                0.00000000e+00, 0.00000000e+00, 0.00000000e+00]])}
    >>> children_cell_type
    {0: array(['6', '8', 'None', 'None', 'None', 'None', 'None', '8', 'None',
               'None', 'None', 'None', 'None', 'None', '6', 'None', 'None', '1',
               'None', 'None', 'None', 'None', '6', 'None', '6', '6', '6', 'None',
               'None', 'None', 'None', '3', '6', '8', 'None', 'None', '6', '4',
               'None', 'None', '6', 'None', 'None', '6', 'None', 'None', 'None',
               'None', '6', '6', 'None', '7', '7', 'None', 'None', '7', 'None',
               'None', '7', 'None', 'None', 'None', 'None', 'None', '6', '8',
               'None', 'None', '8', 'None', 'None', '6', '6', '1', 'None', 'None',
               'None', '1', 'None', '8', 'None', 'None', '7', 'None', '8', '6',
               '6', 'None', '1', 'None', 'None', 'None', 'None', 'None', 'None',
               'None', 'None', 'None', '8', 'None'], dtype='<U21')}

As the name imply, we have more information than just the expression matrix and cell types.
What we called the ``pseudo_time`` matrix describes the time step between two nodes. If it
is 0, it is firmly at the parent node; if it is close to 1, then it is more more its child
node than the parent node. In other words, we can have cell type A that has almost 
differentiated to cell type B. Here, the ``children_cell_type`` matrix come it. It shows
us what the child subtype is. If it's ``None``, then these cells no longer differentiate.

In the case of emulation mode, use the same workflow but add the following line before
sampling:

.. code-block:: python

    cytof_data.generate_cell_graph()

This is all you need!

**Where is PyCytoData?** This must be your burning question! Unfortunately, ``PyCytoData``
currently does not support storing these information. We are investigating the viability
of integrating these into ``PyCytoData``. Don't worry, Cytomulate remains a proud member
of the PyCytoData Alliance.

-------------------

*********************
Using Everything
*********************

Cytomulate is like a buffet because you can get everything all at once. In this case,
all you have to do is to add the effects sequentially after initializing the model and
cell types but before sampling.

Here is an example:

.. code:: python

    >>> from cytomulate import CreaationCytofData

    >>> cytof_data = CreaationCytofData()
    >>> cytof_data.initialize_cell_types()
    >>> cytof_data.generate_overall_batch_effects(variance=0.1)
    >>> cytof_data.generate_local_batch_effects(variance=0.1)
    >>> cytof_data.generate_temporal_effects(variance=0.1)
    >>> cytof_data.generate_cell_graph()
    >>> dataset = cytof_data.sample_to_pycytodata(n_samples = 100)

As usual, the procedure is the same for Emulation Mode.