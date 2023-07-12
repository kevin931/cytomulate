####################
Emulation Mode
####################

Sometimes when we already have data on hand and we can't help but wanting more datasets, then this is the
time to Cytomulate with the **Emulation Mode**. As introduced in the `Quickstart Guide <https://cytomulate.readthedocs.io/en/dev/quickstart.html>`_,
the Emulation Mode allows users to input their own datasets and then simulate data based on the given
dataset. This is useful in a number of scenarios:

1. More datasets are needed to validate a method.
2. The current dataset has only 1 batch, but more than 1 is desired.
3. Biological interpretation with real protein channels and cell types are desired.
4. Analyses call for a basis on real data.

Of course, there are many more reasons to pick Emulation Mode over Creation Mode. The only prerequisite
is that users have to have their own dataset already. If not, they will have to find a dataset online
or actually perform some experiments in the lab. But for the sake of this tutorial, we will assume
that everyone has a dataset available to be emulated.

*******************
Setting Things Up
*******************

Before we start, let's make some necessary imports:

.. code-block:: python

    >>> from cytomulate import EmulationCytofData
    >>> from PyCytoData import DataLoader, FileIO, PyCytoData

Throughout the rest of the tutorial, we will use these functions directly.

.. note::
    ``PyCytoData`` is optional for the demonstration here, and it is not a required dependency
    for ``Cytomulate`` unless the CLI is used. We include examples with ``PyCytoData`` for
    its convenience in Python.

------------------------------

***************
Load Datasets
***************

First things first, we need to load our existing dataset into Python. For this tutorial, we would
like to assume that the data are stored in plain text (*e.g.* txt, csv, etc.) rather than the ``fcs``
format. We understand that the latter is popular in the field, but we do not yet have a good tutorial
on this. We would defer to the documentation of tools such as `fcspasrser <https://github.com/eyurtsev/fcsparser>`_
or `flowCore <https://bioconductor.org/packages/release/bioc/html/flowCore.html>`_. If you do have
``fcs`` files, you can use these tools to load the files, extract the expression matrices and cell
types, and then save the results as plain text files.


In this section, we will showcase a few methods of loading the datasets. For all methods, we will
save the expression matrix as ``dataset`` whereas the cell types as ``cell_types``, both of which
are NumPy arrays. Further, we generally assume that you don't have instrument channels, such as
"Time" and "Bead", in your dataset! If you do, you may need some preprocssing before passing them to
Cytomulate.

The `PyCytoData` Method
------------------------

We highly recommend the usage of PyCytoData to load the data into Python because it is quite easy
to do! If you have the data on hand, you can quickly load the expression matrix and cell types
with the following snippet of code:

.. code-block:: python

    >>> dataset = FileIO.load_delim("<Path_to_Expressions>", skiprows = 1, delim = "\t", dtype = float)
    >>> cell_types = FileIO.load_delim("<Path_to_Cell_Types>", skiprows = 0, delim = "\t", dtype = str)

You will replace the the paths with the actual paths to your files. Note that in the above example, we assume
that the expression matrix has the first row as protein marker names, which will be skipped. We further assume
that the file is saved as a tab-delimited file. You can customize the process to your own files.


If you would like to use the ``PyCytoData`` object, you can use the following example:

.. code-block:: python

    >>> data = FileIO.load_expression("<Path_to_Expressions>", col_names = True, delim = "\t", dtype = float)
    >>> cell_types = FileIO.load_delim("<Path_to_Cell_Types>", skiprows = 0, delim = "\t", dtype = str)
    >>> data.cell_types = cell_types

This has the advantage of storing the channel names in the object via the ``data.channels`` attribute. Now,
the only difference is that you need to use ``data.expression_matrix`` to access your expression
matrix.

PyCytoData has many other features and details which are all helpful in some way, but we cannot rewrite
the documentation here. You can read go ahead and read the `documentation <https://pycytodata.readthedocs.io>`_
if interested. 


The `numpy` Method
---------------------

The PyCytoData's IO methods are actually built on top of the NumPy package. If you wish to use the NumPy directly,
it is also easy to do so:

.. code-block:: python

    >>> import numpy as np
    >>> dataset = np.loadtxt("<Path_to_Expressions>", skiprows = 1, delim = "\t", dtype = float)
    >>> cell_types = np.loadtxt("<Path_to_Cell_Types>", skiprows = 0, delim = "\t", dtype = str)

which has a similar interface to the PyCytoData interface. This is arguably a simpler interface. However, NumPy does
not have the benefits of a PyCytoData object, which you may or may not need.

The `pandas` Method
--------------------

Of course, if you would like to read data using Pandas, you can do so with ease as well. We are going to read the
data and extract the necessary components to be used by Cytomulate. Here, we assume a few things for demonstration
purposes:
 
1. Everything is in one file, including a column with cell types.
2. The cell type column is named "Cell_Types".
3. Column names are available.

.. code-block:: python

    >>> import pandas as pd
    >>> df = pd.read_csv("<Path_to_Expressions>", sep = "\t", header = 0)
    >>> cell_types = df["Cell_Types"].to_numpy()
    >>> dataset = df.drop("Cell_Types", axis = 1).to_numpy()

As with all other methods, there are tons of variations based on the data available! You likely need some
daat wrangling, but the principles are simple!


----------------------------------

*******************
Simple Simulation
*******************

Once you have your datasets, you can run your simplest simulation:

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
                                         labels = cell_types,
                                         max_components = 9,
                                         min_components = 1,
                                         covariance_types = ("full", "tied", "diag", "spherical"))
    >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)


As opposed to the `Quickstart Guide <https://cytomulate.readthedocs.io/en/dev/quickstart.html>`_, here
we are more explicit about the parameters that users can supply or change. Note that all above
listed parameters are defaults, except for the ``expression_matrix`` and ``labels`` that need to be
user supplied. In the above function calls, we are doing the following:

1. Initialize a ``EmulationCytofData`` object with 1 batch.
2. Ask the object to emulate the given ``dataset`` with the given ``cell_types``. The GMM model is allowed to use between 1 and 9 components with one of the four covariance types.
3. We sample 100 cells in total for the dataset.
 
We would like to note that ``n_samples`` indicates the sample size, not the number of batches. In this case, we have
the following outputs:

.. code-block:: python

    >>> expression_matrices
    {0: array([[158.04318784, 143.54615373,  14.55784182, ...,  28.23468037,
                2.053792  ,   0.        ],
               [1.00986273,   6.6360988 ,   0.        , ...,   6.54024249,
                1.7951401 ,  19.13473223],
               [4.30063734,   5.7733213 ,   0.        , ...,  23.30509231,
                8.80240499,  21.81913369],
                ...,
               [191.03235818,  18.41559666,   0.        , ..., 333.3101886 ,
                0.        , 215.16717053],
               [310.09391901, 123.0910225 ,   0.        , ...,   0.        ,
                0.        , 164.39680483],
               [382.54861993,   8.39272611,   3.83182708, ...,  48.24234018,
                3.20327527, 184.68550649]])}
    >>> labels
    {0: array(['Mature CD38lo B', 'Erythroblast', 'Erythroblast',
               'CD11bhi Monocyte', 'Mature CD38lo B', 'NotGated', 'NotGated',
               'CD11b- Monocyte', 'NotGated', 'NotGated', 'NotGated', 'NotGated',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Mature CD4+ T',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Mature CD8+ T',
               'NotGated', 'Mature CD4+ T', 'NotGated', 'Mature CD4+ T',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Mature CD4+ T',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Erythroblast',
               'NotGated', 'NotGated', 'CD11bhi Monocyte', 'NotGated', 'NotGated',
               'CD11bhi Monocyte', 'CD11bmid Monocyte', 'NotGated',
               'Naive CD4+ T', 'Erythroblast', 'Plasmacytoid DC', 'Naive CD8+ T',
               'NotGated', 'Erythroblast', 'NotGated', 'CD11bhi Monocyte',
               'Megakaryocyte', 'NotGated', 'Mature CD4+ T', 'NotGated',
               'Mature CD4+ T', 'Mature CD4+ T', 'NotGated', 'NK', 'NotGated',
               'Naive CD8+ T', 'NotGated', 'NotGated', 'Mature CD8+ T',
               'NotGated', 'NK', 'NotGated', 'Mature CD8+ T', 'NotGated',
               'NotGated', 'NotGated', 'Mature CD4+ T', 'Mature CD38lo B',
               'CD11bhi Monocyte', 'NotGated', 'Mature CD38lo B', 'Naive CD8+ T',
               'Mature CD4+ T', 'NotGated', 'NotGated', 'Erythroblast',
               'NotGated', 'NotGated', 'NotGated', 'NotGated', 'Mature CD4+ T',
               'Mature CD4+ T', 'Erythroblast', 'Mature CD38lo B', 'Erythroblast',
               'NotGated', 'NotGated', 'Naive CD8+ T', 'NotGated',
               'Mature CD8+ T', 'NotGated', 'NotGated', 'Naive CD8+ T',
               'Mature CD8+ T'], dtype='<U17')}

Notice that each of the output is a dictionary of length 1 because we are only simulating 1 batch. We can
see that inside of the dictionary, each is an array. You can extract these arrays or use the
``expression_matrices`` and ``labels`` directly for downstream analyses.

Another thing to note is that the ``cytof_data.sample(n_samples = 100)`` call has 4 returns in a Tuple,
and the example above simply used Tuple unpacking. As the names suggest, the first is expression matrices,
and the second is the labels. In this case, we are not using the third and fourth returns, which are reserved
for pseudo-time and cell hierarchy in the cellular trajectory simulation from the
`Complex Simulation <https://cytomulate.readthedocs.io/en/dev/tutorial/complex.html>`_.

------------------------------

*************************
GMM and Model Selection
*************************

Cytomulate supports numerous model configurations for each dataset. Since the overall model is based on the
Gaussian Mixture Model (GMM), there are two major tuning parameters that are needed:

1. The number of components
2. The covariance matrix structure

As you can see from the example above, you can specify these manually using the simple simulation or go
with the defaults. If you have a specific configuration in mind (e.g. 3 components with a full covariance
matrix), you can simulate with the following:


.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 3,
    ...                                  min_components = 3,
    ...                                  covariance_types = "full")

which tells Cytomulate to use 3 components exactly with a full covariance matrix.

For components, you can choose any you like (except that they should be positive integers! We can't owe you
any components even if you wish!). For covariance matrix, you have four choices: "full", "diag", "tied", or
"spherical". Look at ``sklearn``'s
`documentation here <https://scikit-learn.org/stable/modules/generated/sklearn.mixture.GaussianMixture.html#sklearn.mixture.GaussianMixture>`_
for an explanation of how each of the four options work. In short, if you want everything, the "full" is a
good choice; or if you want independent channels, you need to choose "diag". While we don't recommend the
latter as we saw that correlation between protein channels are present, you can certainly try it if you
see a need.


Model Selection
----------------

As evident from the example above, there are many choices regarding these two crucial parameters. By default,
Cytomulate performs a model selection procedure using the
`Bayesian Information Criterion (BIC) <https://en.wikipedia.org/wiki/Bayesian_information_criterion>`_
for model selection. We perform this for both the number of components and the covariance structure.

For example, if you wish to slect between 1 and 5 components (inclusive) but with the full covariance
matrix, you can run the following code:

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 5,
    ...                                  min_components = 1,
    ...                                  covariance_types = "full")


Or if you wish to fix the number of components while selecting the type of covariance matrices,
run this instead:

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 5,
    ...                                  min_components = 5,
    ...                                  covariance_types = ("full", "tied", "diag", "spherical"))

Note that as long as ``max_component == min_components``, the specified number will be used. Also,
you can specify a subset of teh ``covariance_types`` without having to run all four.

Whether you wish to run model selection or not is up to your dataset and need. By default, Cytomulate
does perform model selection to allow for the best fit for the data. The downside is that this process
does slow down the algorithm because the model needs to be fit multiple times.

Defaults and Recommendations
------------------------------

As mentioned, model selection is performed by Cytomulate. If you would like to select the parameters
yourself, we have a few heuristics that may be helpful:

1. When in doubt, use the "full" covariance matrix to allow for the most flexible model.
2. the default betwene 1 and 9 components per cell type is often sufficient (hence the default).
3. If setting the components manually, perform an exploratory data analysis on your data first to see the data. Increase the numner of components if Cytomulate doesn't fit well.

Of course, there is no one-size-fit-all recommendation! Try it and find out!


----------------------------------


****************************
Cell Abundance
****************************

First, let us define cell abundance: The proportion of each cell type in a given sample. From a first
glance, this is incredibly vague because, and then, it may occur to people that this is somewhat trivial
if we are not dealing with simulated datasets. However, being a good simulator, Cytomulate allows users
to tweak the cell abundance in their dataset or use a real-data-based decision rule for simulation. But
before diving into details, we will discuss a bit of context for Cytomulate.

As defined in our paper, Cytomulate fits a GMM for each cell type, meaning that there are as many GMMs
as cell types. The cell types are user specified naturally. It may be tempting to fit one single GMM
to all the data, which is definitely possible, but the caveate is that a very complex model is needed.
So Cytomulate tries to take advantage of cell types or clusters as given by the user. Then, when the GMMs
are estimated, the next step is to sample from them, which is the ``sample()`` method seen earlier. Recall
that the method takes only 1 argument, which is the number of cells in the output rather than the number
of cells in each cell type. Thus, this section deals with how we estimate the makeup of our dataset.

Default: Data-Based
-----------------------

By default, we learn from the dataset itself. So if you use the following snippet as shown earlier,

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 9,
    ...                                  min_components = 1,
    ...                                  covariance_types = ("full", "tied", "diag", "spherical"))
    >>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 100)

Cytomulate will estimate the proportion of cells from real data. This is generally not affected by other
parameters, such as the number of components and the number of batches. We recommend the default if
you wish the simulated dataset to be as similar to the simulated data as possible.


Random Cell Abundance
-----------------------

There are some situations in which you may want to change the cell abundance, such as

- inducing more variance
- wanting more rare cell types
- needing some difference from the reference dataset.

In such cases, we can allow cell abundance to be randomly generated. There are in general two settings you
can choose. If you wish each cell type to have **equal probability**, then you should run the following
snippet:

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 9,
    ...                                  min_components = 1,
    ...                                  covariance_types = ("full", "tied", "diag", "spherical"))
    >>> cytof_data.generate_cell_abundance(use_observed = False,
                                           is_random = False)

.. note::
    It is uncessary to call ``generate_cell_abundance`` if you wish to use the observed cell abundance
    since it's the default behavior. 


On the other hand, if you wish to add some randomness into the mix, the probability for each cell type
can be generated with the Dirichlet distribution by simply setting ``is_random`` to ``True``:

.. code-block:: python

    >>> cytof_data = EmulationCytofData(n_batches = 1)
    >>> cytof_data.initialize_cell_types(expression_matrix = dataset,
    ...                                  labels = cell_types,
    ...                                  max_components = 9,
    ...                                  min_components = 1,
    ...                                  covariance_types = ("full", "tied", "diag", "spherical"))
    >>> cytof_data.generate_cell_abundance(use_observed = False,
                                           is_random = True)

This last example is especially useful is you wish to have a sample that is slightly different from
real samples.