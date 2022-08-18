#######################
Cytomulate CLI
#######################

We offer a concise and yet powerful Command-Line Interface (CLI) for those
who prefer shell scripting or running some simple simulations without
having to write Python code. Here you will find some examples and
tutorials on how to use ``cytomulate`` in this context and how you can
customize it to your needs.


------------------------

*********************
Arguments and Flags
*********************

Admittedly, the CLI mode is not as feature-rich as the full API, but there
is still plenty you can do. Here, you will find a list of all arguments
and flags at your disposal:

=========================== =============== ============== ============ ===================================================
Argument                      Mode             Inputs       Default         Functionality
--------------------------- --------------- -------------- ------------ ---------------------------------------------------
``--help`` or ``-h``          General          None          None        Print help and documentation
``--version``                 General          None          None        Print the version of cytomulate 
``--creation``                Creation         None          None        Run creation mode
``--emulation``               Emulation        None          None        Run emulation mode
``--n_cells``                 General          ``int``       None        The number of cells per batch.
``--trajectory``              Complex          None          None        Whether to generate trajectories for cell types.
``--batch_effect``            Complex          None          None        Whether to generate batch effect.
``--batch_effect_var``        Complex          ``float``     ``0.1``     The variance for batch effect.
``--temporal_effect``         Complex          None          None        Whether to generate temporal effect.
``--temporal_effect_var``     Complex          ``float``     ``0.1``     The variance for temporal effect.
``--n_batches``               Creation         ``int``       ``1``       The number of batches.
``--n_types``                 Creation         ``int``       ``10``      The number of cell types.
``--n_markers``               Creation         ``int``       ``20``      The number of markers.
``--n_trees``                 Creation         ``int``       ``2``       The number of trees.
``--exprs``                   Emulation        ``str``       None        The path to existing expression matrix.
``--exprs_colnames``          Emulation        None          ``False``   Whether the first row is colnames.
``--exprs_delim``             Emulation        ``str``       ``\t``      The delimiter of the expression matrix.
``--cell_types``              Emulation        ``str``       None        The path to existing cell types.
``--cell_types_colnames``     Emulation        None          ``False``   Whether the first row is colnames.
``--cell_types_delim``        Emulation        ``str``       ``\t``      The delimiter of the cell types.
``--out_dir`` or ``-o``       General          ``str``       None        Directory to save results.
``--make_new_dir``            General          None          ``False``   Whether to create the directory provided.
=========================== =============== ============== ============ ===================================================

If that table look a bit cryptic to you, don't panic! Here a few explanations.
The **Mode** indicates in which context you should use such arguments: the
**General** arguments are for both modes. Some arguments come with defaults,
making them completely optional, but for those that don't (as indicated by None),
you should provide the specifed type.

Now, we're going to explore how you can use these with examples and details.

------------------------------

*************************
Creation Mode
*************************

This is by far the easiest mode. The bare minimum is to specify ``--creation`` to start
along with ``n_cells`` for the number of cells per batch and ``-o`` the output directory:

.. code-block:: shell

    python -m cytomulate \
        --creation \
        --n_cells 1000 \
        -o <your_dir_here>

And the following three files will be saved to the directory you specified:

- exprs.txt
- cell_types.txt
- sample_index.txt

which are ``tsv`` files that are quite self-explanatory.

.. note:: 

    The first row of ``exprs.txt`` is the channel names. But for the cell types
    and sample index, the first row is not header.

You can also have a bit more customization by specifying the details of your samples:

.. code-block:: shell

    python -m cytomulate \
        --creation \
        --n_cells 1000 \
        --n_batches 2 \
        --n_types 5 \
        --n_markers 20 \
        --trajectory \
        -o <your_dir_here>

This will generate two batches with 1000 cells, 5 cell types, and 20 markers with trajectories.
Both batches will be combined and saved to a single file, but the ``sample_index.txt``
will delineate the indices accordingly.

If you wish ``cytomulate`` to create the directory you entered, you can add
the following flag:

.. code-block:: shell

    python -m cytomulate \
        --creation \
        --n_cells 1000 \
        -o <your_dir_here> \
        --make_new_dir


This is how easy creation mode is!

--------------------------

**********************
Emulation Mode
**********************

This mode is slightly more involved because you need to specify the file and
cell types so that the model can emulate it. However, this is not as hard as
it seems:

.. code-block:: shell

    python -m cytomulate \
        --emulation \
        --n_cells 1000 \
        -o <your_dir_here> \
        --exprs <you_path_to_exprssion_matrix> \
        --cell_types <you_path_to_cell_types>

You also have the option to generate trajectories here as well:

.. code-block:: shell

    python -m cytomulate \
        --emulation \
        --n_cells 1000 \
        --trajectory \
        -o <your_dir_here> \
        --exprs <you_path_to_exprssion_matrix> \
        --cell_types <you_path_to_cell_types>


If your reference exppression matrix and cell types are both tab separated
without a header, then this example is **all you need**. However, if your
files are saved differently, you can customize the IO process accordingly:

.. code-block:: shell

    python -m cytomulate \
        --emulation \
        --n_cells 1000 \
        -o <your_dir_here> \
        --exprs <you_path_to_exprssion_matrix> \
        --exprs_colnames \
        --exprs_delim , \
        --cell_types <you_path_to_cell_types> \
        --cell_types_colnames \
        --cell_types_delim ,

Here we've indicated that both files' first rows are the column names
and they are comma separated. These are both pretty standard. If your
files are saved in other formats or your cell types are saved with
your expression matrix, you will need to preprocess them separately
and save to these cytomulate-supported formats accordingly.

------------------------

*********************
Complex Simulations
*********************

As you may have noticed from the `Complex Simulation <https://cytomulate.readthedocs.io/en/dev/tutorial/cli.html>`_
section, we have a few complex simulation options. Above, you've seen the
``--trajectory`` flag. For details on each mode, you can read the linked
tutorial. However, doing this in the CLI is very easy:

.. code-block:: shell

    python -m cytomulate \
        --creation \
        --n_cells 1000 \
        --trajectory \
        --batch_effect \
        --temporal_effect \
        -o <your_dir_here>

This will generate trajectory, temporal effect, and batch effect all
with default settings (var=0.1 and Brownian Bridge for temporal effect).
One thing you can change is the variance:

.. code-block:: shell

    python -m cytomulate \
        --creation \
        --n_cells 1000 \
        --trajectory \
        --batch_effect \
        --batch_effect_var 0.2 \
        --temporal_effect \
        --temporal_effect_var 0.2 \
        -o <your_dir_here>

Here, we have changed the variance to 0.2.

.. note:: We do not support changing the temporal effect model. To do so, use the interactive mode instead.

The interface is the same for the Emulation Mode:

.. code-block:: shell

    python -m cytomulate \
        --emulation \
        --n_cells 1000 \
        --trajectory \
        --batch_effect \
        --batch_effect_var 0.2 \
        --temporal_effect \
        --temporal_effect_var 0.2 \
        -o <your_dir_here> \
        --exprs <you_path_to_exprssion_matrix> \
        --cell_types <you_path_to_cell_types>

For more customizations, read our other tutorials or use the interactive mode.