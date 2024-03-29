![Logo](/assets/cytomulate.jpg)

# cytomulate
> A simulation package for Cytometry by Time-of-Flight (CyTOF)

[![forthebadge](https://forthebadge.com/images/badges/open-source.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)

| Branch | Release | CI/CD | Documentation | Code Coverage |
| --- | --- | --- | --- | --- |
| main | ![Badge1](https://img.shields.io/badge/Version-v0.2.0-success) | ![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=main) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/main/?badge=main) |  [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/dev/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |
| dev | ![Badge1](https://img.shields.io/badge/Version-v0.2.0-success) |![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=dev) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/dev/?badge=dev) | [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/dev/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |


## Installation

You can easily install ``cytomulate`` from either ``PyPI`` or ``conda``. For the former, use the following command:

```shell
pip install cytomulate

```

Or if you are using a conda environment, you can use the following command:

```shell
conda install -c normalizingflow cytomulate

```
If you wish to use ``PyCytoData``, you can install separately with more instructions [here](https://cytomulate.readthedocs.io/en/dev/installation.html).

### Dependencies

Good news: we didn't have to write `cytomulate` from scratch in assembly language! This means that we will need dependencies to install it. Below is a list of packages that you will need:

- numpy
- scipy
- scikit-learn
- networkx
- matplotlib
- tqdm

Most of these are pretty standard! And even better news: the installation instructions should automatically handle all the dependency issues. If you have a problem with installation, let us know and we're happy to help!

While the above are all core dependencies, we highly highly highly highly highly recommend `PyCytoData` as well! You can get all the benefits of an integrated pipeline! But of course, for those of you who don't fancy more dependencies, we understand as well!


## Examples

We have two modes: **Creation Mode** and **Emulation Mode**. The former is probabilistic-model based simulation without the need of datasets; the latter is based on existing datasets to match as much of the existing features as possible. Here, we give two quick examples of how they work.


### Creation Mode

To create your datasets, you can run the following:

```python
>>> from cytomulate import CreationCytofData
>>> cytof_data = CreationCytofData()
>>> cytof_data.initialize_cell_types()
>>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)
```
The ``expression_matrices`` is a dictionary that contains the expression matrix from each sample. Correspondingly, ``labels`` is a dictionary that contains their cell types.


### Emulation Mode

This is a little bit more involved because we need existing data! If you already have your data, congratulations, you are good to go! For this demonstration, we use ``PyCytoData`` to load some example datasets instead (Of course, you will need to install [PyCytoData](https://pycytodata.readthedocs.io/en/latest/index.html) first if you wish to use it):

```python
>>> from cytomulate import EmulationCytoData
>>> from PyCytoData import DataLoader

>>> exprs = DataLoader.load_dataset(dataset="levine13")
>>> exprs.preprocess(arcsinh=True)
>>> cytof_data = EmulationCytofData()
>>> cytof_data.initialize_cell_types(expression_matrix=exprs.expression_matrix,
...                                  labels=exprs.cell_types)
>>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)
```
This is it!

### Working with PyCytoData

![PyCytoData](/assets/pycytodata.jpg)

We're fully compatible with ``PyCytoData``! As you've seen above, you can use ``PyCytoData`` to download datasets! If you're familiar with that interface and in love with its easy workflow, you can have ``cytomulate`` output a ``PyCytoData`` object as well:

```python
>>> from cytomulate import CreationCytofData
>>> cytof_data = CreationCytofData()
>>> cytof_data.initialize_cell_types()
>>> simulation_data = cytof_data.sample_to_pycytodata(n_samples = 1000)
```
This will allow you to use all the downstream capabilities of ``PyCytoData``.

### Command-Line Interface (CLI)

If you prefer to use cytomulate from the command-line, you've got that option as well! One **caveat** is that this mode requires ``PyCytoData`` to be installed. To run the Creation Mode, you can do:

```shell
python -m cytomulate \
	--creation \
	--n_cells 1000 \
	-o <your_dir_here>
```

To run the emulation mode, you can run the following:

```shell
python -m cytomulate \
	--emulation \
	--n_cells 1000 \
	-o <your_dir_here> \
	--exprs <you_path_to_exprssion_matrix> \
	--cell_types <you_path_to_cell_types>
```
Of course, we have much more customization options! For more details, read our [tutorial here](https://cytomulate.readthedocs.io/en/dev/tutorial/cli.html).

## Documentation

For more detailed documentation on ``cytomulate``, please visit our [website](https://cytomulate.readthedocs.io/)! You will find detailed tutorials,
guidelines, development guides, etc.

Our documentation is built automatically on the cloud! If you wish to build locally, check our detailed guide [here](https://cytomulate.readthedocs.io/en/latest/change/build.html)!

## Latest Release: v0.2.0

Welcome to Cytomulate v0.2.0! Hooray! We are not only bringing documentation enhancements, but we
are also introducing a new feature for more accurate simulations!

### Changes and New Features
- The `utilities.univariate_noise_model()` method:
    - Added `half_normal` option to the `noise_distribution` parameter
    - Changed the default `noise_distribution` to `uniform` (This is a **breaking change** because of the benefits to simulated results).
    - A warning is given when no user-specified `noise_distribution` is supplied to warn the breaking change
- Added the `utilities.estimate_noise_model()` method to estimate the noise present in the data
- Added a built-in estimation procedure to match the amount of zeroes observed in the dataset

### Improvements
- Added 4 more detailed tutorials on `our documentation website <https://cytomulate.readthedocs.io>`_
- Improved docstrings with more details on key parameters
- Updated the lastest references and links


## References

If you are cytomulating in your workflow, citing [our paper](https://doi.org/10.1186/s13059-023-03099-1) is appreciated:

```
@article{Yang2023,
	title = {Cytomulate: accurate and efficient simulation of CyTOF data},
    volume = {24},
    ISSN = {1474-760X},
    url = {http://dx.doi.org/10.1186/s13059-023-03099-1},
    DOI = {10.1186/s13059-023-03099-1},
  	number = {1},
	journal = {Genome Biology},
	publisher = {Springer Science and Business Media LLC},
	author = {Yang,  Yuqiu and Wang,  Kaiwen and Lu,  Zeyu and Wang,  Tao and Wang,  Xinlei},
	year = {2023},
	month = nov 
}
```