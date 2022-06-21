![Logo](/assets/cytomulate.jpg)

# cytomulate
> A simulation package for Cytometry by Time-of-Flight (CyTOF)

[![forthebadge](https://forthebadge.com/images/badges/open-source.svg)](https://forthebadge.com)
[![forthebadge](https://forthebadge.com/images/badges/made-with-python.svg)](https://forthebadge.com)

| Branch | Release | CI/CD | Documentation | Code Coverage |
| --- | --- | --- | --- | --- |
| main | ![Badge1](https://img.shields.io/badge/Version-PreRelease-success) | ![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=main) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/main/?badge=main) | [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/main/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |
| dev | ![Badge1](https://img.shields.io/badge/Version-PreRelease-success) |![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=dev) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/dev/?badge=dev) | [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/dev/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |


## Installation

You can easily install ``cytomulate`` from either ``PyPI`` or ``conda``. For the former, use the following command:

```shell

```

Or if you are using a conda environment, you can use the following command:

```shell

```
If you wish to use ``PyCytoData``, you can install separately with more instructions [here](https://cytomulate.readthedocs.io/en/dev/installation.html).

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
The ``expression_matricess`` is a dictionary that contains the expression matrix from each sample. Correspondingly, ``labels`` is a dictionary that contains their cell types.


### Emulation Mode

This is a little bit more involved because we need existing data! If you already have your data, congratulations, you are good to go! For this demonstration, we use ``PyCytoData`` to load some example datasets instead (Of course, you will need to install [PyCytoData](https://pycytodata.readthedocs.io/en/latest/index.html) first if you wish to use it):

```python
>>> from cytomulate import EmulationCytoData
>>> from PyCytoData import DataLoader

>>> exprs = DataLoader.load_dataset(dataset="levine13")
>>> cytof_data = EmulationCytofData()
>>> cytof_data.initialize_cell_types(expression_matrix=exprs.expression_matrix,
...                                  labels=exprs.cell_types)
>>> expression_matrices, labels, _, _ = cytof_data.sample(n_samples = 1000)
```

This is it!

### Working with PyCytoData

As you've seen above, you can use ``PyCytoData`` to download datasets! If you're familiar with that interface and in love with its easy workflow, you can have ``cytomulate`` output a ``PyCytoData`` object as well:

```python
>>> from cytomulate import CreationCytofData
>>> cytof_data = CreationCytofData()
>>> cytof_data.initialize_cell_types()
>>> simulation_data = cytof_data.sample_to_pycytodata(n_samples = 1000)
```
This will allow you to use all the downstream capabilities of ``PyCytoData``.

## Documentation

For more detailed documentation on ``cytomulate``, please visit our [website](https://cytomulate.readthedocs.io/)! You will find detailed tutorials,
guidelines, development guides, etc.

Our documentation is built automatically on the cloud! If you wish to build locally, check our detailed guide [here](https://cytomulate.readthedocs.io/en/latest/change/build.html)!

## Latest Release

This is our latest release:


## References

If you are cytomulating in your workflow, citing [our paper](https://doi.org/10.1101/2022.06.14.496200) is appreciated:

```
@article {Yang2022.06.14.496200,
	author = {Yang, Yuqiu and Wang, Kaiwen and Lu, Zeyu and Wang, Tao and Wang, Xinlei},
	title = {Cytomulate: Accurate and Efficient Simulation of CyTOF data},
	elocation-id = {2022.06.14.496200},
	year = {2022},
	doi = {10.1101/2022.06.14.496200},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2022/06/16/2022.06.14.496200},
	eprint = {https://www.biorxiv.org/content/early/2022/06/16/2022.06.14.496200.full.pdf},
	journal = {bioRxiv}
}
```