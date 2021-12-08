![Logo](/assets/cytomulate.jpg)

# cytomulate
> A simulation package for Cytometry by Time-of-Flight (CyTOF)

| Branch | Release | CI/CD | Documentation | Code Coverage |
| --- | --- | --- | --- | --- |
| main | ![Badge1](https://img.shields.io/badge/Version-PreRelease-success) | ![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=main) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/main/?badge=main) | [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/main/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |
| dev | ![Badge1](https://img.shields.io/badge/Version-PreRelease-success) |![Tests](https://github.com/kevin931/cytomulate/actions/workflows/ci.yml/badge.svg?branch=dev) | [![Documentation Status](https://readthedocs.org/projects/cytomulate/badge/?version=dev)](https://cytomulate.readthedocs.io/en/dev/?badge=dev) | [![codecov](https://codecov.io/gh/kevin931/cytomulate/branch/dev/graph/badge.svg?token=F5H0QTXGMR)](https://codecov.io/gh/kevin931/cytomulate) |

## Notice

This package is currently under development! **All rights reserved** until an official project license is added by our team. For questions and permissions, please contact our team!

## Project Dashboard

Check out our new [project dashboard](https://github.com/kevin931/cytomulate/projects/1). 

## Unit Testing

We will be using ``pytest`` for unit testing (talk to Kevin if this is an issue or if there are justifications for other solutions). Let's strive for 100% test coverage eventually! To get started, you will need to install the following packages:

```
pytest
pytest-mock
pytest-cov
coverage
```

To run tests by installing the package which is recommended, we can do the following:

```shell
python setup.py develop 
pytest --cov cytomulate
```

To run our tests on local files:

```shell
python -m pytest --cov cytomulate
```

If tests don't pass at this stage, it's okay. Kevin will be catching up on writing tests and talking with the team on what works and what not.

## CLI

We do have a CLI now! To test how it works, run the following: 

```shell
python -m cytomulate --version
python -m cytomulate -h
```

## Build Documentation

The documentation is built automatically on the cloud. Upon update, please check [our website](https://cytomulate.readthedocs.io/)!