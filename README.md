# cytomulate
CyTOF Simulation: When in doubt, cytomulate it!

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
pytest . --cov
```

To run our tests on local files:

```shell
python -m pytest . --cov
```

If tests don't pass at this stage, it's okay. Kevin will be catching up on writing tests and talking with the team on what works and what not.

## CLI

We do have a CLI now! To test how it works, run the following: 

```shell
python -m cytomulate --version
python -m cytomulate -h
```

## Build Documentation

To build documentation locally, you will need the following packages:

```
sphinx
sphix-rtd-theme
sphinx-git
sphinxcontrib-autoprogram
sphinx-autodoc-typehints
```
which are both available as python packages in both ``PyPI`` and ``conda``. To build locally, run the following commands:

```shell
cd docs
make html
```

The automatic build process will create the ``/docs/build`` dictory, which contains all the local html files. You can open the documentation by opening ``/docs/build/html/index.html``. This whole directory is designed **not** to be tracked by Git as all such files are considered build artifacts.