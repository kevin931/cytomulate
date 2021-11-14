# cytomulate
CyTOF Simulation: When in doubt, cytomulate it!

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