"""The Command-Line Interface (CLI) of Cytomulate

For convenience, users can access cytomulate without the need of writing a script. Consistent
with python CLI syntax, this ``__main__.py`` can be accessed through ``python -m cytomulate``. 
The CLI has the advantages of convenience while also supporting most functionalities of the
packages.

:Example:

    Check version:

    .. code-block:: bash

        python -m cytomulate --version
        
    Get help:
    
    .. code-block:: bash

        python -m cytomulate -h
        
**Advantage**:
    * Shell-scripting friendly
    * Easy-to-use
    * No python knowledge required
    
**Drawbacks**:
    * Not as configurable as a python program
    * Top level functionality and integration only (i.e. No partial package usage)
"""


import cytomulate
import argparse


parser: "argparse.ArgumentParser" = argparse.ArgumentParser(description="cytomulate: CyTOF Simulation")
parser.add_argument("--version", action="version", version=cytomulate.__version__)


if __name__ == "__main__":
    cmdargs: argparse.Namespace = parser.parse_args()