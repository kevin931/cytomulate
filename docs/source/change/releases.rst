###########
Releases
###########

This is a complete history of ``cytomulate`` releases.

-------------------

**************
v0.2.0
**************

Welcome to Cytomulate v0.2.0! Hooray! We are not only bringing documentation enhancements, but we
are also introducing a new feature for more accurate simulations!

Changes and New Features
--------------------------

- The `utilities.univariate_noise_model()` method:
    - Added `half_normal` option to the `noise_distribution` parameter
    - Changed the default `noise_distribution` to `uniform` (This is a **breaking change** because of the benefits to simulated results).
    - A warning is given when no user-specified `noise_distribution` is supplied to warn the breaking change
- Added the `utilities.estimate_noise_model()` method to estimate the noise present in the data
- Added a built-in estimation procedure to match the amount of zeroes observed in the dataset

Improvements
---------------
- Added 4 more detailed tutorials on `our documentation website <https://cytomulate.readthedocs.io>`_
- Improved docstrings with more details on key parameters
- Updated the lastest references and links

**************
v0.1.1
**************

This is our first maintenance update to be released to v0.1.x,
and we are packing in lots of enhancements! All changes are
regarding documentations!

Improvements
---------------
- Added 4 more detailed tutorials on `our documentation website <https://cytomulate.readthedocs.io>`_
- Improved docstrings with more details on key parameters
- Updated the lastest references and links

**************
v0.1.0
**************

Our **FIRST OFFICIAL RELEASE** is here! From now on, all our
releases will be supported with our standard support cycle.
Here you will find our release notes.

Changes and New Features
--------------------------

- Added Command-Line Interface with support for complex simulations
- Improved docstrings
- Improved documentations with tutorials

From Pre-release
------------------

These are listed for documetation reasons for the first official release.

- Support for ``Emulation Mode`` and ``Creation Mode``
- Support for complex simulations
- Availability on ``PyPI`` and ``conda``


**************
v0.0.2
**************

- This is the first stable pre-release of cytomulate
- A fix for critical installation error from **v0.0.1**.
- Availability on ``PyPI`` and ``conda``.

**********
v0.0.1
**********

- This is official pre-release of cytomulate.
- Introduction of Creation and Emulation mode.
    - Creation Mode for probabilistic model-based simulation
    - Emulation Mode for data-based simulation
    - Support for simulating batches, different cell types, cell differentiations, etc.
- Availability on ``PyPI`` and ``conda``.
