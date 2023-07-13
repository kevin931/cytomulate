Changelog
============

Here are the most recent releases and changes of cytomulate. Currently, we're still under developmet.
Therefore, we don't have any official releases. However, check out our git history to see what we're
doing!

Latest Release
---------------

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

.. toctree::
    :maxdepth: 1

    releases
    recent