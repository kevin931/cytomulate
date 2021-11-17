import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'cytomulate'
copyright = '2021, cytomulate developers'
author = 'cytomulate developers'
release = '0.0.0'

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode"
]


templates_path = ['_templates']
exclude_patterns = []
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
master_doc = "index"