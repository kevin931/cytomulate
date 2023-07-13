import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

project = 'cytomulate'
copyright = '2022-2023, cytomulate developers'
author = 'cytomulate developers'
release = '0.2.0'

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinxcontrib.autoprogram",
    "sphinx_autodoc_typehints",
    "sphinx_git"
]

templates_path = ['_templates']
exclude_patterns = []
html_theme = "sphinx_rtd_theme"
html_logo = "../../assets/cytomulate.jpg"
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': True,
    'logo_only': True,
    'display_version': False
}
master_doc = "index"