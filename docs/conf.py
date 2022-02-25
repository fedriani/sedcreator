# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import matplotlib
from datetime import datetime
import os
import sys
sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------

project = 'sedcreator'
author = 'Ruben Fedriani'
copyright = f'2022-{datetime.utcnow().year}, {author}'

# The full version, including alpha/beta/rc tags
release = '0.6.8'


# -- General configuration ---------------------------------------------------

# The name of an image file (within the static path) to use as favicon
# of the docs. This file should be a Windows icon file (.ico) being 16x16
# or 32x32 pixels large.
html_favicon = os.path.join('_static', 'favicon.ico')

# Configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3/',None),
    'numpy': ('https://numpy.org/doc/stable/',None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/',None),
    'matplotlib': ('https://matplotlib.org/stable/',None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'photutils': ('https://photutils.readthedocs.io/en/stable/', None),
    'h5py': ('https://docs.h5py.org/en/stable/', None)}

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon','matplotlib.sphinxext.plot_directive']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
