# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'HypoMap3D'
copyright = '2026, Max Planck Institute for Metabolism Research'
author = 'MPI-SF | Neurogenomics'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'nbsphinx',
    "myst_parser",
    "sphinx_book_theme",
    'sphinx.ext.githubpages',
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
    "sphinxcontrib.mermaid",
]

# Use Markdown files as source
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['README.md']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_book_theme'
html_title = project
html_logo = '_static/hm3d-logo.png'
html_static_path = ['_static']
html_theme_options = {
    "show_toc_level": 0,
    # "announcement": "The HypoMap3D is still in development - please wait for official release.",
}

# Bibliography Configuration ------------------------------------------------- #
bibtex_bibfiles = ['references.bib']

# Notebook Configuration ('nbsphinx') ---------------------------------------- #
nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
