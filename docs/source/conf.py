# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from pathlib import Path

project = "GaussmarkovSRS"
copyright = "2021-2024 Marco Pizzocaro - Istituto Nazionale di Ricerca Metrologica (INRIM)"
author = "Marco Pizzocaro"
release = "0.1.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "autodoc2",
    "sphinx.ext.napoleon",
]

autodoc2_packages = ["../../src/gauss_markov_srs"]
autodoc2_render_plugin = "myst"


myst_enable_extensions = ["dollarmath", "amsmath"]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
