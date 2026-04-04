# Configuration file for the Sphinx documentation builder.
#
# Full list: https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
import datetime

sys.path.insert(0, os.path.abspath('../..'))

# Mock imports that are unavailable during doc builds (e.g., pywin32 on Linux/CI)
autodoc_mock_imports = [
    'win32com', 'pythoncom', 'pywintypes',
    'win32com.client',
]

# -- Project information -----------------------------------------------------

project = 'pyHMT2D'
year = datetime.date.today().year
copyright = f'2020-{year}, Xiaofeng Liu'
author = 'Xiaofeng Liu'
release = '2.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.extlinks',
    'sphinx.ext.coverage',
    'sphinx_copybutton',
    'myst_parser',
]

# MyST-Parser settings
myst_enable_extensions = [
    'colon_fence',    # ::: directive syntax
]
myst_heading_anchors = 3
suppress_warnings = ['myst.header']

# Source file settings
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
# master_doc can point to .md now — Sphinx resolves by document name (no extension)
master_doc = 'contents'

# Autodoc settings
autoclass_content = "both"
autodoc_inherit_docstrings = True
add_module_names = False

# Pygments
pygments_style = 'friendly'

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_show_sourcelink = False

html_context = {
    'display_github': False,
    'github_user': 'psu-efd',
    'github_repo': 'pyHMT2D',
    'github_version': 'main/docs/',
    'menu_links_name': 'Getting Connected',
    'menu_links': [
        ('<i class="fa fa-github fa-fw"></i> Source Code', 'https://github.com/psu-efd/pyHMT2D'),
        ('<i class="fa fa-gavel fa-fw"></i> Contributing', 'https://github.com/psu-efd/pyHMT2D'),
    ],
}

html_theme_options = {
    'logo_only': True,
}

html_static_path = ['_static']
html_additional_pages = {'index': 'index.html'}

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': '',
    'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'pyHMT2D_API.tex', 'pyHMT2D API',
     'Xiaofeng Liu', 'manual'),
]

# -- Options for Epub output -------------------------------------------------

epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright
epub_exclude_files = ['search.html']

# -- Extension configuration -------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
}

def setup(app):
    app.add_css_file('css/custom.css')
