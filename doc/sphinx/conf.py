from os.path import dirname, isdir, join
import os
import re
import subprocess
def get_version():
    """
    Gets the current version number.
    If in a git repository, it is the current git tag.
    Otherwise it is the one contained in the PKG-INFO file.
    """
    version_re = re.compile('^Version: (.+)$', re.M)
    d = dirname(__file__)+'/../..'

    if isdir(join(d, '.git')):
        # Get the version using "git describe".
        cmd = 'git describe --always --tags --match v[0-9]*'.split()
        try:
            version = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get version number from git tags')
            exit(1)

        # PEP 386 compatibility
        if '-' in version:
            version = '.post'.join(version.split('-')[:2])

        # Don't declare a version "dirty" merely because a time stamp has
        # changed. If it is dirty, append a ".dev1" suffix to indicate a
        # development revision after the release.
        with open(os.devnull, 'w') as fd_devnull:
            subprocess.call(['git', 'status'],
                            stdout=fd_devnull, stderr=fd_devnull)

        cmd = 'git diff-index --name-only HEAD'.split()
        try:
            dirty = subprocess.check_output(cmd).decode().strip()
        except subprocess.CalledProcessError:
            print('Unable to get git index status')
            exit(1)

        # if dirty != '':
            # version += '.dev1'

        # strip the v for pypi
        version = version[1:]

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)

    return version


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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


# -- Project information -----------------------------------------------------

project = 'ngsxfem'
copyright = '2022, Christoph Lehrenfeld and the ngsxfem developer'
author = 'Christoph Lehrenfeld, Fabian Heimann, Henry von Wahl, Janosch Preuß and the ngsxfem community'

# The full version, including alpha/beta/rc tags
release = get_version()

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
# extensions = [ "myst_nb" ]
extensions = ["sphinx.ext.autodoc","sphinx.ext.mathjax","sphinx.ext.todo","sphinx.ext.githubpages",
              "IPython.sphinxext.ipython_console_highlighting", "IPython.sphinxext.ipython_directive",
              # "jupyter_sphinx.execute",
              "jupyter_sphinx",
              "nbsphinx",
              "m2r2",
              "sphinxemoji.sphinxemoji",
              ]
# source_suffix = ['.rst', '.md']

# Run notebook configuration

# The template used when exporting from nbconvert
#   full  - Outputs the full HTML document [Default]
#   basic - Outputs a single div (with no additional resources)
run_notebook_export_template = 'basic'  # Default: 'full'

# Display the source links to the generated evaluated files
run_notebook_display_source_links = False  # Default: True

# Whether or not to evaluate the notebooks prior to embedding them
evaluate_notebooks = True  # Default: True

# START nbsphinx stuff
#increase timeout for cell execution, since some files take long to execute
nbsphinx_timeout = 100000

# If True, the build process is continued even if an exception occurs:
nbsphinx_allow_errors = False

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
.. raw:: html

    <style>
        .p-Widget {
            height: 400px;
        }
        .dg.main {
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li {
            list-style: none;
            margin-left: 0px;
        }
        div.p-Widget div div div div.dg ul li div.dg {
            margin-bottom: 0px;
        }
    </style>

.. only:: html
    .. role:: raw-html(raw)
        :format: html
"""

nbsphinx_widgets_path = ""
# END nbsphinx stuff

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'paper','env']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'alabaster'
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['bullets.css']

# html_theme_options = {
#     'github_user': 'ngsxfem',
#     'github_repo': 'ngsxfem',
#     # 'github_banner':True,
#     # 'travis_button':True,
#     'github_button':True,
#     'fixed_sidebar':True,
#     'font_size':10
#     # 'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
#     # 'analytics_anonymize_ip': False,
# }

master_doc = "contents"

# m2r_parse_relative_links = True

html_sidebars = {
   'index': [
        'about.html',
        'localtoc.html',
        # 'relations.html',
       ],
   '**': ['about.html',
          'localtoc.html',
       ],
}

nbsphinx_allow_errors = True

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# html_js_files = ['webgui_jupyter_widgets.js', 'webgui.js', 'tentswebgui.js']

def setup(app):
    app.add_css_file("custom.css")
