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
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'ngsxfem'
copyright = '2022, Christoph Lehrenfeld and the ngsxfem developer'
author = 'Christoph Lehrenfeld, Fabian Heimann, Henry von Wahl, Janosch Preuß and the ngsxfem community'

# The full version, including alpha/beta/rc tags
release = get_version()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_mdinclude",
    "jupyter_sphinx",
    "sphinx.ext.autodoc",
    "myst_nb",
    ]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

master_doc = "contents"

source_suffix = ['.rst', '.md']

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = ['bullets.css']

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

def setup(app):
    app.add_css_file("custom.css")


#-------nbsphinx (and probably also other) stuff
# The template used when exporting from nbconvert
#   full  - Outputs the full HTML document [Default]
#   basic - Outputs a single div (with no additional resources)
# run_notebook_export_template = 'basic'  # Default: 'full'

# # Display the source links to the generated evaluated files
# run_notebook_display_source_links = False  # Default: True

# # Whether or not to evaluate the notebooks prior to embedding them
# evaluate_notebooks = True  # Default: True

# # START nbsphinx stuff
# #increase timeout for cell execution, since some files take long to execute
# nbsphinx_timeout = 100000

# # If True, the build process is continued even if an exception occurs:
# nbsphinx_allow_errors = False

# # This is processed by Jinja2 and inserted before each notebook
# nbsphinx_prolog = r"""
# .. raw:: html

#     <style>
#         .p-Widget {
#             height: 400px;
#         }
#         .dg.main {
#             margin-left: 0px;
#         }
#         div.p-Widget div div div div.dg ul li {
#             list-style: none;
#             margin-left: 0px;
#         }
#         div.p-Widget div div div div.dg ul li div.dg {
#             margin-bottom: 0px;
#         }
#     </style>

# .. only:: html
#     .. role:: raw-html(raw)
#         :format: html
# """

# nbsphinx_widgets_path = ""
# # END nbsphinx stuff