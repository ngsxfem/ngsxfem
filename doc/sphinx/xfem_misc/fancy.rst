.. ngsxfem documentation master file, created by
   copy+paste and changes from ngsxtrefftz.

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

Fancy stuff
========================================

We give a short introduction to ngsxfem.

.. math::

    \newcommand{\Th}{{\mathcal{T}_h}} 
    \newcommand{\Fh}{\mathcal{F}_h} 
    \newcommand{\dom}{\Omega} 
    \newcommand{\jump}[1]{[\![ #1 ]\!]}
    \newcommand{\tjump}[1]{[\![{#1} ]\!]_\tau}
    \newcommand{\avg}[1]{\{\!\!\{#1\}\!\!\}}
    \newcommand{\nx}{n_\mathbf{x}} 
    \begin{align*} \begin{split}
        \begin{cases}
        -\Delta u = 0 &\text{ in } \Omega, \\
        u=g &\text{ on } \partial \Omega,
        \end{cases}
    \end{split} \end{align*}



Blablablubb
----------------------------------------

Blablablubb

.. jupyter-execute::
    :hide-output:
    :hide-code:

    from ngsolve import *
    from xfem import *


.. jupyter-execute::
    :hide-output:

    from ngsolve import *
    from xfem import *
    print("asdf")

Results in

.. jupyter-execute::
    :hide-code:

    from ngsolve import *
    from xfem import *
    print("asdf")


.. jupyter-execute::

    from ngsolve.webgui import Draw
    mesh = Mesh(unit_square.GenerateMesh(maxh=0.3))
    Draw(mesh)


Documentation
========================================

ABC
----------------------------------------
abc... -> LINK 

Notebooks
========================================

You can run them online ... 

.. toctree::
   :maxdepth: 1

   ../../external_dependencies/jupyter/cutfem.ipynb
