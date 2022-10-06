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


.. mdinclude:: ../../README.md


Introduction to ngsxfem
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


Blabla 1
----------------------------------------
blablabla

.. math::

    \begin{align*} \begin{split}
    \mathbb{T}^p(K):=\big\{
    f\in\mathbb{P}^p(K) \mid \Delta f = 0
    \big\},
    \qquad p\in \mathbb{N}.
    \end{split} \end{align*}



Blablablubb
----------------------------------------

Blablablubb

..
   .. jupyter-execute::
       :hide-output:
       :hide-code:

       from ngsolve import *
       from ngsxfem import *



More
----------------------------------------

..
  .. _notebooks: notebooks/index.html


.. _documentation: 
Documentation
========================================

ABC
----------------------------------------
abc...

Notebooks
========================================

You can run them online ... 

.. toctree::
   :maxdepth: 1

   external_dependencies/jupyter/cutfem.ipynb
