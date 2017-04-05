..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2014-2015 The ChemTools Development Team
    :
    : This file is part of ChemTools.
    :
    : ChemTools is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : ChemTools is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --


.. _usr_installation:

Installation
############

Downloading Code
================

The latest code can be obtained through Github (private at present),

  .. code-block:: bash

     $ git clone https://github.com/QuantumElephant/chemtools.git chemtools


Dependencies
============

The following dependencies will be necessary for ChemTools to build properly,

  .. code-block:: bash


Installation
============

To install ChemTools run:

  .. code-block:: bash

     $ ./setup.py install --user

If you prefer to run ChemTools from source folder, ``PYTHONPATH`` and ``CTDATA`` paths
need to be add into your **~/.bashrc** (Linux) or **~/.profile** (MacOS)

  .. code-block:: bash

     $ export PYTHONPATH=$PYTHONPATH:{path_to_chemtools_repo}/chemtools
     $ export CTDATA={path_to_chemtools_repo}/chemtools/data

Tests
=====

To ensue all the parts of ChemTools work properly, you are suggested to run the nosetests:

  .. code-block:: bash

     $ nosetests -v chemtools


.. _usr_lfs_installation:

Git LFS Installation
====================

`Git Large File Storage (LFS) <https://git-lfs.github.com/>`_ has been used to store the files
and images used in the examples' gallery. In order words, all files in ``doc/examples/*`` are
stored remotely, and only a text pointer to these files exists in the repository.
These files need to be downloaded separately if you would like to run the example scripts, go through
tutorials (using exactly the same files used) or make Chemtools HTML with sphinx.

First, you need to install Git LFS.

Mac OS
~~~~~~

  .. code-block:: bash

     $ sudo port install git-lfs
     $ brew install git-lfs

Linux OS
~~~~~~~~

In your download directory,

  .. code-block:: bash

     $ wget https://github.com/git-lfs/git-lfs/releases/download/v2.0.1/git-lfs-linux-amd64-2.0.1.tar.gz
     $ tar -zxvf git-lfs-linux-amd64-2.0.1.tar.gz
     $ cd git-lfs-2.0.1
     $ ./install.sh


.. _usr_lfs_files:

Downloading Files From Remote Repository
========================================

Having :ref:`Git LFS installed <usr_lfs_installation>`, download the examples files by running:

  .. code-block:: bash

     $ cd data/examples
     $ git lfs pull

To get a list of all the files tracked with Git LFS, use:

  .. code-block:: bash

     $ git lfs ls-files


Documentation Dependencies & Build
==================================

If you are interested in generating the documentation from source, the following
packages are also needed:

* Sphinx >=1.3.1: http://sphinx.pocoo.org/
* Sphinx Read-The-Docs theme >=0.1.8: https://github.com/snide/sphinx_rtd_theme
* sphinxcontrib-bibtex >= 0.3.5: https://pypi.python.org/pypi/sphinxcontrib-bibtex
* :ref:`Git Large File Storage (LFS) <usr_lfs_installation>`

To install the first three dependecies:

* **Ubuntu Linux 16.04**

  .. code-block:: bash

     $ sudo apt-get install python-sphinx python-sphinx-rtd-theme sphinxcontrib-bibtex

* **Ubuntu Linux 15.04 & 14.04 and Mac OS**

  .. code-block:: bash

     pip install --user --upgrade sphinx sphinx_rtd_theme sphinxcontrib-bibtex

The Sphinx Read-The-Docs theme customized for ChemTools can be obtained cloning the repository
as a submodule:

.. code-block:: bash

   $ cd doc/_themes
   $ git submodule update --init --recursive

To automatically generate API documentation and generate HTML (this requires ``data/examples`` files;
to obtain them please refer to :ref:`usr_lfs_files`) use the commands below:

.. code-block:: bash

   $ cd doc
   $ sphinx-apidoc -f -o source ../
   $ make clean
   $ make html

To open the documentation in your default browser, either click on ``doc/_build/html/index.html``
file directly, or run the command below from terminal:

.. code-block:: bash

   $ cd doc
   $ open _build/html/index.html
