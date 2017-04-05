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

     git clone https://github.com/QuantumElephant/chemtools.git chemtools


Dependencies
============

The following dependencies will be necessary for ChemTools to build properly,

* Python >= 2.7, < 3.0: http://www.python.org/ (Also install development files.)
* PIP >= 7.0: https://pip.pypa.io/ (Not required in some scenarios but never bad to have.)
* SciPy >= 0.11.0: http://www.scipy.org/
* NumPy >= 1.9.1: http://www.numpy.org/
* Sympy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* HORTON >= 2.0.1: http://theochem.github.io/horton/2.0.1/index.html
* :ref:`Git Large File Storage (LFS) <usr_lfs_installation>`

See :ref:`Documentation Dependencies <usr_doc>` for the dependencies
required for building the documentations.

Python Dependencies
~~~~~~~~~~~~~~~~~~~

To install the first seven dependecies (Python related dependencies):

* **Ubuntu Linux 16.04**

  .. code-block:: bash

     sudo apt-get install python-dev python-pip python-numpy python-scipy python-sympy \
                  python-matplotlib python-nose

* **Ubuntu Linux 15.04 & 14.04 and Mac OS**

  .. code-block:: bash

     sudo apt-get install python-dev python-pip python-numpy python-scipy python-sympy \
                  python-matplotlib python-nose
     pip install --user --upgrade numpy scipy sympy matplotlib nose

For the latest versions of NumPy, SciPy, SymPy, Matplotlib and Nosetests, ``pip`` can be used to
upgrade existing packages, as shown above.

HORTON
~~~~~~

To install HORTON, follow the instructions from `Download and Install 
<http://theochem.github.io/horton/2.0.1/user_download_and_install.html>`_ of HORTON's documentation.

.. _usr_lfs_installation:

Git LFS
~~~~~~~

`Git Large File Storage (LFS) <https://git-lfs.github.com/>`_ is used to store files that are not
part of the ChemTools code, but are used in some ways, such as generating the examples in the 
documentation.
These files need to be downloaded separately, for example, if you would like to run the example
scripts, go through tutorials (using exactly the same files used) or make Chemtools HTML with
sphinx.

To install Git LFS,

* **Mac OS**

  .. code-block:: bash

     sudo port install git-lfs
     brew install git-lfs

* **Linux OS**

  .. code-block:: bash

     cd your_download_directory
     wget https://github.com/git-lfs/git-lfs/releases/download/v2.0.1/git-lfs-linux-amd64-2.0.1.tar.gz
     tar -zxvf git-lfs-linux-amd64-2.0.1.tar.gz
     cd git-lfs-2.0.1
     ./install.sh


.. _usr_lfs_files:

To download the files,

  .. code-block:: bash

     git lfs pull


To get a list of all the files tracked with Git LFS,

  .. code-block:: bash

     git lfs ls-files


Installation
============

To install ChemTools run:

  .. code-block:: bash

     ./setup.py install --user

If you prefer to run ChemTools from the source folder, ``PYTHONPATH`` and ``CTDATA`` paths
need to be add into your **~/.bashrc** (Linux) or **~/.profile** (MacOS)

  .. code-block:: bash

     export PYTHONPATH=$PYTHONPATH:{path_to_chemtools_repo}/chemtools
     export CTDATA={path_to_chemtools_repo}/chemtools/data

Testing
=======

To ensure that all the parts of ChemTools working properly, nosetests can be used to run ChemTool's
unit tests:

  .. code-block:: bash

     nosetests -v chemtools


.. _usr_doc:

Documentation
=============

If you are interested in generating the documentation from source, the following
packages are also needed:

* Sphinx >=1.3.1: http://sphinx.pocoo.org/
* Sphinx Read-The-Docs theme >=0.1.8: https://github.com/snide/sphinx_rtd_theme
* sphinxcontrib-bibtex >= 0.3.5: https://pypi.python.org/pypi/sphinxcontrib-bibtex
* IPython >= 3.2.1: https://ipython.org/install.html

To install these dependencies,

* **Ubuntu Linux 16.04**

  .. code-block:: bash

     sudo apt-get install python-sphinx python-sphinx-rtd-theme sphinxcontrib-bibtex ipython

* **Ubuntu Linux 15.04 & 14.04 and Mac OS**

  .. code-block:: bash

     pip install --user --upgrade sphinx sphinx_rtd_theme sphinxcontrib-bibtex ipython

The Sphinx Read-The-Docs theme customized for ChemTools can be obtained cloning the repository
as a submodule:

.. code-block:: bash

   git submodule update --init --recursive

To automatically generate API documentation and generate HTML (this requires ``data/examples``
files; to obtain them please refer to :ref:`usr_lfs_files`) use the commands below:

.. code-block:: bash

   cd doc
   sphinx-apidoc -f -o source ../
   make clean
   make html

To open the documentation in your default browser, either click on ``doc/_build/html/index.html``
file directly, or run the command below from terminal:

.. code-block:: bash

   cd doc
   open _build/html/index.html


Quality Assurance
=================
When contributing to the ChemTools repo, the code is remotely tested to see if it meets ChemTools'
standards. To run these tests locally, you must first download and install the quality assurance
code. From the ChemTools main directory,

.. code-block:: bash

   git submodule update --init --recursive
   cd ./tools/inspector
   ./install.sh
   cd ../..

Then, the quality assurance tests can be simulated from the ChemTools main directory with

.. code-block:: bash

   ./tools/inspector/qa/simulate_trapdoor_pr.py

Note that you should be developing on a non-master (feature) branch and merging/rebasing to the
updated master when complete.
