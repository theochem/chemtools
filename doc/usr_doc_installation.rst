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

To install the first seven dependencies (Python related dependencies):

* **Ubuntu Linux 16.04**

  .. code-block:: bash

     sudo apt-get install python-dev python-pip python-numpy python-scipy python-sympy \
                  python-matplotlib python-nose

* **Ubuntu Linux 15.04 & 14.04 and Mac OS**

  .. code-block:: bash

     sudo apt-get install python-dev python-pip python-numpy python-scipy python-sympy \
                  python-matplotlib python-nose
     pip install --user --upgrade numpy scipy sympy matplotlib nose

* **Mac OS (using MacPorts)**

  .. code-block:: bash

     sudo port install python27; sudo port select --set python python27
     sudo port install py27-nose; sudo port select --set nosetests nosetests27
     sudo port install py27-numpy py27-scipy py27-sympy py27-matplotlib
     sudo port install py27-pip; sudo port select --set pip pip27

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

To download the examples files,

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


.. _usr_testing:

Testing
=======

To ensure that all the parts of ChemTools working properly, nosetests can be used to run ChemTool's
automatic tests:

  .. code-block:: bash

     nosetests -v chemtools

At this stage, some ``UserWarning`` messages are printed in between tests which is expected.
However, no test should fail.


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
as a submodule from ChemTools parent directory:

.. code-block:: bash

   git submodule update --init --recursive

Also, make sure that the environment variable ``CTDATA`` is set and
:ref:`examples files are downloaded <usr_lfs_files>`.

To automatically generate API documentation and generate HTML:

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

In case this command did not work, for example on Ubuntu 16.04 you may get a message like **"Couldn't get a
file descriptor referring to the console"**, try:

.. code-block:: bash

   cd doc
   see _build/html/index.html


Quality Assurance
=================

When making a pull request to contribute to the ChemTools repository, the code is remotely tested to see
if it passes all the tests and meets ChemTools' quality standards. To run the tests locally, please refer
to :ref:`Testing <usr_testing>`. If you are interested to run the quality assurance scripts locally, first
install the dependencies below:

* PyLint >= 1.5.0: https://www.pylint.org/
* pycodestyle >= 2.0.0: http://pycodestyle.readthedocs.io/
* pydocstyle >= 1.0.0: http://pydocstyle.readthedocs.io/
* coverage >= 4.1: https://coverage.readthedocs.io/
* Git >= 1.8: https://git-scm.com/
* GitPython >= 2.0.5: http://gitpython.readthedocs.io/

Then, download the quality assurance code by cloning the corresponding submodule:

.. code-block:: bash

   git submodule update --init --recursive

And, run the module's bash script to setup some pre-commit hooks and copy files to run the quality assurance
scripts individually:

.. code-block:: bash

   cd tools/inspector
   ./install.sh
   cd ../..

At this stage, the quality assurance tests can be simulated from the ChemTools parent directory.
For example to run ``pylint`` check,

.. code-block:: bash

   # from ChemTools parent directory
   ./tools/inspector/qa/simulate_trapdoor_pr.py tools/inspector/qa/trapdoor_pylint.py

To run all of the quality assurance scripts,

.. code-block:: bash

   # from ChemTools parent directory
   for i in tools/inspector/qa/trapdoor_*; do tools/inspector/qa/simulate_trapdoor_pr.py $i; done

Note that you should be developing on a feature (not master) branch and merging/rebasing to the
updated master when complete. There should be also no uncommitted changes when running these scripts.
