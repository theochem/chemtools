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


.. _usr_py_depend:

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
* :ref:`Git LFS <usr_lfs_installation>` >= 2.0.1: https://git-lfs.github.com/

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

To install HORTON, **Linux** users can follow the instructions from `Download and Install
<http://theochem.github.io/horton/2.0.1/user_download_and_install_linux.html>`_ of HORTON's documentation.
**MacOS** users can install HORTON with MacPort directly:

.. code-block:: bash

   sudo port install horton


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

  You can install LFS with MacPort

  .. code-block:: bash

     sudo port install git-lfs

  Or Homebrew

  .. code-block:: bash

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


.. _usr_testing:

Testing
=======

To ensure that all the parts of ChemTools working properly, nosetests can be used to run ChemTool's
automatic tests:

  .. code-block:: bash

     nosetests -v chemtools

At this stage, some ``UserWarning`` messages are printed in between tests which is expected.
However, no test should fail.
