..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2016-2019 The ChemTools Development Team
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

Supported System
================
ChemTools is available on ``Linux`` and ``MacOS`` with ``Python 2.7``.
In the future release, we will support ``Windows 10`` and ``Python 3.6+``.


Source Code
===========

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
* SymPy >= 0.7.1: http://code.google.com/p/sympy/
* Matplotlib >= 1.0: http://matplotlib.org/
* Nosetests >= 1.1.2: http://readthedocs.org/docs/nose/en/latest/
* HORTON >= 2.0.1: http://theochem.github.io/horton/2.0.1/index.html
* Git-LFS >= 2.0.1: https://git-lfs.github.com/

See :ref:`Documentation Dependencies <usr_doc>` for the dependencies
required for building the documentations.
All dependencies are recommanded to install through Conda.

Conda
~~~~~
For ``Linux`` and ``MacOS`` users, you can download
`miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ to manage the
dependencies. Either version (2.7 or 3.7) will work.

To activate conda base in system:

.. code-block:: bash

    conda activate

To create a virtual environment for ChemTools:

.. code-block:: bash

    conda create -n chemtools python=2.7

To activate ``ChemTools`` virtual environment:

.. code-block:: bash

    conda activate chemtools

To install HORTON and Git-LFS:

.. code-block:: bash

    conda install -c theochem horton
    conda install git-lfs

To deactivate the virtual environment:

.. code-block:: bash

    conda deactivate


Installation
============
Before installing ChemTools, make sure you are in the ``chemtools`` conda
environment. To install ChemTools to your system, run:

.. code-block:: bash

    pip install .

Or, to install ChemTools inplace as an editable package, run:

.. code-block:: bash

    pip install -e .

.. _usr_testing:

Testing
=======

To ensure that all the parts of ChemTools working properly, nosetests can be used to run ChemTool's
automatic tests:

.. code-block:: bash

  nosetests -v chemtools

At this stage, some ``UserWarning`` messages are printed in between tests which is expected.
However, no test should fail.

Uninstallation
==============

To remove ChemTools from your system, run:

.. code-block:: bash

    pip uninstall chemtools

Intro to VMD
============

Go to `VMD <https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD>`_ ,
select the link based on your system OS. Then you will be brought to the registration
page. Shortly after it's done, install the download file to the proper position on your laptop.

For ``macOS`` and ``Linux`` user, you can add the path to your bash profile
(``.bashrc`` on **Linux**, or ``.bash_profile`` on **macOS**) at your home directory
for quick access if it's not on your path.

.. code-block:: bash

    alias vmd="The/path/to/vmd/command/in/your/system"

Usages
~~~~~~

* To visualize in ``VMD``:

.. code-block:: bash

    vmd -e script.vmd

* To save the state, you can click ``File`` > ``Save Visualizaiton State``.

* If you'd like to save the picture during the visualization, you
can click ``File`` > ``Render`` > ``Filename`` > ``Start Rendering``.
