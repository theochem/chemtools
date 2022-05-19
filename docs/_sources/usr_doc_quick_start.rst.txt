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


Quick Start
###########

ChemTools is a package for interpreting the outputs of molecular quantum chemistry calculations.

To use ChemTools, you should first perform an electronic structure calculation on your system(s) of interest
by using your preferred
software package to generate an output file in a format that is supported by ChemTools. Currently supported data
file formats are:

 .. TODO::
    Table of data file formats supported or a link where this information is provided.

If, you wish to work through the examples in this tutorial, make sure that required packages are installed on your computer.
Please see :ref:`Installation <usr_installation>` for instructions.

How to use ChemTools?
=====================

ChemTools is a Python library for post-processing molecular quantum chemistry calculations. So, like any other library,
it can be directly imported and used in Python scripts and codes.
However, this requires knowledge of Python and programming. For users who are not familiar with programming or Python in
particular, a set of built-in Python scripts are provided.
These can be used on a command-line interface and controlled through command-line options and flags.


ChemTools as a Python Library
-----------------------------

If you have programming experience, we strongly advice you to use
ChemTools as an independent Python library. This provides full access to all features available in ChemTools,
and facilitates implementing new functionality based on the existing tools.

The sample snippet below uses ChemTools as a library in a simple Python script. This can be used as a template for
writing scripts to perform computations with ChemTools.

  .. code-block:: python
     :linenos:

     #!/usr/bin/env python

     # import ChemTools library
     import chemtools

     # print version of ChemTools library
     print chemtools.__version__

Commonly the script is saved in a file with ``.py`` extension, line ``script.py``, and run on the command line.
It is more convenient to make the script executable (which only needs to be done once) before running it:

  .. code-block:: bash

     $ chmod +x script.py
     $ ./script.py

Please see :ref:`Examples Gallery <examples gallery>` to learn more about various functionality of ChemTools library


ChemTools as Python Scripts
---------------------------

The built-in tools allow ChemTools to be run on the command-line, without any programming.
We have tried to make these scripts as versatile as possible, however,
not every feature of ChemTools is accessible through scripts.
You can access these functionalities by ``chemtools`` in terminal with ``--help`` or ``-h``. For example:

  .. code-block:: bash

     $ chemtools --help

Currently, the keywards below are available in ChemTools:
* ``mot``: Molecular Orbital Theory (MOT)

* ``nci``: Non-Covalent Interactions (NCI).

* ``elf``: Electron Localization Function (ELF).

* ``lol``: Localized Orbital Locator (LOL).

* ``lcdft``: Local Conceptual DFT.

* ``gcdft``: Global Conceptual DFT.
