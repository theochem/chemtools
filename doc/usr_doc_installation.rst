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


`Git Large File Storage (LFS) <https://git-lfs.github.com/>`_ has been used to store the files
and images used in the examples' gallery. In order words, all files in ``doc/examples/*`` are
stored remotely, and only a text pointer to these files exists in the repository.
These files need to be downloaded separately, if you would like to run the example scripts or
make Chemtools' HTML with sphinx. To do so, first you need to install Git LFS using:

  .. code-block:: bash

     $ port install git-lfs
     $ brew install git-lfs

Then go to the repository's directory to download the files by running:

  .. code-block:: bash

     $ git lfs pull

To get a list of all the files tracked with Git LFS, use:

  .. code-block:: bash

     $ git lfs ls-files
