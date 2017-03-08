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


.. _tutorial_nci:

Tutorial on Non Covalent Interactions (NCI)
###########################################

The easiest way to calculate the Non Covalent Interactions, by using the default settings is as follows:

Code block::

     import chemtools
     # Build the NCI tool
     nci = NCI.from_file(('h2o_dimer_pbe_sto3g.fchk')
     # dump the files needed for the visualisation by giving a filename for the cube files
     nci.dump_files('h2o_dimer')

this generates three files:
  - h2o_dimer_pbe_sto3g-dens.cube (a density cube file)
  - h2o_dimer_pbe_sto3g-grad.cube (a reduced density gradient cube file)
  - h2o_dimer_pbe_sto3g.vmd  (a vmd script file)


