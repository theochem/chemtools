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


.. _dftbased:

Density Functional Theory (DFT) Based
#####################################

The Hohenberg-Kohn theorem assures us that a system’s ground-state electron density contains all
the information necessary to compute all of its observable properties.
Moreover, the mathematical framework of density-functional theory (DFT) provides a rich framework
for deriving chemically intuitive concepts to compute conceptual properties of chemical systems,
even though those properties are rarely observable.

Unlike the framework of conceptual DFT (which is focused on reactivity), the DFT-based descriptive
tools are primarily used to elucidate molecule's electronic structure and binding, though some
of them are relevant to chemical reactivity as well.
At the most fundamental level, the electron density and its derivatives are of the most fundamental
importance.
At the next level, various energetic quantities associated with DFT have a leading role; these
include the electrostatic potential, the kinetic energy density and various approximations thereto,
etc. Suitably chosen linear combinations of these quantities have proven conceptual significance,
as do several of the “intermediate quantities” that are commonly used in developing new density
functionals in DFT.

 .. Note::
    DFT-based tools can be computed in a spin resolved or non-spin resolved manner. Specifically,
    you can initialize the DFTBased class with `spin="a"`, `spin="b"`, `spin="s"`, or `spin="ab"`
    to use :math:`\alpha`, :math:`\beta`, :math:`\alpha - \beta` or :math:`\alpha + \beta` electron
    density, respectively. For simplicity, we do not specify the spin in mathematical formulas that
    follows.


.. toctree::
   :maxdepth: 1

   sci_doc_dftbased_dens
   sci_doc_dftbased_esp
