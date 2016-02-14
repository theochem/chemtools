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


ChemTools Overview
##################

The topics that should be addressed:

#. Why we are developing ChemTools?
#. What are the main features of ChemTools?
#. Where is it heading?
#. Having a word cloud image of the conceptual DFT terms! :-)

Conceptual tools (mostly density-based tools will be considered here, so this is often called “conceptual DFT”) can be divided into three main categories:

* **Global tools:** There is one number for the entire molecule. Examples: Energy, ionization potential.
* **Local tools:** Every point in space, r, has a value. Examples: electron density, electrostatic potential, Fukui function.
* **Nonlocal tools:** There is a value for pairs (or triples, quadruples, etc.) of points. For example, the linear response function measure the change in electron density at r due to a change in external potential at r’.


