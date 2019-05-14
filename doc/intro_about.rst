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


About ChemTools
###############

ChemTools is a free and open source Python library for interpreting the results of
quantum chemistry calculations. The goal of ChemTools is to provide a toolbox by which
the quantitative output of electronic structure theory calculations can be expressed in
chemical language. ChemTools provides easy-to-use core functionality to compute fundamental
descriptors of conceptual quantum chemistry, together with a flexible set
of utilities allowing scientists to easily test their own discoveries.
ChemTools is designed as a module of the `HORTON package <http://theochem.github.io/horton/>`_,
but can also be used independently to post-process output files of many standard quantum chemistry programs.

Motivated by our interests and our assessment of which portions of the conceptual quantum
chemistry community are most underserved, the current version of ChemTools emphasizes on
conceptual tools associated with, or at least inspired by, density-functional theory (DFT)
and density-matrix theory.
Future developments will include orbital-based tools, information-theoretic methods, and various types
of population and bonding analysis, and quantum chemical topology (including the quantum theory
of atoms in molecules, QTAIM). We also aim to make it easier for theorists to test, implement, and disseminate
new ideas, and to help non-specialists use the most powerful and most recent tools from conceptual
quantum chemistry.

** Main Features **

Conceptual tools (mostly density-based tools will be considered here, so this is often called “conceptual DFT”) can be divided into three main categories:

* **Global tools:** There is one number for the entire molecule.
  Examples: Energy,  ionization potential, electron affinity, chemical potential,
  chemical hardness, chemical softness, hyper-hardness, hyper-softness, electrophilicity, nucleophilicity.
  It provides a quantitative expression for the intrinsic reactivity of reagents in terms of the properties of the isolated molecules.
* **Local tools:** Every point in space, r, has a value. Examples: electron density, electrostatic potential, Fukui function.
* **Nonlocal tools:** There is a value for pairs (or triples, quadruples, etc.) of points. For example, the linear response function measure the change in electron density at r due to a change in external potential at r’.

For canonical ensemble,

 .. image:: ./_static/CDFT.png
     :align: center

For grand-canonical ensemble,

 .. image:: ./_static/CDFT-GCP.png
     :align: center
