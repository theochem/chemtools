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


.. _tutorial_conceptual:

Tutorial on Conceptual Density Functional Theory
################################################

An extensive set of examples are presented to demonstrate how to compute various descriptors with ChemTools.

Code block

  .. code-block:: python
     :linenos:
     :emphasize-lines: 2,4,7
     :caption: sample.py

     import chemtools as ct

     descriptor = ct.Analyze('h2o.fchk', model='quadratic', approx='FMO')

     # 1st derivative of energy w.r.t. N
     print descriptor.glob.chemical_potential
     print descriptor.glob.mu

     # 2nd derivative of energy w.r.t N
     print descriptor.glob.chemical_hardness
     print descriptor.glob.eta


Code block::

     import chemtools
     descriptor = Analyze('h2o.fchk', model='quadratic', approx='FMO')
     # 1st derivative of energy w.r.t. N
     print descriptor.glob.chemical_potential
     print descriptor.glob.mu
     # 2nd derivative of energy w.r.t N
     print descriptor.glob.chemical_hardness
     print descriptor.glob.eta


Objective: :math:`BF_3` is susceptible to nucleophilic attack at the boron cite.
