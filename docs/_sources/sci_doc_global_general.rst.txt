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

.. _general_energy:

General Energy Model
====================

Because it is tedious to implement all the energy models that have been explored in the literature and
impossible to anticipate all the energy models users may wish to explore, ChemTools allows users to define
their own energy model.
As this is a generic models, this model can reproduce the results of
:ref:`linear <linear_energy>`, :ref:`quadratic <quadratic_energy>`, :ref:`exponential <exponential_energy>`,
and :ref:`rational <rational_energy>` energy models as special cases.

The energy expression should be specified symbolically through `Sympy <http://www.sympy.org/en/index.html>`_.
ChemTools then takes the energies for the user-specified electron numbers and solves the (generally nonlinear)
equations for parameters in the symbolically-defined energy model. The energy model is differentiated symbolically
to determine the fundamental global reactivity indicators in the canonical and grand canonical ensemble.
To define derived reactivity indicators, the maximum number of bound electrons, :math:`N_{\text{max}}`, is determined
by numerical minimization,

 .. math:: N_{\text{max}} = \underbrace {\min }_N E(N)

The derived reactivity indicators are then defined using the usual formulas.
