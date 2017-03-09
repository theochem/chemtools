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

.. _general_energy:

General Energy Model :class:`chemtools.tool.globaltool.GeneralGlobalTool`
=========================================================================

In this model, energy is approximated by an user-specified energy model. Given the
known energy values, this model is parametrized and the energy expression can be evaluated
for any number of electrons.
Being a generic models, this model can reproduce the results of
:ref:`linear <linear_energy>`, :ref:`quadratic <quadratic_energy>`, :ref:`exponential <exponential_energy>`,
and :ref:`rational <rational_energy>` energy models as special cases.

The energy expression should be specified symbolically through `Sympy <http://www.sympy.org/en/index.html>`_.


 .. TODO::
    12. Elaborate more on this model.