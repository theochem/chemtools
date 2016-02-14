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


Local Descriptive Tools
#######################

 .. image:: _static/under_construction.gif
    :scale: 75 %
    :align: center

Local descriptive tools assign a value to every point in space. In conceptual DFT,
local descriptors arise as functional derivatives of global descriptors with respect
to local quantities, typically the external potential. Ergo,

 .. math::

    \rho(r)    &= {\left( \frac{\delta E}{\delta v(r)} \right)_N} \\
    f(r)       &= {\left( \frac{\delta}{\delta v(r)} {\left( \frac{\partial E}{\partial N} \right)_{v(r)}} \right)_N} \\
               &= {\left( \frac{\delta \mu}{\delta v(r)} \right)_N} \\
    f^{(2)}(r) &= {\left( \frac{\delta}{\delta v(r)} {\left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)}} \right)_N} \\
               &= {\left( \frac{\delta \eta}{\delta v(r)} \right)_N} \\
    f^{(n)}(r) &= {\left( \frac{\delta}{\delta v(r)} {\left( \frac{\partial^n E}{\partial N^n} \right)_{v(r)}} \right)_N} \\
               &= {\left( \frac{\delta \eta^{(n-1)}}{\delta v(r)} \right)_N}

