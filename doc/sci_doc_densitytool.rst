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


.. _density_tools:

Density-Based Local Descriptors
###############################

All the tools for calculating wich the electron density :math:`\rho\left(\mathbf{r}\right)`, gradianet and hessian
of the :math:`N` electron reference state is enough.

**Electron density** :math:`\rho\left(\mathbf{r}\right)` represents ...

**Gradient of electron density** :math:`\nabla \rho\left(\mathbf{r}\right)` represents the first-order partial
derivatives of electron density with respect to coordinates:

 .. math:: \nabla \rho\left(\mathbf{r}\right) =
           \left( \frac{\partial}{\partial x}\mathbf{i}, \frac{\partial}{\partial y}\mathbf{j}, \frac{\partial}{\partial z}\mathbf{k}\right) \rho\left(\mathbf{r}\right)

**Hessian of electron density** :math:`\nabla^2 \rho\left(\mathbf{r}\right)` represents the second-order
partial derivative of electron density with respect to coordinates:


