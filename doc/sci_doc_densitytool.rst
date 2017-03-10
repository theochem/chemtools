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

All the tools for calculations based on the electron density :math:`\rho\left(\mathbf{r}\right)`,
gradiant and hessian of the :math:`N` electron reference state.

**Electron density** :math:`\rho\left(\mathbf{r}\right)` represents the probability of finding an electron within a certain volume element :math:`d\boldsymbol{r}`. It is found
by integrating the square of the wave function over all variables except the spatial coordinates of one electron:

 .. math:: \rho(\boldsymbol{r}) = N\int \ldots \int \vert 
           \Psi(\boldsymbol{x}_1 , \boldsymbol{x}_2 , \ldots , \boldsymbol{x}_N) \vert^2 
           d\sigma_1 , d\boldsymbol{x}_2 , \ldots , d\boldsymbol{x}_N ,

where the vectors :math:`\boldsymbol{x}_i` include the space coordinates, :math:`\boldsymbol{r}_i`, and the spin coordinate, :math:`\sigma_i`, of the i-th electron.

**Gradient of electron density** :math:`\nabla \rho\left(\mathbf{r}\right)` represents the first-order partial
derivatives of the electron density with respect tothe  coordinates:

 .. math:: \nabla \rho\left(\mathbf{r}\right) =
           \left( \frac{\partial}{\partial x} \rho\left(\mathbf{r}\right), \frac{\partial}{\partial y} \rho\left(\mathbf{r}\right), \frac{\partial}{\partial z} \rho\left(\mathbf{r}\right)\right) 

**Hessian of electron density** :math:`\nabla^2 \rho\left(\mathbf{r}\right)` represents the second-order
partial derivative of the electron density with respect to coordinates:

 .. math::  \nabla^2 \rho\left(\mathbf{r}\right) = \begin{bmatrix}
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial x^2}  & 
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial x \partial y}  &
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial x \partial z}  \\
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial y \partial x}  & 
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial y^2}  &
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial y \partial z}  \\
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial z \partial x}  & 
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial z \partial y}  &
                \frac{\partial^2 \rho\left(\mathbf{r}\right)}{\partial z^2}  \\
            \end{bmatrix}


