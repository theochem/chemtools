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


.. _dftbased_dens:

Electron Density and Its Derivatives
====================================


**Electron Density** :math:`\rho\left(\mathbf{r}\right)` represents the probability of observing an
electron (possibly with specified spin) within a certain volume element :math:`d\boldsymbol{r}` at
:math:`\mathbf{r}`. It is most fundamentally defined through the wave-function:

 .. math:: \rho(\boldsymbol{r}) = N\int \ldots \int \vert
           \Psi(\boldsymbol{x}_1 , \boldsymbol{x}_2 , \ldots , \boldsymbol{x}_N) \vert^2
           d\sigma_1 , d\boldsymbol{x}_2 , \ldots , d\boldsymbol{x}_N ,

where the vectors :math:`\boldsymbol{x}_i` include the space coordinates, :math:`\boldsymbol{r}_i`,
and the spin coordinate, :math:`\sigma_i`, of the i-th electron.


**Shanon Information** :math:`\sigma\left(\mathbf{r}\right)` represents the electron density
per particle:

 .. math::
    \sigma\left(\mathbf{r}\right) = \rho\left(\mathbf{r}\right) \ln \rho\left(\mathbf{r}\right)


**Electron Density Gradient** :math:`\nabla\rho\left(\mathbf{r}\right)` represents the 1st-order
partial derivatives of the electron density with respect to the coordinates which is a vector:

 .. math::
    \nabla \rho\left(\mathbf{r}\right) = \begin{bmatrix}
    \frac{\partial \rho\left(\mathbf{r}\right)}{\partial x} \\
    \frac{\partial \rho\left(\mathbf{r}\right)}{\partial y} \\
    \frac{\partial \rho\left(\mathbf{r}\right)}{\partial z}\end{bmatrix}

The gradient points in the direction of steepest increase of the electron density at a point, and
is zero at a position of a critical point of the electron density.
The :math:`\nabla \rho\left(\mathbf{r}\right)` is therefore a key ingredient for the topological
analysis of the electron density.


**Electron Density Gradient Norm** :math:`\nabla\rho\left(\mathbf{r}\right)` represents the norm
of the electron density gradient vector:

 .. math::
    \lvert \nabla\rho\left(\mathbf{r}\right) \rvert &=
    \sqrt{\nabla\rho\left(\mathbf{r}\right) \cdot \nabla\rho\left(\mathbf{r}\right)} \\ &=
    \sqrt{\left(\frac{\partial \rho\left(\mathbf{r}\right)}{\partial x}\right)^2 +
          \left(\frac{\partial \rho\left(\mathbf{r}\right)}{\partial y}\right)^2 +
          \left(\frac{\partial \rho\left(\mathbf{r}\right)}{\partial z}\right)^2}


**Electron Density Hessian** :math:`\nabla\nabla^{T}\rho\left(\mathbf{r}\right)` represents the
2nd-order partial derivatives of the electron density with respect to coordinates which is a
symmetric matrix:

 .. math::
    \nabla \nabla^{T} \rho\left(\mathbf{r}\right) = \begin{bmatrix}
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

The eigenvalues of the Hessian matrix are often used to classify critical points.


**Electron Density Laplacian** :math:`\nabla^2\rho\left(\mathbf{r}\right)` represents the
trace of the electron density Hessian matrix which is equal to sum of its :math:`\left(\lambda_1,
\lambda_2, \lambda_3\right)` eigenvalues:

 .. math::
    \nabla^2 \rho\left(\mathbf{r}\right) = \nabla\cdot\nabla\rho\left(\mathbf{r}\right) =
    \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial x^2} +
    \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial y^2} +
    \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial z^2} =
    \lambda_1 + \lambda_2 + \lambda_3

It is perhaps the first quantity that was used to visualize shells in molecules. It has been used
and interpreted as a electron-pair-region-locator at least since Richard Bader nailed a picture
of the Laplacian of the electron density to Ron Gillespie’s office door!)
The Laplacian of the electron density tends to be negative where electron density accumulates
and positive elsewhere.

The gradient and Laplacian of the electron density have units of :math:`\text{(length)}^{–4}` and
:math:`\text{(length)}^{–5}`, respectively, and they increase systematically as the number of
electrons increases.
For this reason, one often defines the dimensionless gradient/Laplacian, which are often called the
electron density reduced gradient/Laplacian.
These are often useful for building density functionals in DFT and/or visualizing purposes that are
suitable for both core and valence regions.


**Reduced Density Gradient** :math:`s\left(\mathbf{r}\right)` represents the dimensionless
form of the electron density gradient norm, where the pre-factor is chosen by convention, and is
omitted in some other implementations.

 .. math::
    s\left(\mathbf{r}\right) = \frac{1}{3\left(2\pi ^2 \right)^{1/3}}
    \frac{\lvert \nabla \rho\left(\mathbf{r}\right) \rvert}{\rho\left(\mathbf{r}\right)^{4/3}}

This is used, for example, in the non-covalent interactions analysis (NCI).
The traditional pre-factor one uses is different depending on whether one uses spin-resolved or
spin-unresolved treatments.


**Reduced Density Hessian** :math:`q\left(\mathbf{r}\right)` represents the dimensionless form of
the electron density Laplacian,  where the pre-factor is chosen by convention, and is
omitted in some other implementations.

 .. math::
    q\left(\mathbf{r}\right) =

This is used, for example, in analyzing the shell structure and electron delocalization/metallicity.
The traditional pre-factor one uses is different depending on whether one uses spin-resolved or
spin-unresolved treatments.
