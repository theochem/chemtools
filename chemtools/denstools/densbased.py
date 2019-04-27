# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
# pragma pylint: disable=invalid-name
"""Density-Based Local Reactivity Tools."""


import numpy as np


__all__ = ['DensityBasedTool']


class DensityBasedTool(object):
    """Class of density-based local descriptive tools."""

    def __init__(self, dens, grad, lap=None, kin=None):
        """Initialize class with density and gradient.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.
        lap : np.ndarray, optional
            Laplacian of electron density evaluated on a set of grid points.
        kin : np.ndarray, optional
            Positive definite kinetic energy density evaluated on a set of grid points.

        """
        if dens.ndim != 1:
            raise ValueError('Argument density should be a 1-dimensional array.')
        if grad.shape != (dens.size, 3):
            raise ValueError('Argument gradient should have same shape as dens array.'
                             ' {0}!={1}'.format(grad.shape, dens.shape))
        if lap is not None and lap.shape != dens.shape:
            raise ValueError('Argument laplacian should have same shape as dens array.'
                             ' {0}!={1}'.format(lap.shape, dens.shape))
        self._dens = dens
        self._grad = grad
        self._lap = lap
        self._kin = kin

    @property
    def density(self):
        r"""Electron density :math:`\rho\left(\mathbf{r}\right)`."""
        return self._dens

    @property
    def gradient(self):
        r"""Gradient of electron density :math:`\nabla \rho\left(\mathbf{r}\right)`.

        This is the first-order partial derivatives of electron density w.r.t. coordinate
        :math:`\mathbf{r} = \left(x\mathbf{i}, y\mathbf{j}, z\mathbf{k}\right)`,

         .. math::
            \nabla\rho\left(\mathbf{r}\right) =
            \left(\frac{\partial}{\partial x}\mathbf{i}, \frac{\partial}{\partial y}\mathbf{j},
                  \frac{\partial}{\partial z}\mathbf{k}\right) \rho\left(\mathbf{r}\right)
        """
        return self._grad

    @property
    def laplacian(self):
        r"""Laplacian of electron density :math:`\nabla ^2 \rho\left(\mathbf{r}\right)`.

        This is defined as the trace of Hessian matrix of electron density which is equal to
        the sum of its :math:`\left(\lambda_1, \lambda_2, \lambda_3\right)` eigen-values:

        .. math::
           \nabla^2 \rho\left(\mathbf{r}\right) = \nabla\cdot\nabla\rho\left(\mathbf{r}\right) =
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial x^2} +
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial y^2} +
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial z^2} =
                     \lambda_1 + \lambda_2 + \lambda_3
        """
        return self._lap

    @property
    def kinetic_energy_density_positive_definite(self):
        r"""Positive definite kinetic energy density.

        .. math::
           \tau \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2
        """
        return self._kin

    @property
    def shannon_information(self):
        r"""Shannon information defined as :math:`\rho(r) \ln \rho(r)`."""
        # masking might be needed
        value = self._dens * np.log(self._dens)
        return value

    @property
    def gradient_norm(self):
        r"""Norm of the gradient of electron density.

        .. math::
           \lvert \nabla \rho\left(\mathbf{r}\right) \rvert = \sqrt{
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial x}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial y}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial z}\right)^2 }
        """
        norm = np.linalg.norm(self._grad, axis=1)
        return norm

    @property
    def reduced_density_gradient(self):
        r"""Reduced density gradient.

        .. math::
           s\left(\mathbf{r}\right) = \frac{1}{2\left(3\pi ^2 \right)^{1/3}}
           \frac{\lvert \nabla\rho\left(\mathbf{r}\right) \rvert}{\rho\left(\mathbf{r}\right)^{4/3}}
        """
        # Mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self._dens, 1.0e-30)
        mdens.filled(1.0e-30)
        # Compute reduced density gradient
        prefactor = 0.5 / (3.0 * np.pi**2)**(1.0 / 3.0)
        rdg = prefactor * self.gradient_norm / mdens**(4.0 / 3.0)
        return rdg

    @property
    def kinetic_energy_density_weizsacker(self):
        r"""Weizsacker kinetic energy density.

        .. math::
           \tau_\text{W} \left(\mathbf{r}\right) =
           \frac{\lvert \nabla\rho\left(\mathbf{r}\right) \rvert^2}{8 \rho\left(\mathbf{r}\right)}
        """
        # Mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self._dens, 1.0e-30)
        mdens.filled(1.0e-30)
        # Compute Weizsacker kinetic energy
        kinetic = self.gradient_norm**2.0 / (8.0 * mdens)
        return kinetic

    @property
    def kinetic_energy_density_thomas_fermi(self):
        r"""Thomas-Fermi kinetic energy density.

        .. math::
           \tau_\text{TF} \left(\mathbf{r}\right) = \frac{3}{10} \left(6 \pi^2 \right)^{2/3}
                  \left(\frac{\rho\left(\mathbf{r}\right)}{2}\right)^{5/3}
        """
        # Compute Thomas-Fermi kinetic energy
        prefactor = 0.3 * (3.0 * np.pi**2.0)**(2.0 / 3.0)
        kinetic = prefactor * self._dens ** (5.0 / 3.0)
        return kinetic
