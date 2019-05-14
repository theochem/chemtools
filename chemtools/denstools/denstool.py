# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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


__all__ = ['DensTool', 'DensGradBasedTool', 'DensGradLapBasedTool', 'DensGradLapKedBasedTool']


class DensTool(object):
    """Local descriptive tools based on density."""

    def __init__(self, dens):
        """Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.

        """
        if dens.ndim != 1:
            raise ValueError('Argument dens should be a 1-dimensional array.')
        self._dens = dens

    @property
    def density(self):
        r"""Electron density :math:`\rho\left(\mathbf{r}\right)`."""
        return self._dens

    @property
    def shannon_information(self):
        r"""Shannon information defined as :math:`\rho(r) \ln \rho(r)`."""
        # TODO: masking might be needed
        value = self.density * np.log(self.density)
        return value

    @property
    def kinetic_energy_density_thomas_fermi(self):
        r"""Thomas-Fermi kinetic energy density.

        .. math::
           \tau_\text{TF} \left(\mathbf{r}\right) = \frac{3}{10} \left(6 \pi^2 \right)^{2/3}
                  \left(\frac{\rho\left(\mathbf{r}\right)}{2}\right)^{5/3}
        """
        # compute Thomas-Fermi kinetic energy
        prefactor = 0.3 * (3.0 * np.pi**2.0)**(2.0 / 3.0)
        kinetic = prefactor * self.density ** (5.0 / 3.0)
        return kinetic


class DensGradBasedTool(DensTool):
    """Local descriptive tools based on density & gradient."""

    def __init__(self, dens, grad):
        """Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.

        """
        super(DensGradBasedTool, self).__init__(dens)
        if grad.shape != (dens.size, 3):
            raise ValueError('Argument grad should be of {0} shape.'.format((dens.shape, 3)))
        self._grad = grad

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
    def gradient_norm(self):
        r"""Norm of the gradient of electron density.

        .. math::
           \lvert \nabla \rho\left(\mathbf{r}\right) \rvert = \sqrt{
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial x}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial y}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial z}\right)^2 }
        """
        norm = np.linalg.norm(self.gradient, axis=1)
        return norm

    @property
    def reduced_density_gradient(self):
        r"""Reduced density gradient.

        .. math::
           s\left(\mathbf{r}\right) = \frac{1}{2\left(3\pi ^2 \right)^{1/3}}
           \frac{\lvert \nabla\rho\left(\mathbf{r}\right) \rvert}{\rho\left(\mathbf{r}\right)^{4/3}}
        """
        # Mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self.density, 1.0e-30)
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
        # mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self.density, 1.0e-30)
        mdens.filled(1.0e-30)
        # compute Weizsacker kinetic energy
        kinetic = self.gradient_norm**2.0 / (8.0 * mdens)
        return kinetic


class DensGradLapBasedTool(DensGradBasedTool):
    """Local descriptive tools based on density, gradient & Laplacian."""

    def __init__(self, dens, grad, lap):
        """Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.
        lap : np.ndarray
            Laplacian of electron density evaluated on a set of grid points.

        """
        super(DensGradLapBasedTool, self).__init__(dens, grad)
        if lap.shape != dens.shape:
            raise ValueError('Argument lap should be of {0} shape.'.format(dens.shape))
        self._lap = lap

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
    def kinetic_energy_density_gradient_expansion(self):
        r"""Gradient expansion approximation of kinetic energy density.

        .. math::
           \tau_\text{GEA} \left(\mathbf{r}\right) =
           \tau_\text{TF} \left(\mathbf{r}\right) +
           \tfrac{1}{9} \tau_\text{W} \left(\mathbf{r}\right) +
           \tfrac{1}{6} \nabla^2 \rho\left(\mathbf{r}\right)
        """
        return self.kinetic_energy_density_gradient_expansion_general(9., 6.)

    @property
    def kinetic_energy_density_gradient_expansion_empirical(self):
        r"""Empirical gradient expansion approximation of kinetic energy density.

        .. math::
           \tau_\text{empGEA} \left(\mathbf{r}\right) =
           \tau_\text{TF} \left(\mathbf{r}\right) +
           \tfrac{1}{5} \tau_\text{W} \left(\mathbf{r}\right) +
           \tfrac{1}{6} \nabla^2 \rho\left(\mathbf{r}\right)
        """
        return self.kinetic_energy_density_gradient_expansion_general(5., 6.)

    def kinetic_energy_density_gradient_expansion_general(self, alpha, beta):
        r"""General gradient expansion approximation of kinetic energy density.

        .. math::
           \tau_\text{genGEA} \left(\mathbf{r}\right) =
           \tau_\text{TF} \left(\mathbf{r}\right) +
           \tfrac{1}{\alpha} \tau_\text{W} \left(\mathbf{r}\right) +
           \tfrac{1}{\beta} \nabla^2 \rho\left(\mathbf{r}\right)

        Parameters
        ----------
        alpha : float
            Value of parameter :math:`\alpha`.
        beta : float
            Value of parameter :math:`\beta`.
        """
        value = self.kinetic_energy_density_thomas_fermi
        value += self.kinetic_energy_density_weizsacker / alpha
        value += self.laplacian / beta
        return value


class DensGradLapKedBasedTool(DensGradLapBasedTool):
    """Local descriptive tools based on density, gradient, Laplacian & kinetic energy density."""

    def __init__(self, dens, grad, lap, ked):
        """Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.
        lap : np.ndarray
            Laplacian of electron density evaluated on a set of grid points.
        ked : np.ndarray
            Kinetic energy density evaluated on a set of grid points.

        """
        super(DensGradLapKedBasedTool, self).__init__(dens, grad)
        if lap.shape != ked.shape:
            raise ValueError('Argument ked should be of {0} shape.'.format(dens.shape))
        self._ked = ked

    @property
    def kinetic_energy_density_positive_definite(self):
        r"""Positive definite kinetic energy density.

        .. math::
           \tau_\text{PD} \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2
        """
        return self._ked

    def kinetic_energy_density_general(self, alpha):
        r"""Return general(ish) kinetic energy density.

        .. math::
           \tau_\text{G} \left(\mathbf{r}, \alpha\right) =
               \tau_\text{PD} \left(\mathbf{r}\right) +
               \tfrac{1}{4} (\alpha - 1) \nabla^2 \rho\left(\mathbf{r}\right)

        Parameters
        ----------
        alpha : float
            Value of parameter :math:`\alpha`.
        """
        return self.kinetic_energy_density_positive_definite + self.laplacian * (alpha - 1) / 4.
