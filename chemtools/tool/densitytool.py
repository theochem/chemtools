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
'''Density-Based Local Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


import numpy as np


class DensityLocalTool(object):
    '''
    Class of desnity-based local descriptive tools.
    '''
    def __init__(self, density, gradient, hessian=None):
        '''
        Parameters
        ----------
        density : np.ndarray
            Electron density of the system evaluated on a grid
        gradient : np.ndarray
            Gradient vector of electron density evaluated on a grid
        hessian : np.ndarray
            Hessian matrix of electron density evaluated on a grid
        '''
        if density.ndim != 1:
            raise ValueError('Argument desnity should be a 1-dimensioanl array.')
        if gradient.shape != (density.size, 3):
            raise ValueError('Argument gradient should have same shape as density arrary. {0}!={1}'.format(gradient.shape, density.shape))
        if (hessian is not None) and hessian.shape != (density.size, 3, 3):
                raise ValueError('Argument hessian\'s shape is not consistent with the density array. {0}!={1}'.format(hessian.shape, (density.size, 3, 3)))

        self._density = density
        self._gradient = gradient
        self._hessian = hessian

    @property
    def density(self):
        r'''
        Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on a grid.
        '''
        return self._density

    @property
    def gradient(self):
        r'''
        Gradient vector of electron :math:`\nabla \rho\left(\mathbf{r}\right)`
        defined as the first-order partial derivatives of electron density w.r.t. coordinate
        :math:`\mathbf{r} = \left(x\mathbf{i}, y\mathbf{j}, z\mathbf{k}\right)`:

         .. math::

            \nabla\rho\left(\mathbf{r}\right) =
            \left(\frac{\partial}{\partial x}\mathbf{i}, \frac{\partial}{\partial y}\mathbf{j},
                  \frac{\partial}{\partial z}\mathbf{k}\right) \rho\left(\mathbf{r}\right)
        '''
        return self._gradient

    @property
    def hessian(self):
        r'''
        Hessian matrix of electron density :math:`\nabla^2 \rho\left(\mathbf{r}\right)`
        defined as the second-order partial derivatives of electron density w.r.t coordinate
        :math:`\mathbf{r} = \left(x\mathbf{i}, y\mathbf{j}, z\mathbf{k}\right)`:
        '''
        return self._hessian

    @property
    def laplacian(self):
        r'''
        Laplacian of electron density :math:`\nabla ^2 \rho\left(\mathbf{r}\right)` defined
        as the trace of Hessian matrix of electron desnity which is equal to the sum of
        :math:`\left(\lambda_1, \lambda_2, \lambda_3\right)` eigen-values of Hessian matrix:

        .. math::
           \nabla^2 \rho\left(\mathbf{r}\right) = \nabla\cdot\nabla\rho\left(\mathbf{r}\right) =
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial x^2} +
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial y^2} +
                     \frac{\partial^2\rho\left(\mathbf{r}\right)}{\partial z^2} =
                     \lambda_1 + \lambda_2 + \lambda_3
        '''
        # This is not a local tool!
        pass

    @property
    def shanon_information(self):
        r'''
        Shanon information defined as :math:`\rho(r) \ln \rho(r)`.
        '''
        # masking might be needed
        value = self._density * np.log(self._density)
        return value

    @property
    def gradient_norm(self):
        r'''
        Gradinet norm representing the norm of the gradient vector at every point:

        .. math::
           \lvert \nabla \rho\left(\mathbf{r}\right) \rvert = \sqrt{
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial x}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial y}\right)^2 +
                  \left(\frac{\partial\rho\left(\mathbf{r}\right)}{\partial z}\right)^2 }
        '''
        norm = np.linalg.norm(self._gradient, axis=1)
        return norm

    @property
    def reduced_density_gradient(self):
        r'''
        Reduced density gradient (RDG) defined as:

        .. math::
           s\left(\mathbf{r}\right) = \frac{1}{3\left(2\pi ^2 \right)^{1/3}}
           \frac{\lvert \nabla \rho\left(\mathbf{r}\right) \rvert}{\rho\left(\mathbf{r}\right)^{4/3}}
        '''
        # Mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self._density, 1.0e-30)
        mdens.filled(1.0e-30)
        # Compute reduced density gradient
        prefactor = 0.5 / (3.0 * np.pi**2)**(1.0/3.0)
        rdg = prefactor * self.gradient_norm / mdens**(4.0/3.0)
        return rdg

    @property
    def weizsacker_kinetic_energy_density(self):
        r'''
        Weizsacker kinetic energy/local steric energy/Fisher information density defined as:

        .. math::
           T\left(\mathbf{r}\right) =
           \frac{\lvert \nabla \rho\left(\mathbf{r}\right) \rvert ^2}{8 \rho\left(\mathbf{r}\right)}
        '''
        # Mask density values less than 1.0d-30 to avoid diving by zero
        mdens = np.ma.masked_less(self._density, 1.0e-30)
        mdens.filled(1.0e-30)
        # Compute Weizsacker kinetic energy
        wke = self.gradient_norm**2.0 / (8.0 * mdens)
        return wke

    @property
    def thomas_fermi_kinetic_energy_density(self):
        r'''
            Thomas-Fermi kinetic energy density defined as:

        .. math::
            T\left(\mathbf{r}\right) = \frac{3}{10} \left( 6 \pi ^2 \right)^{2/3}
            \left( \frac{\rho\left(\mathbf{r}\right)}{2} \right)^{5/3}
            '''
        # Compute Thomas-Fermi kinetic energy
        prefactor = 0.3 * (3.0 * (np.pi**2.0))**(2.0/3.0)
        fivethird = 5.0 / 3.0
        tfke =  prefactor * (self._density**fivethird)
        return tfke


    def electrostatic_potential(self, numbers, coordinates, int_weights, int_points, points):
        r'''
        Electrostatic potential defined as:

        .. math::
           \Phi\left(\mathbf{r}\right) = - v \left(\mathbf{r}\right) - \int \frac{\rho\left(\mathbf{r}'\right)}{|\mathbf{r} - \mathbf{r}'|} d \mathbf{r}'

        Parameters
        ----------
        numbers : np.ndarray
            The atomic numbers of the system
        coordinates : np.ndarray
            The coordinates of the (nuclear) charges of the system
        int_weights : np.ndarray
            The integration weights.
            This should have the same dimension as the density used to initialize the class.
        int_points : np.ndarray
            The coordinates of the integration points.
        points : np.ndarray
            The coordinates of the point(s) on which to calculate the electrostatic potential.
        '''
        # check consistency of arrays
        if len(coordinates) != len(numbers):
            raise ValueError('Argument numbers & coordinates should have the same length. {0}!={1}'.format(len(coordinates), len(numbers)))
        if len(int_weights) != len(self._density):
            raise ValueError('Argument int_weights & density should have the same shape. {0}!={1}'.format(int_weights.shape, self._density.shape))
        if len(int_points) != len(self._density):
            raise ValueError('Argument int_points & density should have the same shape. {0}!={1}'.format(int_weights.shape, self._density.shape))

        # array to store esp
        esp = np.zeros(points.shape[0])

        # compute esp at every point
        for n, point in enumerate(points):
            # esp at a given point from nuclei
            for index, number in enumerate(numbers):
                delta = point - coordinates[index, :]
                distance = np.linalg.norm(delta)
                esp[n] += number / distance

            # esp at a given point from electron density multiplied by integration weight
            for index, density in enumerate(self._density):
                delta = point - int_points[index, :]
                distance = np.linalg.norm(delta)
                # avoid computing esp, if points are very close
                if distance > 1.e-6:
                   esp[n] -= density * int_weights[index] / distance

        return esp
