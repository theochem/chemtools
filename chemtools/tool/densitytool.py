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
    def reduced_density_gradient(self):
        r'''
        Reduced density gradient (RDG) defined as:

        .. math::
            \frac{1}{3\left( 3\pi ^2 \right)^{1/3}} \frac{\lvert \nabla \rho \rvert}{\rho^{4/3}}

        '''
        # prefactor and 4/3:
        factor = 2.0*((3.0*(np.pi**2.0))**(1.0/3.0))
        fourtird = 4.0/3.0
        # masking density value's less than 1.0d-30 so we won't devide by 0:
        mdens = np.ma.masked_less(self._density, 1.0e-30)
        mdens.filled(1.0e-30)

        return np.divide(np.linalg.norm(self._gradient, axis=1),(factor*(mdens**fourtird)))

    @property
    def electrostatic_potential(self):
        r'''
        Electrostatic potential defined as:

        .. math::
           \Phi\left(\mathbf{r}\right) = - v \left(\mathbf{r}\right) - \int \frac{\rho\left(\mathbf{r}'\right)}{|\mathbf{r} - \mathbf{r}'|} d \mathbf{r}'
        '''
        # You need v(r)
        pass
