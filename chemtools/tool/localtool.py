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
'''Local Conceptual Density Functional Theory (DFT) Reactivity Tools.

   This module contains various local tool classes corresponding to
   linear, quadratic, exponential, and rational energy models.
'''


import numpy as np
from chemtools.utils import doc_inherit


class BaseLocalTool(object):
    '''
    Base class of local conceptual DFT reactivity descriptors.
    '''
    def __init__(self, density_zero, density_plus, density_minus):
        r'''
        Parameters
        ----------
        density_zero : np.ndarray
            Electron density of :math:`N_0`-electron system, i.e. :math:`\rho_{N_0}\left(\mathbf{r}\right)`.
        density_plus : np.ndarray
            Electron density of :math:`(N_0 + 1)`-electron system, i.e. :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)`.
        density_minus : np.ndarray
            Electron density of :math:`(N_0 - 1)`-electron system, i.e. :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`.
        '''
        if np.any(density_zero < 0):
            raise ValueError('Argument density_zero should be all positive!')
        if np.any(density_plus < 0):
            raise ValueError('Argument density_plus should be all positive!')
        if np.any(density_minus < 0):
            raise ValueError('Argument density_minus should be all positive!')
        self._density_zero = density_zero
        self._density_plus = density_plus
        self._density_minus = density_minus

    @property
    def density_zero(self):
        r'''
        Electron density of :math:`N_0`-electron system, i.e.
        :math:`\rho_{N_0}\left(\mathbf{r}\right)`.
        '''
        return self._density_zero

    @property
    def density_plus(self):
        r'''
        Electron density of :math:`(N_0 + 1)`-electron system, i.e.
        :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)`.
        '''
        return self._density_plus

    @property
    def density_minus(self):
        r'''
        Electron density of :math:`(N_0 - 1)`-electron system, i.e.
        :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`.
        '''
        return self._density_minus

    # def __getattr__(self, attr):
    #     '''
    #     '''
    #     # Identify the global property and the type of Fukui Function
    #     global_prop, ff_type = attr.rsplit('_', 1)

    #     # Check for availability of GlobalConceptualTool instance
    #     if self._global_instance is None:
    #         raise ValueError('The argument global_instance is None!')

    #     # Check for availability of global property
    #     if global_prop not in dir(self._global_instance):
    #         raise ValueError('The global property={0} is not known.'.format(global_prop))

    #     # Check for availability of the type of Fukui Function
    #     if ff_type not in ['plus', 'minus', 'zero']:
    #         raise ValueError('The attribute ff_{0} is not known.'.format(ff_type))

    #     # Get the global property & the Fukui Function
    #     global_descriptor = getattr(self._global_instance, global_prop)
    #     fukui_function = getattr(self, 'ff_' + ff_type)
    #     if fukui_function is None:
    #         raise ValueError('The ff_{0} is None!'.format(ff_type))

    #     # Compute the local property
    #     local_descriptor = global_descriptor * fukui_function

    #     return local_descriptor


class LinearLocalTool(BaseLocalTool):
    '''
    Class of local conceptual DFT reactivity descriptors based on the linear energy model.
    '''
    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus):
        super(self.__class__, self).__init__(density_zero, density_plus, density_minus)

    @property
    def ff_plus(self):
        r'''
        Fukui Function from above defined as:

        .. math::
           f^+\left(\mathbf{r}\right) = \rho_{N_0 + 1}\left(\mathbf{r}\right) -
                                        \rho_{N_0}\left(\mathbf{r}\right)
        '''
        ff = self._density_plus - self._density_zero
        return ff

    @property
    def ff_minus(self):
        r'''
        Fukui Function from below define as:

        .. math::
           f^-\left(\mathbf{r}\right) = \rho_{N_0}\left(\mathbf{r}\right) -
                                        \rho_{N_0 - 1}\left(\mathbf{r}\right)
        '''
        ff = self._density_zero - self._density_minus
        return ff

    @property
    def ff_zero(self):
        r'''
        Fukui Function from center defined as the average of :attr:`ff_plus` and :attr:`ff_minus`:

        .. math::
           f^0\left(\mathbf{r}\right) = \frac{f^+\left(\mathbf{r}\right) + f^-\left(\mathbf{r}\right)}{2} =
                    \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}
        '''
        ff = 0.5 * (self._density_plus - self._density_minus)
        return ff

    @property
    def dual_descriptor(self):
        r'''
        Dual descriptor defined as the difference of :attr:`ff_plus` and :attr:`ff_minus`:

        .. math::
           d\left(\mathbf{r}\right) = f^+\left(\mathbf{r}\right) - f^-\left(\mathbf{r}\right) =
            \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0 - 1}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)
        '''
        value = self._density_plus - 2 * self._density_zero + self._density_minus
        return value


class QuadraticLocalTool(BaseLocalTool):
    '''
    Class of local conceptual DFT reactivity descriptors based on the quadratic energy model.
    '''
    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus):
        super(self.__class__, self).__init__(density_zero, density_plus, density_minus)

    @property
    def fukui_function(self):
        r'''
        Fukui function ...

        .. math::

           f\left(\mathbf{r}\right) = f^{(1)}\left(\mathbf{r}\right) &=
                                      \frac{\rho_{N_0+1}\left(\mathbf{r}\right) -
                                      \rho_{N_0-1}\left(\mathbf{r}\right)}{2}
        '''
        ff = 0.5 * (self._density_plus - self._density_minus)
        return ff

    @property
    def hardness(self):
        r'''
        Hardness defined as ...

        .. math::
           \eta\left(\mathbf{r}\right) = f^{(2)}\left(\mathbf{r}\right) =
               \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
               \rho_{N_0 - 1}\left(\mathbf{r}\right)
        '''
        ff2 = self._density_plus - 2 * self._density_zero + self._density_minus
        return ff2

    def softness(self, global_hardness):
        r'''
        Softness defined as ...

        .. math::

           s\left(\mathbf{r}\right) = S \cdot f\left(\mathbf{r}\right) =
                  \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2 \eta} =
                  \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2 \left(IP - EA\right)}

        Parameters
        ----------
        global_hardness : float
            The value of gloabl hardness.
        '''
        s_value = self.fukui_function / global_hardness
        return s_value

    def hyper_softness(self, global_hardness):
        r'''
        Hyper-softness defined as ...

        .. math::

           s^{(2)}\left(\mathbf{r}\right) =
               \frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
               \rho_{N_0 + 1}\left(\mathbf{r}\right)}{\eta^2} =
               \frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
               \rho_{N_0 + 1}\left(\mathbf{r}\right)}{\left(IP - EA\right)^2}

        Parameters
        ----------
        global_hardness : float
            The value of gloabl hardness.
        '''
        s_value = self.hardness / global_hardness**2
        return s_value

    # def compute_fukui_function(self, order):
    #     r'''Fukui Function

    #     .. math::

    #        f^{(0)}\left(\mathbf{r}\right) &= \rho_{N_0}\left(\mathbf{r}\right) \\
    #        f\left(\mathbf{r}\right) = f^{(1)}\left(\mathbf{r}\right) &=
    #                                   \frac{\rho_{N_0+1}\left(\mathbf{r}\right) -
    #                                   \rho_{N_0-1}\left(\mathbf{r}\right)}{2} \\
    #        f^{(2)}\left(\mathbf{r}\right) &= \rho_{N_0 - 1}\left(\mathbf{r}\right) -
    #                                        2 \rho_{N_0}\left(\mathbf{r}\right) +
    #                                          \rho_{N_0 + 1}\left(\mathbf{r}\right) \\
    #        f^{(n)}\left(\mathbf{r}\right) &= 0 \text{ for } n \geq 3
    #     '''
    #     if not (isinstance(order, int) and order >= 0):
    #         raise ValueError('Argument order should be a positive integer.')

    #     if order == 0:
    #         ff = self._density_zero
    #     elif order == 1:
    #         ff = 0.5 * (self._density_plus - self._density_minus)
    #     elif order == 2:
    #         ff = self._density_plus - 2 * self._density_zero + self._density_minus
    #     else:
    #         ff = 0

    #     return ff
