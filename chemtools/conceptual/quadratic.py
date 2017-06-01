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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Quadratic Energy Model.

This module contains the global and local tool classes corresponding to quadratic energy models.
"""

from horton import log
from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool

__all__ = ['QuadraticGlobalTool', 'QuadraticLocalTool']


class QuadraticGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the quadratic energy model.

    The energy is approximated as a quadratic function of the number of electrons,

    .. math:: E(N) = a + b N + c N^2

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values, the unknown parameters
    of the energy model are obtained by interpolation.

    First, second and higher order derivatives of the quadratic energy model with respect to
    the number of electrons at fixed external potential are given by:

    .. math::
       \left(\frac{\partial E}{\partial N}\right)_{v(\mathbf{r})} &= b + 2 c N \\
       \left(\frac{\partial^2 E}{\partial N^2}\right)_{v(\mathbf{r})} &= 2 c \\
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} &= 0
             \quad \text{for} \quad n \geq 2
    """

    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # calculate parameters a, b, c of quadratic energy model
        c = 0.5 * (energy_minus - 2 * energy_zero + energy_plus)
        b = 0.5 * (energy_plus - energy_minus) - 2 * c * n0
        a = energy_zero - b * n0 - c * (n0**2)
        self._params = [a, b, c]
        # calculate Nmax (number of electrons for which energy is minimum)
        n_max = - b / (2 * c)
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        """Parameter :math:`a`, :math:`b` and :math:`c` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        # evaluate energy
        value = self._params[0] + self._params[1] * n_elec + self._params[2] * n_elec**2
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy derivative evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if order == 1:
            deriv = self._params[1] + 2 * n_elec * self._params[2]
        elif order == 2:
            deriv = 2 * self._params[2]
        elif order >= 2:
            deriv = 0.0
        return deriv


class QuadraticLocalTool(BaseLocalTool):
    r"""
    Class of local conceptual DFT reactivity descriptors based on the quadratic energy model.

    Considering the interpolated quadratic energy expression,

    .. math::
       E\left(N\right) = E\left(N_0\right) &+ \left(\frac{E\left(N_0 + 1\right) -
                         E\left(N_0 - 1\right)}{2}\right) \left(N - N_0\right) \\
                &+ \left(\frac{E\left(N_0 - 1\right) - 2 E\left(N_0\right) +
                   E\left(N_0 + 1\right)}{2}\right) \left(N - N_0\right)^2 \\

    and its first and second derivatives with respect to the number of electrons at constant
    external potential,

    .. math::
       \mu\left(N\right) &= \left(\frac{\partial E\left(N\right)}{\partial N}
                            \right)_{v(\mathbf{r})} \\
         &= \left(\frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2}\right) +
            \left[E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)
            \right] \left(N - N_0\right) \\
       \eta\left(N\right) &= \left(\frac{\partial^2 E\left(N\right)}{\partial^2 N}
                             \right)_{v(\mathbf{r})}
         = E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)

    the quadratic local tools are obtained by taking the functional derivative of these expressions
    with respect to external potential :math:`v(\mathbf{r})` at fixed number of electrons.
    """

    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus, n0):
        super(self.__class__, self).__init__(density_zero, density_plus, density_minus, n0)
        # Fukui function and dual descriptor of N0-electron system
        self._ff0 = 0.5 * (self._density_plus - self._density_minus)
        self._df0 = self._density_plus - 2 * self._density_zero + self._density_minus

    def density(self, number_electrons=None):
        r"""
        Return quadratic electron density of :math:`N`-electron system :math:`\rho_{N}(\mathbf{r})`.

        This is defined as the functional derivative of quadratic energy model w.r.t.
        external potential at fixed number of electrons,

        .. math::
           \rho_{N}(\mathbf{r}) = \rho_{N_0}\left(\mathbf{r}\right)
             &+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right) \\
             &+ \left(\frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
                \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)}{2}\right)
                \left(N - N_0\right)^2

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`\rho_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._density_zero
        else:
            dN = (number_electrons - self._n0)
            rho = self._density_zero + self._ff0 * dN + 0.5 * self._df0 * dN**2
            return rho

    def fukui_function(self, number_electrons=None):
        r"""
        Return quadratic Fukui function of :math:`N`-electron system, :math:`f_{N}(\mathbf{r})`.

        This is defined as the functional derivative of quadratic chemical potential w.r.t.
        external potential at fixed number of electrons,

        .. math::
           f_{N}(\mathbf{r}) = \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                 \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2} \right) +
                 \left[\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
                 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)
                 \right] \left(N - N_0\right)

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`f_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._ff0
        else:
            ff = self._ff0 + self._df0 * (number_electrons - self.n0)
            return ff

    def dual_descriptor(self):
        r"""
        Quadratic dual descriptor of :math:`N`-electron system, :math:`\Delta f_{N}(\mathbf{r})`.

        This is defined as the functional derivative of quadratic chemical hardness
        w.r.t. external potential at fixed number of electrons,

        .. math::
           \Delta f_{N}(\mathbf{r}) = \rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
            \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)

        The quadratic dual descriptor is independent of the number electrons.
        """
        return self._df0

    def softness(self, global_softness, number_electrons=None):
        r"""
        Return quadratic softness of :math:`N`-electron system, :math:`s_{N}(\mathbf{r})`.

        .. math::
           s_N\left(\mathbf{r}\right) = S \cdot f_N\left(\mathbf{r}\right)

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`s_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        s_value = global_softness * self.fukui_function(number_electrons)
        return s_value

    def hyper_softness(self, global_softness):
        r"""
        Quadratic hyper-softness of :math:`N`-electron system, :math:`s_N^{(2)}(\mathbf{r})`.

        .. math::
           s_N^{(2)}\left(\mathbf{r}\right) = S^2 \cdot \Delta f_N\left(\mathbf{r}\right)

        The quadratic hyper-softness is independent of the number electrons.

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        """
        s_value = global_softness**2 * self.dual_descriptor()
        return s_value
