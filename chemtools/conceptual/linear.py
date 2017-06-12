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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Linear Energy Model.

This module contains the global and local tool classes corresponding to linear energy models.
"""

from horton import log
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool
from chemtools.utils.utils import doc_inherit

__all__ = ['LinearGlobalTool', 'LinearLocalTool']


class LinearGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the linear energy model.

    The energy is approximated as a piece-wise linear function of the number of electrons,

    .. math:: E(N) = a + b N

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values, the unknown parameters
    of the energy model are obtained by interpolation.

    .. math::
       E(N) = \begin{cases}
               (N_0 - N) \times E(N_0 - 1) + (N - N_0 + 1) \times E(N_0) & N \leqslant N_0 \\
               (N_0 + 1 - N) \times E(N_0) + (N - N_0) \times E(N_0 + 1) & N \geqslant N_0 \\
              \end{cases}

    Because the energy model is not differentiable at integer number of electrons, the first
    derivative of the energy w.r.t. number if electrons is calculated from above, from below
    and averaged:

    .. math::
       \mu^{-} &= -IP \\\
       \mu^{0} &= \frac{\mu^{+} + \mu^{-}}{2} \\
       \mu^{+} &= -EA \\
    """

    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        if energy_zero < energy_plus:
            n_max = n0
        else:
            n_max = None
        super(LinearGlobalTool, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def mu_minus(self):
        r"""
        Chemical potential from below.

        .. math:: \mu^{-} = E\left(N_0\right) - E\left(N_0 - 1\right)  = -IP
        """
        return -1 * self._ip

    @property
    def mu_plus(self):
        r"""
        Chemical potential from above.

        .. math:: \mu^{+} = E\left(N_0 + 1\right) - E\left(N_0\right)  = -EA
        """
        return -1 * self._ea

    @property
    def mu_zero(self):
        r"""
        Chemical potential averaged of :math:`N_0^{+}` and :math:`N_0^{-}`.

        .. math::
           \mu^{0} = \frac{\mu^{+} + \mu^{-}}{2}
                   = \frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2} = - \frac{IP + EA}{2}
        """
        return -0.5 * (self._ea + self._ip)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec is None:
            return None
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        # evaluate energy
        value = self._energy_zero
        if n_elec < self._n0:
            value += (self._n0 - n_elec) * self._ip
        elif n_elec > self._n0:
            value += (self._n0 - n_elec) * self._ea
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec is None:
            return None
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy derivative evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if n_elec == self._n0:
            deriv = None
        elif order >= 2:
            deriv = 0.0
        elif n_elec < self._n0:
            deriv = - self._ip
        elif n_elec > self._n0:
            deriv = - self._ea
        return deriv


class LinearLocalTool(BaseLocalTool):
    r"""
    Class of local conceptual DFT reactivity descriptors based on the linear energy model.

    Considering the interpolated linear energy expression,

    .. math::
       E\left(N\right) =
        \begin{cases}
         \left(N - N_0 + 1\right) E\left(N_0\right) - \left(N - N_0\right) E\left(N_0 - 1\right)
            &  N \leqslant N_0 \\
         \left(N - N_0\right) E\left(N_0 + 1\right) - \left(N - N_0 - 1\right) E\left(N_0\right)
            &  N \geqslant N_0 \\
        \end{cases} \\

    and its derivative with respect to the number of electrons at constant external potential,

    .. math::
       \mu\left(N\right) =
        \begin{cases}
         \mu^- &= E\left(N_0\right) - E\left(N_0 - 1\right) = - IP &&  N < N_0 \\
         \mu^0 &= 0.5 \left(E\left(N_0 + 1\right) - E\left(N_0 - 1\right)\right) = -0.5 (IP + EA)
               && N = N_0 \\
         \mu^+ &= E\left(N_0 + 1\right) - E\left(N_0\right) = - EA &&  N > N_0 \\
        \end{cases}

    the linear local tools are obtained by taking the functional derivative of these expressions
    with respect to external potential :math:`v(\mathbf{r})` at fixed number of electrons.
    """

    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus, n0):
        super(LinearLocalTool, self).__init__(density_zero, density_plus, density_minus, n0)
        self._ff_plus = self._density_plus - self._density_zero
        self._ff_minus = self._density_zero - self._density_minus
        self._ff_zero = 0.5 * (self._density_plus - self._density_minus)

    @property
    def ff_plus(self):
        r"""
        Fukui Function from above, :math:`f^+(\mathbf{r})`.

        .. math::
           f^+\left(\mathbf{r}\right) = \rho_{N_0 + 1}\left(\mathbf{r}\right) -
                                        \rho_{N_0}\left(\mathbf{r}\right)
        """
        return self._ff_plus

    @property
    def ff_minus(self):
        r"""
        Fukui Function from below, :math:`f^-(\mathbf{r})`.

        .. math::
           f^-\left(\mathbf{r}\right) = \rho_{N_0}\left(\mathbf{r}\right) -
                                        \rho_{N_0 - 1}\left(\mathbf{r}\right)
        """
        return self._ff_minus

    @property
    def ff_zero(self):
        r"""
        Fukui Function from center, :math:`f^0(\mathbf{r})`.

        This is defined as the average of :attr:`ff_plus` and :attr:`ff_minus`,

        .. math::
           f^0\left(\mathbf{r}\right) =
           \frac{f^+\left(\mathbf{r}\right) + f^-\left(\mathbf{r}\right)}{2} =
           \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}
        """
        return self._ff_zero

    def density(self, number_electrons=None):
        r"""
        Return linear electron density of :math:`N`-electron system, :math:`\rho_{N}(\mathbf{r})`.

        This is defined as the functional derivative of linear energy model w.r.t.
        external potential at fixed number of electrons, i.e.,

        .. math::
           \rho_{N}(\mathbf{r}) =
            \begin{cases}
             \rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r})
                       \right] \left(N - N_0\right) & \text{ for } N \leqslant N_0 \\
             \rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r})
                       \right] \left(N - N_0\right) & \text{ for } N \geqslant N_0 \\
            \end{cases}

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`\rho_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        rho = self._density_zero.copy()
        if (number_electrons is not None) and (number_electrons != self._n0):
            if number_electrons < self._n0:
                rho += self._ff_minus * (number_electrons - self._n0)
            elif number_electrons > self._n0:
                rho += self._ff_plus * (number_electrons - self._n0)
        return rho

    def fukui_function(self, number_electrons=None):
        r"""
        Return linear Fukui function of :math:`N`-electron system, :math:`f_{N}(\mathbf{r})`.

        This is defined as the functional derivative of linear chemical potential w.r.t.
        external potential at fixed number of electrons,

        .. math::
           f_{N}(\mathbf{r}) =
             \begin{cases}
               f^-(\mathbf{r}) &= \rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r}) && N < N_0 \\
               f^0\left(\mathbf{r}\right) &= 0.5 \left(\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                        \rho_{N_0 - 1}\left(\mathbf{r}\right)\right) && N = N_0 \\
               f^+(\mathbf{r}) &= \rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r}) && N > N_0 \\
             \end{cases}

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`f^0\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._ff_zero
        elif number_electrons < self._n0:
            return self._ff_minus
        elif number_electrons > self._n0:
            return self._ff_plus
        else:
            raise ValueError('Argument number_electrons={0} is not valid.'.format(number_electrons))

    def softness(self, global_softness, number_electrons=None):
        r"""
        Return linear softness of :math:`N`-electron system, :math:`s_N(\mathbf{r})`.

        .. math::
           s_N\left(\mathbf{r}\right) = S \cdot f_N\left(\mathbf{r}\right) =
             \begin{cases}
               S \cdot f^-(\mathbf{r}) & N < N_0 \\
               S \cdot f^0\left(\mathbf{r}\right) & N = N_0 \\
               S \cdot f^+(\mathbf{r}) &  N > N_0 \\
             \end{cases}

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`S \cdot f^0\left(\mathbf{r}\right)`
            is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return global_softness * self._ff_zero
        elif number_electrons < self._n0:
            return global_softness * self._ff_minus
        elif number_electrons > self._n0:
            return global_softness * self._ff_plus
        else:
            raise ValueError('Argument number_electrons={0} is not valid.'.format(number_electrons))
