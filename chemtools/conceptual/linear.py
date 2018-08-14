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
from chemtools.conceptual.utils import check_dict_values, check_number_electrons
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

    def __init__(self, dict_energy):
        r"""Initialize linear energy model to compute global reactivity descriptors.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): E(N_0 - 1), N_0: E(N_0),
            (N_0 + 1): E(N_0 + 1)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.

        """
        # check number of electrons & energy values
        n_ref, energy_m, energy_0, energy_p = check_dict_values(dict_energy)
        # calculate parameters a, b, a' and b' of linear energy model
        param_b = energy_0 - energy_m
        param_a = energy_0 - n_ref * param_b
        param_b_prime = energy_p - energy_0
        param_a_prime = energy_0 - n_ref * param_b_prime
        self._params = [param_a, param_b, param_a_prime, param_b_prime]
        # calculate N_max
        if energy_0 < energy_p:
            n_max = n_ref
        else:
            n_max = None
        super(LinearGlobalTool, self).__init__(n_ref, n_max)
        self.dict_energy = dict_energy

    @property
    def params(self):
        r"""Parameter :math:`a`, :math:`b`, :math:`a^\prime` & :math:`b^\prime` of energy model."""
        return self._params

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
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # evaluate energy
        if n_elec <= self._n0:
            value = self._params[0] + n_elec * self._params[1]
        elif n_elec > self._n0:
            value = self._params[2] + n_elec * self._params[3]
        else:
            raise ValueError("Argument n_elec should be a positive real number!")
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")
        # evaluate derivative
        if n_elec == self._n0:
            deriv = None
        elif order >= 2:
            deriv = 0.0
        elif n_elec < self._n0:
            deriv = self._params[1]
        elif n_elec > self._n0:
            deriv = self._params[3]
        else:
            raise ValueError("Argument n_elec should be a positive real number!")
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
    def __init__(self, dict_density):
        # check number of electrons & density values
        n_ref, dens_m, dens_0, dens_p = check_dict_values(dict_density)
        # compute ff+, ff- & ff0
        self._ff_plus = dens_p - dens_0
        self._ff_minus = dens_0 - dens_m
        self._ff_zero = 0.5 * (dens_p - dens_m)
        super(LinearLocalTool, self).__init__(dict_density, n_ref)

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

    def density(self, n_elec):
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
        n_elec : float
            Number of electrons.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute density
        rho = self.density_zero.copy()
        if n_elec != self._n0:
            if n_elec < self._n0:
                rho += self._ff_minus * (n_elec - self._n0)
            elif n_elec > self._n0:
                rho += self._ff_plus * (n_elec - self._n0)
        return rho

    def fukui_function(self, n_elec):
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
        n_elec : float
            Number of electrons.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute fukui function
        if n_elec == self._n0:
            ff = self._ff_zero
        elif n_elec < self._n0:
            ff = self._ff_minus
        elif n_elec > self._n0:
            ff = self._ff_plus
        return ff

    def softness(self, n_elec, global_softness):
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
        n_elec : float
            Number of electrons.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute softness
        if n_elec == self._n0:
            softness = global_softness * self._ff_zero
        elif n_elec < self._n0:
            softness = global_softness * self._ff_minus
        elif n_elec > self._n0:
            softness = global_softness * self._ff_plus
        return softness
