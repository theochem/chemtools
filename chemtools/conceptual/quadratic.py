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


from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons


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

    def __init__(self, dict_energy):
        r"""Initialize quadratic energy model to compute global reactivity descriptors.

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
        # calculate parameters a, b, c of quadratic energy model
        energy_m, energy_0, energy_p = [dict_energy[n] for n in sorted(dict_energy.keys())]
        param_c = 0.5 * (energy_m - 2 * energy_0 + energy_p)
        param_b = 0.5 * (energy_p - energy_m) - 2 * param_c * n_ref
        param_a = energy_0 - param_b * n_ref - param_c * (n_ref**2)
        self._params = [param_a, param_b, param_c]
        # calculate N_max (number of electrons for which energy is minimum)
        n_max = - param_b / (2 * param_c)
        super(QuadraticGlobalTool, self).__init__(n_ref, n_max)
        self.dict_energy = dict_energy

    @property
    def params(self):
        """Parameter :math:`a`, :math:`b` and :math:`c` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # evaluate energy
        value = self._params[0] + self._params[1] * n_elec + self._params[2] * n_elec**2
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")
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
    def __init__(self, dict_density):
        # check number of electrons & density values
        n_ref, dens_m, dens_0, dens_p = check_dict_values(dict_density)
        # compute fukui function & dual descriptor of N0-electron system
        self._ff0 = 0.5 * (dens_p - dens_m)
        self._df0 = dens_p - 2 * dens_0 + dens_m
        super(QuadraticLocalTool, self).__init__(dict_density, n_ref)

    def density(self, n_elec):
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
        n_elec : float
            Number of electrons.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute density
        rho = self.density_zero.copy()
        if n_elec != self._n0:
            dN = n_elec - self._n0
            rho += self._ff0 * dN + 0.5 * self._df0 * dN**2
        return rho

    def fukui_function(self, n_elec):
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
        n_elec : float
            Number of electrons.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute fukui function
        if n_elec == self._n0:
            return self._ff0
        else:
            ff = self._ff0 + self._df0 * (n_elec - self.n0)
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

    def softness(self, n_elec, global_softness):
        r"""
        Return quadratic softness of :math:`N`-electron system, :math:`s_{N}(\mathbf{r})`.

        .. math::
           s_N\left(\mathbf{r}\right) = S \cdot f_N\left(\mathbf{r}\right)

        Parameters
        ----------
        n_elec : float
            Number of electrons.
        global_softness : float
            The value of global softness.
        """
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute softness
        s_value = global_softness * self.fukui_function(n_elec)
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
