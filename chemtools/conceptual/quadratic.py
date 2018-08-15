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
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool, BaseCondensedTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons


__all__ = ["QuadraticGlobalTool", "QuadraticLocalTool", "QuadraticCondensedTool"]


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
        else:
            deriv = 0.
        return deriv


class QuadraticLocalTool(BaseLocalTool):
    r"""
    Class of local conceptual DFT reactivity descriptors based on the quadratic energy model.

    Considering the interpolated :class:`quadratic energy model <QuadraticGlobalTool>` and its
    derivatives, the quadratic local tools are obtained by taking the functional derivative
    of these expressions with respect to external potential :math:`v(\mathbf{r})` at fixed
    number of electrons :math:`N`.

    Given the electron density corresponding to energy values used for interpolating the
    energy model, i.e., :math:`\rho_{N_0 - 1}(\mathbf{r})`, :math:`\rho_{N_0}(\mathbf{r})`
    and :math:`\rho_{N_0 + 1}(\mathbf{r})`, the :func:`density <QuadraticLocalTool.density>`
    of the :math:`N` electron system :math:`\rho_{N}(\mathbf{r})` is given by:

    .. math::
       \rho_{N}(\mathbf{r}) = \rho_{N_0}\left(\mathbf{r}\right)
         &+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) -
            \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right) \\
         &+ \left(\frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
            \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)}{2}\right)
            \left(N - N_0\right)^2

    The :func:`density derivative <QuadraticLocalTool.density_derivative>` with respect to the
    number of electrons at fixed external potential is given by:

    .. math::
      \left(\frac{\partial \rho_N(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} &=
      \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}
      \right) + \left[\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
      \rho_{N_0 - 1}\left(\mathbf{r}\right) \right] \left(N - N_0\right) \\
      \left(\frac{\partial^2 \rho_N(\mathbf{r})}{\partial N^2}\right)_{v(\mathbf{r})} &=
      \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
      \rho_{N_0 - 1}\left(\mathbf{r}\right) \\
      \left(\frac{\partial^n \rho_N(\mathbf{r})}{\partial N^n}\right)_{v(\mathbf{r})} &= 0
      \text{ for } n \geqslant 3
    """

    def __init__(self, dict_density, n_max=None, global_softness=None):
        r"""Initialize quadratic density model to compute local reactivity descriptors.

        Parameters
        ----------
        dict_density : dict
            Dictionary of number of electrons (keys) and corresponding density array (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): \rho_{N_0 - 1}\left(\mathbf{
            r}\right), N_0: \rho_{N_0}\left(\mathbf{r}\right), (N_0 + 1): \rho_{N_0 + 1}\left(
            \mathbf{r}\right)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.
        """
        # check number of electrons & density values
        n_ref, dens_m, dens_0, dens_p = check_dict_values(dict_density)
        # compute fukui function & dual descriptor of N0-electron system
        self._ff0 = 0.5 * (dens_p - dens_m)
        self._df0 = dens_p - 2 * dens_0 + dens_m
        super(QuadraticLocalTool, self).__init__(n_ref, n_max, global_softness)
        self.dict_density = dict_density

    @doc_inherit(BaseLocalTool)
    def density(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute density
        rho = self.dict_density[self.n0].copy()
        rho += self._ff0 * (n_elec - self._n0) + 0.5 * self._df0 * (n_elec - self._n0)**2
        return rho

    @doc_inherit(BaseLocalTool)
    def density_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")
        if order == 1:
            deriv = self._ff0 + self._df0 * (n_elec - self.n0)
        elif order == 2:
            deriv = self._df0
        else:
            deriv = 0.
        return deriv


class QuadraticCondensedTool(BaseCondensedTool):
    r"""Condensed conceptual DFT reactivity descriptors class based on the quadratic energy model.

    This class contains the atom-condensed equivalent of :class:`QuadraticLocalTool` reactivity
    indicators.
    """

    def __init__(self, dict_population, n_max=None, global_softness=None):
        r"""Initialize quadratic population model to compute condensed reactivity descriptors.

        Parameters
        ----------
        dict_population : dict
            Dictionary of number of electrons (keys) and corresponding atomic populations array
            (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): {N_A \left(N_0 - 1\right)},
            N_0: {N_A \left(N_0\right)}, (N_0 + 1): {N_A \left(N_0 + 1\right)}`.
            The :math:`N_0` value is considered as the reference number of electrons.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.
        """
        # check number of electrons & density values
        n_ref, pop_m, pop_0, pop_p = check_dict_values(dict_population)
        # compute condensed fukui function & dual descriptor of N0-electron system
        self._ff0 = 0.5 * (pop_p - pop_m)
        self._df0 = pop_p - 2 * pop_0 + pop_m
        super(QuadraticCondensedTool, self).__init__(n_ref, n_max, global_softness)
        self.dict_population = dict_population

    @doc_inherit(BaseCondensedTool)
    def population(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # compute density
        pop = self.dict_population[self.n_ref].copy()
        pop += self._ff0 * (n_elec - self._n0) + 0.5 * self._df0 * (n_elec - self._n0)**2
        return pop

    @doc_inherit(BaseCondensedTool)
    def population_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")
        if order == 1:
            deriv = self._ff0 + self._df0 * (n_elec - self.n_ref)
        elif order == 2:
            deriv = self._df0
        else:
            deriv = 0.
        return deriv
