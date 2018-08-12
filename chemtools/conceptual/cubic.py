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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Cubic Energy Model.

This module contains the global and local tool classes corresponding to cubic energy models.
"""


import numpy as np

from horton import log

from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool


__all__ = ["CubicGlobalTool"]


class CubicGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the cubic energy model.

    The energy is approximated as a quadratic function of the number of electrons,

    .. math:: E(N) = a + b N + c N^2 + d N^3

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values, the unknown parameters
    of the energy model are obtained by interpolation.

    First, second and higher order derivatives of the cubic energy model with respect to
    the number of electrons at fixed external potential are given by:

    .. math::
       \left(\frac{\partial E}{\partial N}\right)_{v(\mathbf{r})} &= b + 2 c N  + 3 d N^2\\
       \left(\frac{\partial^2 E}{\partial N^2}\right)_{v(\mathbf{r})} &= 2 c  + 6 d N\\
       \left(\frac{\partial^3 E}{\partial N^3}\right)_{v(\mathbf{r})} &= 6 d \\
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} &= 0
             \quad \text{for} \quad n \geq 3
    """

    def __init__(self, dict_energy, omega=0.5):
        """Initialize a cubic energy model.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
        omega : float
            Value of omega parameter in the energy model.

        """
        # check energy values
        if len(dict_energy) != 3 or not all([key >= 0 for key in dict_energy.keys()]):
            raise ValueError('Cubic model requires 3 energy values corresponding '
                             'to positive number of electrons!')
        # find reference number of electrons
        n_ref = sorted(dict_energy.keys())[1]
        if n_ref < 1:
            raise ValueError('The n_ref cannot be less than one! Given n_ref={0}'.format(n_ref))
        if sorted(dict_energy.keys()) != [n_ref - 1, n_ref, n_ref + 1]:
            raise ValueError('Number of electrons should differ by one!')
        energy_m, energy_0, energy_p = [dict_energy[n] for n in sorted(dict_energy.keys())]

        self._omega = omega
        param_a = energy_0
        param_b = -omega * energy_m + 2. * omega * energy_0 - omega * energy_p
        param_b += energy_p - energy_0
        param_c = (energy_m - 2. * energy_0 + energy_p) / 2.
        param_d = (2. * omega - 1.) * (energy_m - 2. * energy_0 + energy_p) / 2.
        self._params = np.array([param_a, param_b, param_c, param_d])
        super(CubicGlobalTool, self).__init__(dict_energy, n_ref, None)

    @property
    def omega(self):
        r"""Parameter :math:`\omega` in the energy model."""
        return self._omega

    @property
    def params(self):
        r"""Parameters :math:`a`, :math:`b`, :math:`c`, and :math:`d` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        # compute the change in the number of electrons w.r.t. N0
        delta_n = n_elec - self._n0
        # compute energy
        result = self._params[0] + self._params[1] * delta_n + self._params[2] * delta_n ** 2
        result += self._params[3] * delta_n ** 3.
        return result

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy derivative evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # compute the change in the number of electrons w.r.t. N0
        delta_n = n_elec - self._n0
        # compute derivative of energy
        if order == 1:
            result = self._params[1] + 2. * self._params[2] * delta_n + \
                     3. * self._params[3] * delta_n ** 2.
        elif order == 2:
            result = 2. * self._params[2] + 6. * self._params[3] * delta_n
        elif order == 3:
            result = 6. * self._params[3]
        else:
            result = 0
        return result
