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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Rational Energy Model.

This module contains the global tool class corresponding to rational energy models.
"""


import math
import numpy as np

from chemtools.conceptual.base import BaseGlobalTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons
from chemtools.utils.utils import doc_inherit


__all__ = ['RationalGlobalTool']


class RationalGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the rational energy model.

    The energy is approximated as a 3-parameter rational function of the number of electrons,

    .. math:: E(N) = \frac{a_0 + a_1 N}{1 + b_1 N}

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values, the unknown parameters
    of the energy model are obtained by interpolation.

    The :math:`n^{\text{th}}`-order derivatives of the rational energy model with respect to
    the number of electrons at fixed external potential is given by:

    .. math::
       \left(\frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})} =
             \frac{b_1^{n - 1} (a_1 - a_0 b_1) n!}{(1 + b_1 N)^{2n}}
    """

    def __init__(self, dict_energy):
        r"""Initialize rational energy model to compute global reactivity descriptors.

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
        # check energy values
        if not energy_m > energy_0 >= energy_p:
            energies = [energy_m, energy_0, energy_p]
            raise ValueError("For rational model, the energy values for consecutive number of "
                             "electrons should be monotonic! E={0}".format(energies))
        # calculate parameters a0, a1 and b1 of rational energy model
        param_b1 = - (energy_p - 2 * energy_0 + energy_m)
        param_b1 /= ((n_ref + 1) * energy_p - 2 * n_ref * energy_0 + (n_ref - 1) * energy_m)
        param_a1 = (1 + param_b1 * n_ref) * (energy_p - energy_0) + (param_b1 * energy_p)
        param_a0 = - param_a1 * n_ref + energy_0 * (1 + param_b1 * n_ref)
        self._params = [param_a0, param_a1, param_b1]
        # calculate N_max
        n_max = float('inf')
        super(RationalGlobalTool, self).__init__(n_ref, n_max)
        self.dict_energy = dict_energy

    @property
    def params(self):
        """Parameter :math:`a_0`, :math:`a_1` and :math:`b_1` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # evaluate energy
        if np.isinf(n_elec):
            # limit of E(N) as N goes to infinity equals a1/b1
            value = self._params[1] / self._params[2]
        else:
            value = (self._params[0] + self._params[1] * n_elec) / (1 + self._params[2] * n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")
        # evaluate derivative
        if np.isinf(n_elec):
            # limit of E(N) derivatives as N goes to infinity equals zero
            deriv = 0.0
        else:
            deriv = (-self._params[2])**(order - 1)
            deriv *= (self._params[1] - self._params[0] * self._params[2]) * math.factorial(order)
            deriv /= (1 + self._params[2] * n_elec)**(order + 1)
        return deriv
