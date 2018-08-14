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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Exponential Energy Model.

This module contains the global tool class corresponding to exponential energy models.
"""


import math
import numpy as np

from chemtools.conceptual.base import BaseGlobalTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons
from chemtools.utils.utils import doc_inherit


__all__ = ["ExponentialGlobalTool"]


class ExponentialGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the exponential energy model.

    The energy is approximated as a exponential function of the number of electrons,

    .. math:: E(N) = A \exp(-\gamma(N-N_0)) + B

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values, the unknown parameters
    of the energy model are obtained by interpolation.

    The :math:`n^{\text{th}}`-order derivative of the rational energy model with respect to
    the number of electrons at fixed external potential is given by:

    .. math::
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} =
              A (-\gamma)^n \exp(-\gamma (N - N_0))
    """

    def __init__(self, dict_energy):
        r"""Initialize exponential energy model to compute global reactivity descriptors.

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
        if not energy_m > energy_0 > energy_p:
            energies = [energy_m, energy_0, energy_p]
            raise ValueError("For exponential model, the energy values for consecutive number of "
                             "electrons should be monotonic! E={0}".format(energies))
        # calculate parameters A, B, gamma parameters of the exponential model
        param_a = (energy_m - energy_0) * (energy_0 - energy_p)
        param_a /= (energy_m - 2 * energy_0 + energy_p)
        param_b = energy_0 - param_a
        param_g = (energy_m - 2 * energy_0 + energy_p) / (energy_p - energy_0)
        param_g = math.log(1. - param_g)
        self._params = [param_a, param_g, param_b]
        # calculate N_max
        n_max = float("inf")
        super(ExponentialGlobalTool, self).__init__(n_ref, n_max)
        self.dict_energy = dict_energy

    @property
    def params(self):
        r"""Parameter :math:`A`, :math:`\gamma` and :math:`B` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)
        # evaluate energy
        if np.isinf(n_elec):
            # limit of E(N) as N goes to infinity equals B
            value = self._params[2]
        else:
            dn = n_elec - self._n0
            value = self._params[0] * math.exp(- self._params[1] * dn) + self._params[2]
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
            dn = n_elec - self._n0
            deriv = self._params[0] * (- self._params[1])**order * math.exp(- self._params[1] * dn)
        return deriv
