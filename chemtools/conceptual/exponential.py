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
from horton import log
from chemtools.conceptual.base import BaseGlobalTool
from chemtools.utils.utils import doc_inherit

__all__ = ['ExponentialGlobalTool']


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

    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # check energy values are monotonic, i.e. E(N-1) > E(N) > E(N+1)
        if not energy_minus > energy_zero > energy_plus:
            energies = [energy_minus, energy_zero, energy_plus]
            n_values = [n0 - 1, n0, n0 + 1]
            raise ValueError('To interpolate exponential energy model, E values vs. N should be ' +
                             'monotonic! Given E={0} for N={1}.'.format(energies, n_values))

        # calculate the A, B, gamma parameters of the model and N_max
        param_a = (energy_minus - energy_zero) * (energy_zero - energy_plus)
        param_a /= (energy_minus - 2 * energy_zero + energy_plus)
        param_b = energy_zero - param_a
        gamma = (energy_minus - 2 * energy_zero + energy_plus) / (energy_plus - energy_zero)
        gamma = math.log(1. - gamma)
        self._params = [param_a, gamma, param_b]
        # calulate N_max
        n_max = float('inf')
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        r"""Parameter :math:`A`, :math:`\gamma` and :math:`B` of energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
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
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy derivative evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if np.isinf(n_elec):
            # limit of E(N) derivatives as N goes to infinity equals zero
            deriv = 0.0
        else:
            dn = n_elec - self._n0
            deriv = self._params[0] * (- self._params[1])**order * math.exp(- self._params[1] * dn)
        return deriv
