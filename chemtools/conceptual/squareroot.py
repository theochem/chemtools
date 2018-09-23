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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Square Root Energy Model.

This module contains the global tool class corresponding to square root energy models.
"""

import numpy as np

from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons

__all__ = ['SquareRootGlobalTool']


class SquareRootGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the square-root energy model.

    Given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy.

    The energy is approximated as a square-root function of the number of electrons,
    and the three unknown parameters are obtained by interpolating to the given values of energy.

    .. math:: E(N) = a + b \sqrt{N}  + c N

    First, second and higher order derivatives of the square-root energy model with respect to
    the number of electrons at fixed external potential are given by:

    .. math::
       \left(\frac{\partial E}{\partial N}\right)_{v(\mathbf{r})} &= \frac{b}{2\sqrt{N}} + c \\
       \left(\frac{\partial^2 E}{\partial N^2}\right)_{v(\mathbf{r})} &=
       \frac{-b}{4} N^{-\frac{3}{2}}\\
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} &=
       \frac{\Pi_{i=1}^{n-1} (2i - 1) a_1}
       {2^n (-1)^{n - 1} } N^{-(n-1)} \sqrt{N}
             \quad \text{for} \quad n \geq 2
    """
    def __init__(self, dict_energy):
        # check number of electrons & energy values
        n_ref, energy_m, energy_0, energy_p = check_dict_values(dict_energy)

        # Compute The Coefficients
        a1 = (energy_p - 2. * energy_0 + energy_m) / \
             (np.sqrt(n_ref + 1) - 2. * np.sqrt(n_ref) + np.sqrt(n_ref - 1))
        a2 = (energy_p - energy_m) / 2. - a1 * (np.sqrt(n_ref + 1) - np.sqrt(n_ref - 1)) / 2.
        a0 = energy_0 - a1 * np.sqrt(n_ref) - a2 * n_ref
        self._params = [a0, a1, a2]
        n_max = self._compute_nmax()
        super(SquareRootGlobalTool, self).__init__(n_ref, n_max)

    @property
    def params(self):
        """
        Parameters :math:`a`, :math:`b` and :math:`c` of energy model.
        """
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)

        # Square Root Model goes to infinity as N goes to infinity.
        if np.isinf(n_elec):
            return np.inf

        return self._params[0] + self._params[1] * np.sqrt(n_elec) + self._params[2] * n_elec

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self._n0 - 1, self._n0 + 1)

        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")

        # Evaluate Derivative
        if order == 1:
            if n_elec == np.inf:
                # The limit as N goes to infinity on the first order derivative is a2
                return self._params[2]
            deriv_value = self._params[2] + self._params[1] / (2. * np.sqrt(n_elec))
        else:
            if n_elec == np.inf:
                # Limit as N goes to infinity on the higher order derivative is zero
                return 0
            coefficient_factor = np.prod(2. * np.arange(1, order) - 1) * self._params[1]
            coefficient_factor /= (2.**order * (-1)**(order - 1))
            deriv_value = coefficient_factor * n_elec**(-(order - 1)) * np.sqrt(n_elec**(-1))
        return deriv_value

    def _compute_nmax(self):
        # Compute the local minimum, n_max
        _, a1, a2 = self.params

        if a1 == 0.:
            # If a1 is zero then it is just a linear model.
            raise ValueError("Coefficient a1 cannot be zero or else it is a linear model.")
        # a1 requires to be negative to ensure minimum is concave up
        if a1 < 0.:
            if a2 <= 0.:
                # Energy goes to negative infinity as it is monotonically decreasing.
                n_max = np.inf
            elif a2 > 0.:
                # Local Minima exists.
                n_max = (a1 / (2. * a2)) ** 2
        elif a1 > 0.:
            if a2 >= 0.:
                # Minima occurs at the origin as it is monotonically increasing.
                n_max = 0.
            elif a2 < 0.:
                # Energy goes to -infinity.
                n_max = np.inf
        return n_max
