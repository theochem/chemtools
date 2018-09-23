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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Least Norm.

This module contains the global tool class corresponding to least norm energy models.
"""

import numpy as np
from scipy.misc import factorial

from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool
from chemtools.conceptual.utils import check_dict_values, check_number_electrons


__all__ = ['LeastNormGlobalTool']


class LeastNormGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the least norm energy model.

    It is motivated by the desire for a model based on the n-th degree polynomial that is smooth
    and normalized.

    .. math::
        Taylor Series for Energy Model
        E(N) = E(N_0) - (\omega I + (1 - \omega)A) N + \sum_{j=2}^{\inf}\frac{a_j(I - A}{j!} N^j

        given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy and
        given \omega is a real number.
        The coefficient :math: 'a_j' is given by: \\
            a_j = \[
                    \begin{cases}
                    \frac{1.98619}{j!} & j is even\\
                    \frac{17.9551(2. \omega 0 1)}{j!} & j is odd
                    \end{cases}
                   \]

            a_j = \frac{1}{(j!)^(1-w)}\[
                    \begin{cases}
                    \frac{1.}{\alpha_x + \beta_X} & j is even\\
                    \frac{(2. \omega 0 1)}{\alpha_x + \beta_x} & j is odd
                    \end{cases}
                \], for the unweighted and weighted coefficients respectively.

        The Computation for the hyper-parameters for the weighted coefficients are given as:
            \alpha_x = \sum_{n=2}^{\inf}\frac{1}{(n!)^(2 - x)}
            \beta_x = \sum_{n=2}^{\inf}\frac{(-1)^n}{(n!)^(2 - x)}

    Parameters are computed by the formulas above and approximated in the weighted case by the
    the use of the epsilon parameter.

    For more information on this model, please see the reference.

    References
    ----------
    .. [1] When is the Fukui Function Not Normalized? The Danger of Inconsistent Energy
        Interpolation Models in Density Functional Theory.
        Farnaz Heidar-Zadeh, Ramón Alain Miranda-Quintana, Toon Verstraelen, Patrick Bultinck,
        and Paul W. Ayers. Journal of Chemical Theory and Computation 2016 12 (12), 5777-5787 DOI:
        10.1021/acs.jctc.6b00494
    """

    def __init__(self, dict_energy, omega, nth_order=5, weight=0., eps=1e-3):
        r"""Initialize the least-norm energy model to compute global reactivity descriptors.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): E(N_0 - 1), N_0: E(N_0),
            (N_0 + 1): E(N_0 + 1)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        omega : float, default=0.5
            The parameter :math:`\omega` in the energy model. Reasonable choices are between 0 and 1
            see reference for more details.
        nth_order : int, default=5
            The nth degree for the taylor polynomial.
        weight : float, default=0.
            The weight between 0 and 1 for the weighted least norm energy model. The weight with
            value zero gives the un-weighted least norm model.
        eps : float, default=1e-3
            The accuracy for computing the hyper-parameters :math: '\alpha_x' and :math: '\beta_x',
            only needed for the weighted model.
        """
        if weight > 1. or weight < 0.:
            raise ValueError("Weights have to be between 0 and 1. It is {0}".format(weight))
        # check number of electrons & energy values
        n_ref, energy_m, energy_0, energy_p = check_dict_values(dict_energy)

        self._omega = omega
        self._weight = weight
        self._nth_order = nth_order
        # Default coefficients for the first two terms.
        self.dict_energy = dict_energy

        # Need these to initially define the parameters
        ionization = energy_m - energy_0
        affinity = energy_0 - energy_p
        diff = ionization - affinity
        self._params = [energy_0, -(omega * ionization + (1 - omega) * affinity)]

        # Compute the coefficients up-to nth order, specified by `nth_order` argument.
        if weight != 0.:
            # If weighted model, compute the weighted coefficients.
            alpha, beta = self._hyperparameters(eps)
            self._params += [self._coefficients_weighted(j, alpha, beta) * diff / factorial(j)
                             for j in range(2, self._nth_order + 1)]
        else:
            self._params += [self._coefficients_unweighted(j) * diff / factorial(j)
                             for j in range(2, self._nth_order + 1)]

        self._n_max = self._compute_n_max()
        super(LeastNormGlobalTool, self).__init__(n_ref, self._n_max)

    @property
    def omega(self):
        r"""Constant for computing the coefficients."""
        return self._omega

    @property
    def weight(self):
        r"""Weight for computing coefficients for the weighted least-norm model."""
        return self._weight

    @property
    def nth_order(self):
        r"""Degree of the highest polynomial, truncated from the taylor polynomial."""
        return self._nth_order

    @property
    def params(self):
        r"""Parameters :math:`a_j` of the least-norm energy model."""
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        energy_val = 0.
        for i, coefficient in enumerate(self.params):
            energy_val += coefficient * n_elec ** i
        return energy_val

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        # check n_elec argument
        check_number_electrons(n_elec, self.n0 - 1., self.n0 + 1.)

        # check order
        if not (isinstance(order, int) and order > 0):
            raise ValueError("Argument order should be an integer greater than or equal to 1.")

        deriv = 0
        # Evaluate the derivative of each term of the energy model, evaluated at n_elec.
        for term in range(order, self._nth_order + 1):
            diff = term - order
            deriv += n_elec ** diff * self._params[term] * factorial(term) / factorial(diff)
        return deriv

    def _compute_n_max(self):
        # Leading coefficient dictates whether it is bounded or not.
        if self.params[-1] < 0.:
            return np.inf

        # Solve for the roots of the derivative
        deriv = np.flip(self.params[1:] * np.arange(1, self._nth_order + 1), 0)
        roots = np.roots(deriv)

        # Use the second derivative test.
        sec_deriv_params = np.flip(deriv[::-1][1:] * np.arange(2, self._nth_order + 1), 0)
        roots_output = []
        for i, root in enumerate(roots):
            sec_deriv = np.polyval(sec_deriv_params, root)
            if sec_deriv > 0. and np.isreal(root):
                roots_output.append(roots[i])
        roots_output = np.array(roots_output)

        if len(roots_output) == 0.:
            # No Maximum were found.
            return None

        # All of the roots were negative, then 0 electron has smallest value.
        if np.all(np.array(roots_output) < 0.):
            return 0.

        # Compare Energy Levels of different roots and get the one with maximum value.
        energy_level = np.array([self.energy(x) for x in roots_output[roots_output >= 0.]])
        return roots_output[roots_output >= 0.][energy_level.argmin()]

    def _coefficients_unweighted(self, j):
        # Get Coefficients/Parameters of the un-weighted energy model.
        if j % 2 == 0:
            return 1.98619 / factorial(j)
        return 17.9551 * (2. * self.omega - 1) / factorial(j)

    def _coefficients_weighted(self, j, alpha, beta):
        # Get Coefficients/Parameters of the weighted energy model.
        assert self._weight != 0.
        factor = 1. / (factorial(j) ** (1. - self._weight) * (alpha + beta))
        if j % 2 == 0:
            return factor
        return (2. * self._omega - 1.) * factor

    def _alpha_term(self, j):
        # Only needed for the weighted model to calculate the coefficients.
        return 1. / factorial(j) ** (2. - self._weight)

    def _beta_term(self, j):
        # Only needed for the weighted model to calculate the coefficients.
        return (-1.) ** j / factorial(j) ** (2. - self._weight)

    def _hyperparameters(self, eps):
        r"""Hyper-parameters needed to define the coefficients for the weighted model."""
        # The series, to calculate hyper-parameters (alpha, beta), start at term=2.
        prev_alpha, alpha = 1e10, self._alpha_term(2)
        prev_beta, beta = 1e10, self._beta_term(2)

        term = 3
        while np.abs(alpha - prev_alpha) > eps or np.abs(beta - prev_beta) > eps:
            prev_alpha = alpha.copy()
            prev_beta = beta.copy()
            alpha += self._alpha_term(term)
            beta += self._beta_term(term)
            term += 1
        return alpha, beta
