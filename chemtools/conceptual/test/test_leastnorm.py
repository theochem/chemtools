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
# pragma pylint: disable=protected-access
"""Test chemtools.conceptual.leastnorm Module."""

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import sympy as sp

from chemtools.conceptual.leastnorm import LeastNormGlobalTool


def make_symbolic_least_norm_model(omega, nth_order, weight=None, return_coeffs=False):
    """Return symbolic quadratic energy, energy derivative for unweighted and weighted model."""
    n_elec, j = sp.symbols('n_elec j')
    # Define Default Parameters used in all tests.
    n_ref = 10
    energy_zero0 = 100.
    energy_minus0 = 25.3
    energy_plus0 = 50.5
    ion = energy_minus0 - energy_zero0
    aff = energy_zero0 - energy_plus0
    dict_energy = {n_ref - 1: 25.3, n_ref: 100., n_ref + 1: 50.5}

    # Make Symbolic Parameters then Define the energy.
    param0 = energy_zero0
    param1 = - (omega * ion + (1 - omega) * aff)
    if weight is None:
        coefficient_even = 1.98619 * (ion - aff) / sp.factorial(j) ** 2
        coefficient_odd = 17.9551 * (2. * omega - 1) * (ion - aff) / sp.factorial(j) ** 2
    else:
        alpha_func = 1. / sp.factorial(n_elec) ** (2. - weight)
        beta_func = (-1.) ** n_elec / sp.factorial(n_elec) ** (2. - weight)

        alpha, beta = sp.symbols("alpha beta")
        weighted_coeff_even = 1. / (sp.factorial(j) ** (1. - weight) * (alpha + beta))
        weighted_coeff_odd = (2. * omega - 1) / (sp.factorial(j) ** (1. - weight) * (alpha + beta))
        coefficient_even = weighted_coeff_even * (ion - aff) / sp.factorial(j)
        coefficient_odd = weighted_coeff_odd * (ion - aff) / sp.factorial(j)

    params = [param0, param1]
    energy = param0 + param1 * n_elec
    for i in range(2, nth_order + 1):
        if i % 2 == 0.:
            coeff = coefficient_even.subs(j, i)
            params += [coeff]
        else:
            coeff = coefficient_odd.subs(j, i)
            params += [coefficient_odd.subs(j, i)]
        energy += coeff * n_elec**i

    # Define the first five derivatives.
    first_deriv = energy.diff(n_elec)
    sec_deriv = first_deriv.diff(n_elec)
    third_deriv = sec_deriv.diff(n_elec)
    fourth_deriv = third_deriv.diff(n_elec)
    fifth_deriv = fourth_deriv.diff(n_elec)
    expr = [energy, first_deriv, sec_deriv, third_deriv, fourth_deriv, fifth_deriv]
    if weight is not None:
        if return_coeffs:
            return n_ref, dict_energy, params, expr, alpha_func, beta_func, weighted_coeff_even,\
                weighted_coeff_odd
        return n_ref, dict_energy, params, expr, alpha_func, beta_func
    return n_ref, dict_energy, params, expr


def test_nth_order():
    # Test the number of terms for the unweighted least norm model.
    n_ref = 10
    dict_energy = {n_ref - 1: 25.3, n_ref: 100., n_ref + 1: 50.5}
    unweighted = LeastNormGlobalTool(dict_energy, 0.5, nth_order=5)
    assert_equal(len(unweighted.params), 6)

    weighted = LeastNormGlobalTool(dict_energy, 0.5, nth_order=10, weight=0.5)
    assert_equal(len(weighted.params), 11)


def test_parameters_unweighted():
    # Test parameters of the unweighted least norm.
    for omega in np.arange(-10., 10.):
        for order in range(0, 10):
            _, dict_energy, params, _ = make_symbolic_least_norm_model(omega, order)

            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
            assert_almost_equal(params, unweighted.params)


def test_parameters_weighted():
    # Test parameters of the weighted least norm.
    omega = 5.
    for n_degree in range(0, 5):
        # Test weighted model.
        for weight in np.linspace(0.01, 1., 5):
            n_elec = sp.symbols("n_elec")
            _, dict_energy, params, _, alpha, beta = make_symbolic_least_norm_model(omega, n_degree,
                                                                                    weight)
            alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
            beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
            weighted = LeastNormGlobalTool(dict_energy, omega, n_degree, weight, eps=1e-9)
            for i in range(2, n_degree):
                assert_almost_equal(params[i].subs([("alpha", alpha), ("beta", beta)]),
                                    weighted.params[i])


def test_energy_least_norm():
    # Test energy model of unweighted least norm.
    for omega in np.arange(-10., 10.):
        for order in range(0, 6):
            _, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            desired_values = [expr[0].subs([("n_elec", x)]) for x in range(1, 20)]
            actual_values = [unweighted.energy(x) for x in range(1, 20)]
            assert_almost_equal(actual_values, desired_values)


def test_energy_weighted_leastnorm():
    # Test energy model of weighted least norm.
    omega = -5.
    for order in range(0, 3):
        # Test weighted model.
        for weight in np.linspace(0.01, 1., 3):
            n_elec = sp.symbols("n_elec")
            _, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega, order,
                                                                                  weight)
            alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
            beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
            weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

            actual = [weighted.energy(n_elec) for n_elec in range(1, 20)]
            desired = [expr[0].subs([("n_elec", n_elec), ('alpha', alpha), ('beta', beta)]).evalf()
                       for n_elec in range(1, 20)]
            assert_almost_equal(actual, desired, decimal=4)


def test_derivative_energy():
    # Test derivative of energy model of unweighted least norm model.
    for omega in np.arange(-2., 2.):
        # Take first five derivatives so I need atleast five terms.
        for order in range(5, 8):
            _, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            # Go through each derivative up to fifth derivative.
            for d_order in range(1, 6):
                desired_values = [expr[d_order].subs([("n_elec", x)]) for x in range(1, 20)]
                actual_values = [unweighted.energy_derivative(x, d_order) for x in range(1, 20)]
                assert_almost_equal(actual_values, desired_values)


def test_derivative_energy_weighted():
    # Test derivative of energy model of weighted least norm model.
    omega = 2.
    for order in range(5, 7):
        # Test weighted model.
        weight = 0.5
        n_elec = sp.symbols("n_elec")
        _, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega, order,
                                                                              weight)
        alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
        beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
        weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

        # Go through each derivative up to fifth derivative.
        for d_order in range(1, 6):
            desired_values = [expr[d_order].subs([("n_elec", x), ("alpha", alpha),
                                                  ("beta", beta)]).evalf()
                              for x in range(1, 5)]
            actual_values = [weighted.energy_derivative(x, d_order) for x in range(1, 5)]
            assert_almost_equal(actual_values, desired_values, decimal=4)


def test_chemical_potential():
    # Test chemical potential of the least norm model.
    for omega in np.arange(-2., 2., step=0.5):
        for order in range(5, 10):
            n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            assert_almost_equal(unweighted.chemical_potential, expr[1].subs([("n_elec", n_ref)]))
            assert_almost_equal(unweighted.mu, expr[1].subs([("n_elec", n_ref)]))

            # Create weighted least norm.
            weight = 0.5
            n_elec = sp.symbols("n_elec")
            n_ref, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega,
                                                                                      order, weight)
            alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
            beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
            weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)
            desired = expr[1].subs([("n_elec", n_ref), ("alpha", alpha), ("beta", beta)]).evalf()
            assert_almost_equal(weighted.chemical_potential, desired)
            assert_almost_equal(weighted.mu, desired)


def test_chemical_hardness():
    # Test chemical hardness of least norm model.
    for omega in np.arange(-2., 2., step=0.5):
        for order in range(5, 10):
            n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            assert_almost_equal(unweighted.chemical_hardness, expr[2].subs([("n_elec", n_ref)]))
            assert_almost_equal(unweighted.eta, expr[2].subs([("n_elec", n_ref)]))

            # Test weight least norm model.
            weight = 0.5
            n_elec = sp.symbols("n_elec")
            n_ref, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega,
                                                                                      order, weight)
            alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
            beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
            weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

            desired = expr[2].subs([("n_elec", n_ref), ("alpha", alpha), ("beta", beta)]).evalf()
            assert_almost_equal(weighted.chemical_hardness, desired)
            assert_almost_equal(weighted.eta, desired)


def test_hyper_hardness():
    # Test hyperhardness for the least norm model.
    for omega in np.arange(-2., 2., step=0.5):
        for order in range(5, 10):
            n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            assert_almost_equal(unweighted.hyper_hardness(2), expr[3].subs([("n_elec", n_ref)]))
            assert_almost_equal(unweighted.hyper_hardness(3), expr[4].subs([("n_elec", n_ref)]))

            # Test weight least norm model.
            weight = 0.5
            n_elec = sp.symbols("n_elec")
            n_ref, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega,
                                                                                      order, weight)
            alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
            beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
            weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

            sub = [("n_elec", n_ref), ("alpha", alpha), ("beta", beta)]
            assert_almost_equal(weighted.hyper_hardness(2), expr[3].subs(sub).evalf())
            assert_almost_equal(weighted.hyper_hardness(3), expr[4].subs(sub).evalf())


def test_chemical_softness():
    # Test chemical softness for the least norm model.
    omega = -2.
    for order in range(5, 8):
        n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
        unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

        assert_almost_equal(unweighted.softness, 1. / expr[2].subs([("n_elec", n_ref)]))

        weight = 0.1
        n_elec = sp.symbols("n_elec")
        n_ref, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega,
                                                                                  order, weight)
        alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
        beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
        weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

        sub = [("n_elec", n_ref), ("alpha", alpha), ("beta", beta)]
        assert_almost_equal(weighted.softness, 1. / expr[2].subs(sub))


def test_hyper_softness():
    # Test hyper softness for the least norm model.
    omega = 5.
    for order in range(5, 8):
        n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
        unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
        substitute = [("n_elec", n_ref)]

        assert_almost_equal(unweighted.hyper_softness(2),
                            -expr[3].subs(substitute) / expr[2].subs(substitute) ** 3)

        weight = 0.2
        n_elec = sp.symbols("n_elec")
        n_ref, dict_energy, _, expr, alpha, beta = make_symbolic_least_norm_model(omega, order,
                                                                                  weight)
        alpha = sp.Sum(alpha, (n_elec, 2, 13)).doit()
        beta = sp.Sum(beta, (n_elec, 2, 13)).doit()
        weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

        sub = [("n_elec", n_ref), ("alpha", alpha), ("beta", beta)]
        desired = -expr[3].subs(sub) / expr[2].subs(sub)**3.
        assert_almost_equal(weighted.hyper_softness(2), desired)


def test_electron_aff_ionization():
    # Test electron affinity and ionization for the least norm model.
    for omega in np.arange(-5., 5.):
        for order in range(5, 10):
            n_ref, dict_energy, _, expr = make_symbolic_least_norm_model(omega, order)
            unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

            substitute = [("n_elec", n_ref)]
            assert_almost_equal(unweighted.electron_affinity,
                                expr[0].subs(substitute) - expr[0].subs([("n_elec", n_ref + 1.)]))
            assert_almost_equal(unweighted.ionization_potential,
                                expr[0].subs([("n_elec", n_ref - 1.)]) - expr[0].subs(substitute))


def test_alpha_and_beta():
    # Test the alpha and beta coefficients for the weighted least norm model.
    omega = 5.
    order = 7
    n_elec = sp.symbols("n_elec")

    for weight in np.arange(0.01, 1., 0.5):
        _, dict_energy, _, _, alpha, beta = make_symbolic_least_norm_model(omega, order, weight)
        weighted = LeastNormGlobalTool(dict_energy, omega, order, weight, eps=1e-9)

        # Test Beta
        actual = [weighted._beta_term(j) for j in range(1, 20)]
        desired = [beta.subs([(n_elec, j)]) for j in range(1, 20)]
        assert_almost_equal(actual, desired)

        # Test Alpha
        actual = [weighted._alpha_term(j) for j in range(1, 20)]
        desired = [alpha.subs([(n_elec, j)]) for j in range(1, 20)]
        assert_almost_equal(actual, desired)


def test_coefficients_weighted():
    # Test the coefficients for the weighted least norm model.
    omega = 5.
    order = 6
    n_elec = sp.symbols("n_elec")

    for weight in np.arange(0.01, 1., 0.4):
        _, energy, _, _, alpha, beta, even, odd = make_symbolic_least_norm_model(omega, order,
                                                                                 weight, True)
        weighted = LeastNormGlobalTool(energy, omega, order, weight, eps=1e-9)

        alpha = sp.Sum(alpha, (n_elec, 2, 15)).doit()
        beta = sp.Sum(beta, (n_elec, 2, 15)).doit()

        # Test Odd Weighted Coefficients
        actual = [weighted._coefficients_weighted(j, alpha, beta) for j in range(3, 10, 2)]
        desired = [odd.subs([("j", j), ("alpha", alpha), ("beta", beta)]).evalf()
                   for j in range(3, 10, 2)]
        assert_almost_equal(actual, desired)

        # Test Even Weighted Coefficients
        actual = [weighted._coefficients_weighted(j, alpha, beta) for j in range(2, 10, 2)]
        desired = [even.subs([("j", j), ("alpha", alpha), ("beta", beta)]).evalf()
                   for j in range(2, 10, 2)]
        assert_almost_equal(actual, desired)


def test_n_max():
    # Test getting the maximum for the unweighted least norm model.
    omega = 0.5
    order = 2
    _, dict_energy, _, _ = make_symbolic_least_norm_model(omega, order)
    unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)

    # Quadratic is unbounded below when a < 0. for a * x**2.
    for third_coefficient in np.arange(-10., -0.01):
        unweighted._params = [5., 2., third_coefficient]
        assert_equal(unweighted._compute_n_max(), np.inf)

    # Quadratic is bounded when a < 0.
    unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
    for third_coefficient in np.arange(0.01, 10.):
        unweighted._params = [5., 2., third_coefficient]
        actual = unweighted._compute_n_max()
        minima = -2. / (2. * third_coefficient)
        if minima < 0.:
            minima = 0.
        assert_equal(actual, minima)

    # Cubic, order=3
    order = 3
    _, dict_energy, _, _ = make_symbolic_least_norm_model(omega, order)
    unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
    unweighted._params = [0., 0., -3., 5.]
    assert_equal(unweighted._compute_n_max(), 0.4)

    # Quadric
    order = 4
    _, dict_energy, _, _ = make_symbolic_least_norm_model(omega, order)
    unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
    unweighted._params = [-3., 0, -3., 5., 5.]
    assert_almost_equal(unweighted._compute_n_max(), 0.2887, decimal=2)

    # Quantic
    order = 5
    _, dict_energy, _, _ = make_symbolic_least_norm_model(omega, order)
    unweighted = LeastNormGlobalTool(dict_energy=dict_energy, omega=omega, nth_order=order)
    unweighted._params = [-3., 0, -3., 5., -5., 10.]
    assert_almost_equal(unweighted._compute_n_max(), 0.4)

    # Leading coefficient is greater than zero, hence goes to infinity
    unweighted._params = [3., 0, 3., -5., 5., 110.]
    assert_almost_equal(unweighted._compute_n_max(), 0.)

    # The first coefficient is negative and thus energy goes to -infinity.
    unweighted._params = [0., 60., -32, -11, 8, -1.]
    assert_almost_equal(unweighted._compute_n_max(), np.inf, decimal=4)

    # Test getting global maxima when having multiple local maxima
    unweighted._params = [0., -60., 32, 11, -8, 1.]
    assert_almost_equal(unweighted._compute_n_max(), 4.35476, decimal=4)
