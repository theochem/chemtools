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
"""Test chemtools.conceptual.exponential Module."""

import sympy as sp
from numpy.testing import assert_raises, assert_equal, assert_almost_equal
from chemtools.conceptual.exponential import ExponentialGlobalTool


def make_symbolic_exponential_model(a, b, c, n0):
    """Return symbolic quadratic energy, energy derivative & grand potential expressions."""
    ne = sp.symbols('ne')
    energy = a * sp.exp(b * (ne - n0)) + c
    expr = (lambda n: float(energy.subs({'ne': n})),
            lambda n, r: float(sp.diff(energy, 'ne', r).subs({'ne': n})),
            lambda n: float(energy.subs({'ne': n}) - n * sp.diff(energy, 'ne', 1).subs({'ne': n})))
    return expr


def make_analytical_grand_derivatives(deriv):
    """Return analytical 2nd, 3rd, 4th and 5th derivatives of grand potential w.r.t. mu."""
    expr = (lambda n: -1.0/deriv(n, 2),
            lambda n: deriv(n, 3)/deriv(n, 2)**3,
            lambda n: deriv(n, 4)/deriv(n, 2)**4 - (3*deriv(n, 3)**2)/deriv(n, 2)**5,
            lambda n: (deriv(n, 5)/deriv(n, 2)**5 - 10*deriv(n, 3)*deriv(n, 4)/deriv(n, 2)**6 +
                       15*deriv(n, 3)**3/deriv(n, 2)**7))
    return expr


def test_global_exponential_raises():
    # check invalid energy values
    assert_raises(ValueError, ExponentialGlobalTool, {5.: 15.0, 6.: 16.5, 4.: 18.1})
    assert_raises(ValueError, ExponentialGlobalTool, {6.: -15.0, 7.: -16.5, 5.: -18.1})
    assert_raises(ValueError, ExponentialGlobalTool, {10: -15.0, 11: -14.5, 9: -16.0})
    assert_raises(ValueError, ExponentialGlobalTool, {8: -15.0, 9: -14.9, 7: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {8: 15.0, 9: 14.9, 7: 14.0, 10: 15.0})
    # check invalid N0
    assert_raises(ValueError, ExponentialGlobalTool, {0: -15.0, 1: -14.4, -1: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {0.3: -15.0, 1.3: -14.4, -0.7: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {0.98: -15.0, 1.98: -14.4, -0.02: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {-1.: -15.0, 0.: -14.9, -2.: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {-2: -15.0, -1: -14.9, -3: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {-2: -15.0, -1: -14.9, -3: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {0.0: -15.0, 0.5: -14.9, 1.0: -14.0})
    assert_raises(ValueError, ExponentialGlobalTool, {2.: -15.0, 3.5: -14.9, 4.0: -14.0})
    # check invalid N
    model = ExponentialGlobalTool({5.: 5.2, 6.: 4.8, 4.: 6.0})
    assert_raises(ValueError, model.energy, -0.005)
    assert_raises(ValueError, model.energy, -1.35)
    assert_raises(ValueError, model.energy, -2.45)
    assert_raises(ValueError, model.energy_derivative, -0.05, 1)
    assert_raises(ValueError, model.energy_derivative, -1.91, 2)
    # check invalid derivative order
    assert_raises(ValueError, model.energy_derivative, 5.0, 1.)
    assert_raises(ValueError, model.energy_derivative, 5.0, 0.2)
    assert_raises(ValueError, model.energy_derivative, 5.0, -1)
    assert_raises(ValueError, model.energy_derivative, 5.0, -3)
    assert_raises(ValueError, model.energy_derivative, 5, '1')
    assert_raises(ValueError, model.energy_derivative, 5, [1])
    assert_raises(ValueError, model.energy_derivative, 3, 1.1)


def test_global_exponential_energy():
    # E(N) = 5.0 * exp(-0.1 * (N - 10)) + 3.0
    energy, deriv, _ = make_symbolic_exponential_model(5.0, -0.1, 3.0, 10.)
    # build exponential global tool
    model = ExponentialGlobalTool({10: 8.0, 11: 7.524187090179797, 9: 8.525854590378238})
    assert_almost_equal(model.params[0], 5.0, decimal=6)
    assert_almost_equal(model.params[2], 3.0, decimal=6)
    assert_almost_equal(model.params[1], 0.1, decimal=6)
    assert_almost_equal(model.n0, 10, decimal=6)
    # check E(N)
    assert_almost_equal(model.energy(20), energy(20), decimal=6)
    assert_almost_equal(model.energy(10), energy(10), decimal=6)
    assert_almost_equal(model.energy(8), energy(8), decimal=6)
    assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    assert_almost_equal(model.energy_derivative(18.1), deriv(18.1, 1), decimal=6)
    assert_almost_equal(model.energy_derivative(10), deriv(10, 1), decimal=6)
    assert_almost_equal(model.energy_derivative(8.5), deriv(8.5, 1), decimal=6)
    assert_almost_equal(model.energy_derivative(12.25), deriv(12.25, 1), decimal=6)
    # check d2E(N)
    assert_almost_equal(model.energy_derivative(17.3, 2), deriv(17.3, 2), decimal=6)
    assert_almost_equal(model.energy_derivative(10, 2), deriv(10, 2), decimal=6)
    assert_almost_equal(model.energy_derivative(9.1, 2), deriv(9.1, 2), decimal=6)
    # check d^nE(N) for n > 2
    assert_almost_equal(model.energy_derivative(13.7, 3), deriv(13.7, 3), decimal=6)
    assert_almost_equal(model.energy_derivative(10, 4), deriv(10, 4), decimal=6)
    assert_almost_equal(model.energy_derivative(4.5, 5), deriv(4.5, 5), decimal=6)
    assert_almost_equal(model.energy_derivative(6.5, 10), deriv(6.5, 10), decimal=6)


def test_global_exponential_energy_reactivity():
    # E(N) = 5.0 * exp(-0.1 * (N - 10)) + 3.0
    energy, deriv, _ = make_symbolic_exponential_model(5.0, -0.1, 3.0, 10)
    ip = energy(9) - energy(10)
    ea = energy(10) - energy(11)
    # build exponential global tool
    model = ExponentialGlobalTool({10: 8.0, 11: 7.524187090179797, 9: 8.525854590378238})
    # check ionization potential and electron affinity
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electronegativity, -deriv(10, 1), decimal=6)
    # check chemical potential, chemical hardness, and related tools
    assert_almost_equal(model.chemical_potential, deriv(10, 1), decimal=6)
    assert_almost_equal(model.chemical_hardness, deriv(10, 2), decimal=6)
    assert_almost_equal(model.mu, deriv(10, 1), decimal=6)
    assert_almost_equal(model.eta, deriv(10, 2), decimal=6)
    assert_almost_equal(model.hyper_hardness(2), deriv(10, 3), decimal=6)
    assert_almost_equal(model.hyper_hardness(3), deriv(10, 4), decimal=6)
    assert_almost_equal(model.hyper_hardness(4), deriv(10, 5), decimal=6)
    assert_almost_equal(model.softness, 1.0 / deriv(10, 2), decimal=6)
    assert_almost_equal(model.softness, 1.0 / (5. * 0.1**2), decimal=6)
    # check n_max and related descriptors (expected values are computed symbolically)
    assert_equal(model.n_max, float('inf'))
    assert_almost_equal(model.energy(model.n_max), 3.0, decimal=6)
    assert_almost_equal(model.energy_derivative(model.n_max), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(model.n_max, 2), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(model.n_max, 3), 0.0, decimal=6)
    assert_almost_equal(model.electrophilicity, 5.0, decimal=6)
    assert_almost_equal(model.nucleofugality, -4.52418709, decimal=6)
    assert_almost_equal(model.electrofugality, 5.52585459, decimal=6)
    # assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)


def test_global_exponential_grand_potential_n():
    # E(N) = 5.0 * exp(-0.1 * (N - 10)) + 3.0
    _, deriv, grand = make_symbolic_exponential_model(5.0, -0.1, 3.0, 10)
    # build exponential global tool
    model = ExponentialGlobalTool({10: 8.0, 11: 7.524187090179797, 9: 8.525854590378238})
    # check grand potential (as a function of N)
    assert_almost_equal(model.grand_potential(9.), grand(9.), decimal=6)
    assert_almost_equal(model.grand_potential(10), grand(10), decimal=6)
    assert_almost_equal(model.grand_potential(11.), grand(11.), decimal=6)
    assert_almost_equal(model.grand_potential(9.123), grand(9.123), decimal=6)
    assert_almost_equal(model.grand_potential(8.5), grand(8.5), decimal=6)
    assert_almost_equal(model.grand_potential(10.001), grand(10.001), decimal=6)
    # check grand potential derivative (as a function of N)
    assert_almost_equal(model.grand_potential_derivative(9., 1), -9, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.5, 1), -10.5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(11.001, 1), -11.001, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(9.15, 1), -9.15, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(model.n0, 1), -10, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10., 2), -1/(5. * 0.1**2), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10., 3), -1/(5.**2 * 0.1**3), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10., 4), -2/(5.**3 * 0.1**4), decimal=6)
    # get analytical expression of 2nd, 3rd, 4th and 5th derivatives
    d2oemga, d3omega, d4omega, d5omega = make_analytical_grand_derivatives(deriv)
    assert_almost_equal(model.grand_potential_derivative(9.6, 2), d2oemga(9.6), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.76, 2), d2oemga(10.76), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(9.23, 3), d3omega(9.23), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(11.56, 3), d3omega(11.56), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(9.23, 3), d3omega(9.23), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(8.67, 4), d4omega(8.67), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.4, 4), d4omega(10.4), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(12.6, 4), d4omega(12.6), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(11.5, 5), d5omega(11.5), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(13., 5), d5omega(13.), decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.02, 5), d5omega(10.02), decimal=6)


def test_global_exponential_grand_potential_mu():
    # E(N) = 5.0 * exp(-0.1 * (N - 10)) + 3.0
    _, deriv, grand = make_symbolic_exponential_model(5.0, -0.1, 3.0, 10)
    # build exponential global tool
    model = ExponentialGlobalTool({10: 8.0, 11: 7.524187090179797, 9: 8.525854590378238})
    # check mu to N conversion
    assert_almost_equal(model.convert_mu_to_n(-0.5), 10., decimal=6)
    assert_almost_equal(model.convert_mu_to_n(-0.5476345026), 9.09, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(-0.4230574978), 11.671, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(-0.6086285202), 8.034, decimal=6)
    # check grand potential derivative (as a function of mu)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81, 1), 1), -5.81, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(10.5, 1), 1), -10.5, decimal=6)
    # get analytical expression of 2nd, 3rd, 4th and 5th derivatives
    d2oemga, d3omega, d4omega, d5omega = make_analytical_grand_derivatives(deriv)
    mu = deriv(11.2, 1)
    assert_almost_equal(model.grand_potential_mu_derivative(mu, 2), d2oemga(11.2), decimal=6)
    mu = deriv(11.2, 1)
    assert_almost_equal(model.grand_potential_mu_derivative(mu, 3), d3omega(11.2), decimal=6)
    mu = deriv(9.34, 1)
    assert_almost_equal(model.grand_potential_mu_derivative(mu, 4), d4omega(9.34), decimal=6)
    mu = deriv(10.7, 1)
    assert_almost_equal(model.grand_potential_mu_derivative(mu, 5), d5omega(10.7), decimal=6)
    # check grand potential (as a function of mu)
    assert_almost_equal(model.grand_potential_mu(deriv(9., 1)), grand(9.), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(10.1, 1)), grand(10.1), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(11., 1)), grand(11.), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(12.3, 1)), grand(12.3), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(20.4, 1)), grand(20.4), decimal=6)


def test_global_exponential_grand_potential_reactivity():
    # E(N) = 5.0 * exp(-0.1 * (N - 10)) + 3.0
    # build exponential global tool
    model = ExponentialGlobalTool({10: 8.0, 11: 7.524187090179797, 9: 8.525854590378238})
    # check hyper-softnesses
    assert_almost_equal(model.hyper_softness(2), 1.0 / (5.**2 * 0.1**3), decimal=6)
    assert_almost_equal(model.hyper_softness(3), 2.0 / (5.**3 * 0.1**4), decimal=6)
