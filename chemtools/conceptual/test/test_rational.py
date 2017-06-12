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
"""Test chemtools.conceptual.rational Module."""

import math
import numpy as np
import sympy as sp
from chemtools.conceptual.rational import RationalGlobalTool


def test_global_rational_pnpp_energy():
    # E(N) = (0.5 - 2.2 N) / (1 + 0.7 N)
    ne = sp.symbols('ne')
    energy_model = (0.5 - 2.2 * ne) / (1 + 0.7 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    # Build rational global tool instance
    model = RationalGlobalTool(-1.6250, -1.96774193, -1.0, 2.0)
    # check parameters
    np.testing.assert_almost_equal(model.n0, 2.0, decimal=6)
    np.testing.assert_almost_equal(model.params[0], 0.5, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -2.2, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 0.7, decimal=6)
    # check energy values (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy(0), energy(0), decimal=6)
    np.testing.assert_almost_equal(model.energy(1), energy(1), decimal=6)
    np.testing.assert_almost_equal(model.energy(2), energy(2), decimal=6)
    np.testing.assert_almost_equal(model.energy(3), energy(3), decimal=6)
    np.testing.assert_almost_equal(model.energy(4), energy(4), decimal=6)
    np.testing.assert_almost_equal(model.energy(5), energy(5), decimal=6)
    np.testing.assert_almost_equal(model.energy(6), energy(6), decimal=6)
    np.testing.assert_almost_equal(model.energy(1.5), energy(1.5), decimal=6)
    np.testing.assert_almost_equal(model.energy(0.8), energy(0.8), decimal=6)
    # check energy derivatives (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy_derivative(0), deriv(0, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1), deriv(1, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(2), deriv(2, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3), deriv(3, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4), deriv(4, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(5), deriv(5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(6), deriv(6, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.5), deriv(1.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.8), deriv(0.8, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.5, 2), deriv(1.5, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.8, 2), deriv(0.8, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.1, 3), deriv(1.1, 3), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(2.5, 4), deriv(2.5, 4), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.65, 5), deriv(0.65, 5), decimal=5)
    np.testing.assert_almost_equal(model.energy_derivative(1.90, 6), deriv(1.90, 6), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3.20, 3), deriv(3.20, 3), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4.05, 7), deriv(4.05, 7), decimal=6)


def test_global_rational_pnpp_energy_reactivity():
    # E(N) = (0.5 - 2.2 N) / (1 + 0.7 N)
    ne = sp.symbols('ne')
    energy_model = (0.5 - 2.2 * ne) / (1 + 0.7 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    # Build rational global tool instance
    model = RationalGlobalTool(-1.6250, -1.96774193, -1.0, 2.0)
    # check global descriptors (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.ip, energy(1.0) - energy(2.0), decimal=6)
    np.testing.assert_almost_equal(model.ea, energy(2.0) - energy(3.0), decimal=6)
    np.testing.assert_almost_equal(model.mu, deriv(2.0, 1), decimal=6)
    np.testing.assert_almost_equal(model.eta, deriv(2.0, 2), decimal=6)
    np.testing.assert_almost_equal(model.ionization_potential, energy(1.0) - energy(2.0), decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, energy(2.0) - energy(3.0), decimal=6)
    np.testing.assert_almost_equal(model.chemical_potential, deriv(2.0, 1), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, deriv(2.0, 2), decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, -deriv(2.0, 1), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), deriv(2.0, 3), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), deriv(2.0, 4), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), deriv(2.0, 5), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(5), deriv(2.0, 6), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(6), deriv(2.0, 7), decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0 / deriv(2.0, 2), decimal=6)
    # check n_max and related descriptors (expected values are computed symbolically)
    np.testing.assert_equal(model.n_max, float('inf'))
    np.testing.assert_almost_equal(model.energy(model.n_max), -3.14285714, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, 1.51785714, decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, -1.17511520, decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, 2.14285714, decimal=6)


def test_global_rational_pnpp_grand_potential():
    # E(N) = (0.5 - 2.2 N) / (1 + 0.7 N)
    ne = sp.symbols('ne')
    energy_model = (0.5 - 2.2 * ne) / (1 + 0.7 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    grand = lambda n: energy(n) - n * deriv(n, 1)
    # Build rational global tool instance
    model = RationalGlobalTool(-1.6250, -1.96774193, -1.0, 2.0)
    # check grand potential (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential(1.), grand(1.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(2.), grand(2.0), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(3.), grand(3.0), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(2.78), grand(2.78), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(5.2), grand(5.2), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(0.), grand(0.), decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential(model.n_max), , decimal=6)
    # check grand potential derivative (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.), -2.0, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(1.4, 1), -1.4, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.86), -2.86, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2., 2), -3.87226890, decimal=6)
    # expected values based on derived formulas
    a0, a1, b1 = 0.5, -2.2, 0.7
    dE = lambda n, r: (-b1)**(r-1) * (a1 - a0*b1) * math.factorial(r) / (1 + b1*n)**(r+1)
    d2omega = lambda n: -1. / dE(n, 2)
    d3omega = lambda n: dE(n, 3) / dE(n, 2)**3
    # d4omega = lambda n: dE(n, 4) / dE(n, 2)**4 - (3 * dE(n, 3)**2) / dE(n, 2)**5
    # d5omega = lambda n: (dE(n, 5)/dE(n, 2)**5 - (10 * dE(n, 3) * dE(n, 4))/dE(n, 2)**6 +
    #                      (15 * dE(n, 3)**3)/dE(n, 2)**7)
    domega_n = model.grand_potential_derivative
    np.testing.assert_almost_equal(domega_n(4.67, 1), -4.67, decimal=6)
    np.testing.assert_almost_equal(domega_n(3.5, 2), d2omega(3.5), decimal=6)
    np.testing.assert_almost_equal(domega_n(4.1, 2), d2omega(4.1), decimal=6)
    np.testing.assert_almost_equal(domega_n(4.67, 2), d2omega(4.67), decimal=6)
    np.testing.assert_almost_equal(domega_n(2.9, 3), d3omega(2.9), decimal=5)
    np.testing.assert_almost_equal(domega_n(4.67, 3), d3omega(4.67), decimal=4)
    # np.testing.assert_almost_equal(domega_n(1.6, 4), d4omega(1.6), decimal=6)
    # np.testing.assert_almost_equal(domega_n(2.92, 4), d4omega(2.92), decimal=6)
    # np.testing.assert_almost_equal(domega_n(5.01, 5), d5omega(5.01), decimal=6)
    # np.testing.assert_almost_equal(domega_n(4.101, 5), d5omega(4.101), decimal=6)
    # np.testing.assert_almost_equal(domega_n(model.n_max, 4), decimal=6)
    # check mu to N conversion
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.4427083333), 2., decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.5799422391), 1.567, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.9515745573), 0.91, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.2641542934), 3.01, decimal=6)
    # check grand potential (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.125925925), -1.70370370, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.442708333), -0.73958333, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.232747054), -1.27423079, decimal=6)
    # check grand potential derivative (as a function of mu)
    domega_mu = model.grand_potential_mu_derivative
    np.testing.assert_almost_equal(domega_mu(dE(5.81, 1), 1), -5.81, decimal=6)
    np.testing.assert_almost_equal(domega_mu(dE(4.67, 1), 2), d2omega(4.67), decimal=5)
    # np.testing.assert_almost_equal(domega_mu(dE(6.45, 1), 3), d3omega(6.45), decimal=6)
    # np.testing.assert_almost_equal(domega_mu(dE(5.12, 1), 4), d4omega(5.12), decimal=6)


def test_global_rational_pnpp_grand_potential_reactivity():
    # E(N) = (0.5 - 2.2 N) / (1 + 0.7 N)
    n0, a0, a1, b1 = 2.0, 0.5, -2.2, 0.7
    # build global tool
    model = RationalGlobalTool(-1.6250, -1.96774193, -1.0, 2.0)
    # check hyper-softnesses
    expected = 3.0 * (1 + b1 * n0)**5 / (4 * b1 * (a1 - a0 * b1)**2)
    np.testing.assert_almost_equal(model.hyper_softness(2), expected, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.0, 3), -expected, decimal=6)
    expected = -15 * (1 + b1 * n0)**7 / (8 * b1 * (a1 - a0 * b1)**3)
    np.testing.assert_almost_equal(model.hyper_softness(3), expected, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.0, 4), -expected, decimal=6)


def test_global_rational_nnpp_energy():
    # E(N) = (-0.15 - 4.2 N) / (1 + 0.45 N)
    ne = sp.symbols('ne')
    energy_model = (-0.15 - 4.2 * ne) / (1.0 + 0.45 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    # build global tool
    model = RationalGlobalTool(-6.99363057, -7.23428571, -6.69064748, 6.5)
    # check parameters
    np.testing.assert_almost_equal(model.n0, 6.5, decimal=6)
    np.testing.assert_almost_equal(model.params[0], -0.15, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -4.2, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 0.45, decimal=6)
    # check energy values (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy(6.5), energy(6.5), decimal=6)
    np.testing.assert_almost_equal(model.energy(7.5), energy(7.5), decimal=6)
    np.testing.assert_almost_equal(model.energy(5.5), energy(5.5), decimal=6)
    np.testing.assert_almost_equal(model.energy(5.0), energy(5.0), decimal=6)
    np.testing.assert_almost_equal(model.energy(8.0), energy(8.0), decimal=6)
    # check energy derivatives (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy_derivative(6.5), deriv(6.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(7.5), deriv(7.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(5.5), deriv(5.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4, 2), deriv(4.0, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10., 3), deriv(10., 3), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(9.5, 4), deriv(9.5, 4), decimal=6)


def test_global_rational_nnpp_energy_reactivity():
    # E(N) = (-0.15 - 4.2 N) / (1 + 0.45 N)
    ne = sp.symbols('ne')
    energy_model = (-0.15 - 4.2 * ne) / (1.0 + 0.45 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    # build global tool
    model = RationalGlobalTool(-6.99363057, -7.23428571, -6.69064748, 6.5)
    # check global descriptors (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.ip, energy(5.5) - energy(6.5), decimal=6)
    np.testing.assert_almost_equal(model.ea, energy(6.5) - energy(7.5), decimal=6)
    np.testing.assert_almost_equal(model.mu, deriv(6.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.eta, deriv(6.5, 2), decimal=6)
    np.testing.assert_almost_equal(model.ionization_potential, energy(5.5) - energy(6.5), decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, energy(6.5) - energy(7.5), decimal=6)
    np.testing.assert_almost_equal(model.chemical_potential, deriv(6.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, deriv(6.5, 2), decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, -deriv(6.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), deriv(6.5, 3), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), deriv(6.5, 4), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), deriv(6.5, 5), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(5), deriv(6.5, 6), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(6), deriv(6.5, 7), decimal=6)
    np.testing.assert_almost_equal(model.softness, 1. / deriv(6.5, 2), decimal=6)
    # check n_max and related descriptors (expected values are computed symbolically)
    np.testing.assert_equal(model.n_max, float('inf'))
    np.testing.assert_almost_equal(model.energy(model.n_max), -9.33333333, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, 2.33970276, decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, -2.099047619, decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, 2.64268585, decimal=6)


def test_global_rational_nnpp_grand_potential():
    # E(N) = (-0.15 - 4.2 N) / (1 + 0.45 N)
    ne = sp.symbols('ne')
    energy_model = (0.5 - 2.2 * ne) / (1 + 0.7 * ne)
    energy = lambda n: energy_model.subs({'ne': n})
    deriv = lambda n, r: sp.diff(energy_model, 'ne', r).subs({'ne': n})
    model = RationalGlobalTool(-6.99363057, -7.23428571, -6.69064748, 6.5)
    # check grand potential (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential(6.5), -5.2500304, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(7.91), -5.7468530, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(0.), -0.15, decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential(model.n_max), , decimal=6)
    # check grand potential derivative (as a function of N)
    # expected values based on derived formulas
    a0, a1, b1 = -0.15, -4.2, 0.45
    dE = lambda n, r: (-b1)**(r-1) * (a1 - a0*b1) * math.factorial(r) / (1 + b1*n)**(r+1)
    d2omega = lambda n: -1. / dE(n, 2)
    d3omega = lambda n: dE(n, 3) / dE(n, 2)**3
    # d4omega = lambda n: dE(n, 4) / dE(n, 2)**4 - (3 * dE(n, 3)**2) / dE(n, 2)**5
    # d5omega = lambda n: (dE(n, 5)/dE(n, 2)**5 - (10 * dE(n, 3) * dE(n, 4))/dE(n, 2)**6 +
    #                      (15 * dE(n, 3)**3)/dE(n, 2)**7)
    domega_n = model.grand_potential_derivative
    np.testing.assert_almost_equal(domega_n(6.5, 1), -6.5, decimal=6)
    np.testing.assert_almost_equal(domega_n(7.1, 1), -7.1, decimal=6)
    np.testing.assert_almost_equal(domega_n(5.8, 2), d2omega(5.8), decimal=6)
    np.testing.assert_almost_equal(domega_n(0.0, 3), d3omega(0.0), decimal=6)
    # np.testing.assert_almost_equal(domega_n(8.01, 4), d4omega(8.01), decimal=6)
    # np.testing.assert_almost_equal(domega_n(6.901, 5), d5omega(6.901), decimal=6)
    # check mu to N conversion
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.2682461763), 6.5, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.2345757894), 7.105, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.1956803972), 7.99, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(-0.3568526811), 5.34, decimal=6)
    # check grand potential (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.26824617), -5.2500304, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.19153203), -5.8048876, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.20521256), -5.6965107, decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(model.n_max, 4), , decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.268246176), -5.2500304, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(-0.198782625), -5.7468530, decimal=6)
    # check grand potential derivative (as a function of mu)
    domega_mu = model.grand_potential_mu_derivative
    np.testing.assert_almost_equal(domega_mu(dE(6.301, 1), 1), -6.301, decimal=6)
    np.testing.assert_almost_equal(domega_mu(dE(5.55, 1), 2), d2omega(5.55), decimal=6)
    np.testing.assert_almost_equal(domega_mu(dE(6.99, 1), 3), d3omega(6.99), decimal=6)
    # np.testing.assert_almost_equal(domega_mu(dE(7.1, 1), 4), d4omega(7.1), decimal=6)
    # np.testing.assert_almost_equal(domega_mu(dE(7.6, 1), 5), d5omega(7.6), decimal=6)


def test_global_rational_nnpp_grand_potential_reactivity():
    # E(N) = (-0.15 - 4.2 N) / (1 + 0.45 N)
    n0, a0, a1, b1 = 6.5, -0.15, -4.2, 0.45
    # build global tool
    model = RationalGlobalTool(-6.99363057, -7.23428571, -6.69064748, 6.5)
    # check hyper-softnesses
    expected = 3.0 * (1 + b1 * n0)**5 / (4 * b1 * (a1 - a0 * b1)**2)
    np.testing.assert_almost_equal(model.hyper_softness(2), expected, decimal=5)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5, 3), -expected, decimal=5)
    expected = -15 * (1 + b1 * n0)**7 / (8 * b1 * (a1 - a0 * b1)**3)
    np.testing.assert_almost_equal(model.hyper_softness(3), expected, decimal=4)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5, 4), -expected, decimal=4)
