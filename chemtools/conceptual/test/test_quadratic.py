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
# pragma pylint: disable=invalid-name
"""Test chemtools.conceptual.quadratic Module."""

import numpy as np
from numpy.testing import assert_raises, assert_equal, assert_almost_equal
from chemtools.conceptual.quadratic import QuadraticGlobalTool, QuadraticLocalTool
from chemtools.conceptual.quadratic import QuadraticCondensedTool


def make_symbolic_quadratic_model(a, b, c):
    """Return quadratic energy, energy derivative & grand potential expressions."""
    symbolic_expr = (lambda n: a * n * n + b * n + c,
                     lambda n: 2. * a * n + b,
                     lambda n: a * n * n + b * n + c - n * (2. * a * n + b))
    return symbolic_expr


def test_global_quadratic_raises():
    # check invalid N0
    assert_raises(ValueError, QuadraticGlobalTool, {0: -5.5, 1: -6.0, -1: -7.0})
    assert_raises(ValueError, QuadraticGlobalTool, {0.1: -5.5, 1.1: -6.0, -0.9: -7.0})
    assert_raises(ValueError, QuadraticGlobalTool, {0.99: -5.5, 1.99: -6.0, -0.1: -7.0})
    assert_raises(ValueError, QuadraticGlobalTool, {-1.: -5.5, 0.: -6.0, -2.: -7.0})
    assert_raises(ValueError, QuadraticGlobalTool, {-2: -5.5, -1: -6.0, -3: -7.0})
    assert_raises(ValueError, QuadraticGlobalTool, {0.0: -5.5, 0.5: -6.0, 1.0: -7.0})
    # check invalid N
    model = QuadraticGlobalTool({5.0: 5.0, 6.0: 10.0, 4.0: 8.0})
    assert_raises(ValueError, model.energy, -0.0001)
    assert_raises(ValueError, model.energy, -1.7)
    assert_raises(ValueError, model.energy, -2.5)
    assert_raises(ValueError, model.energy_derivative, -0.025, 1)
    assert_raises(ValueError, model.energy_derivative, -1.91, 2)
    # check invalid derivative order
    assert_raises(ValueError, model.energy_derivative, 5.0, 1.)
    assert_raises(ValueError, model.energy_derivative, 5.0, 0.2)
    assert_raises(ValueError, model.energy_derivative, 5.0, -1)
    assert_raises(ValueError, model.energy_derivative, 5.0, -3)
    assert_raises(ValueError, model.energy_derivative, 5, '1')
    assert_raises(ValueError, model.energy_derivative, 5, [1])
    assert_raises(ValueError, model.energy_derivative, 3, 1.1)


def test_global_quadratic_nnp_energy():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    energy, deriv, _ = make_symbolic_quadratic_model(1.0, -25.0, -9.0)
    # build global tool
    model = QuadraticGlobalTool({15: -159.0, 16: -153.0, 14: -163.0})
    # check parameters
    assert_almost_equal(model.params[0], -9.0, decimal=6)
    assert_almost_equal(model.params[1], -25.0, decimal=6)
    assert_almost_equal(model.params[2], 1.0, decimal=6)
    assert_almost_equal(model.n0, 15, decimal=6)
    assert_almost_equal(model.n_max, 12.5, decimal=6)
    # check E(N)
    assert_almost_equal(model.energy(20), energy(20), decimal=6)
    assert_almost_equal(model.energy(10), energy(10), decimal=6)
    assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    assert_almost_equal(model.energy_derivative(20, 2), 2.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 2), 2.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 2), 2.0, decimal=6)
    # check d^nE(N) for n > 2
    assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)


def test_global_quadratic_nnp_energy_reactivity():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    energy, _, _ = make_symbolic_quadratic_model(1.0, -25.0, -9.0)
    # build global tool
    model = QuadraticGlobalTool({15: -159.0, 16: -153.0, 14: -163.0})
    # check ionization potential and electron affinity
    ip = energy(14) - energy(15)
    ea = energy(15) - energy(16)
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.eta, ip - ea, decimal=6)
    assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)


def test_global_quadratic_nnp_grand_potential_n():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    _, _, grand = make_symbolic_quadratic_model(1.0, -25.0, -9.0)
    # build global tool
    model = QuadraticGlobalTool({15: -159.0, 16: -153.0, 14: -163.0})
    # check grand potential (as a function of N)
    assert_almost_equal(model.grand_potential(15), grand(15), decimal=6)
    assert_almost_equal(model.grand_potential(14), grand(14), decimal=6)
    assert_almost_equal(model.grand_potential(16.), grand(16.), decimal=6)
    assert_almost_equal(model.grand_potential(15.2), grand(15.2), decimal=6)
    assert_almost_equal(model.grand_potential(14.62), grand(14.62), decimal=6)
    assert_almost_equal(model.grand_potential(11.5), grand(11.5), decimal=6)
    # check grand potential derivative (as a function of N)
    assert_almost_equal(model.grand_potential_derivative(14.), -14, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(15.), -15, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(16.), -16, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(15.001), -15.001, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(14.67), -14.67, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(16.91), -16.91, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(model.n0, 1), -15, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(16., 1), -16., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(17.125, 1), -17.125, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(model.n0, 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(15.89, 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(14.03, 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(16.51, 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(model.n0, 3), 0.0, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(15, 4), 0.0, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(17.5, 5), 0.0, decimal=6)


def test_global_quadratic_nnp_grand_potential_mu():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    _, deriv, grand = make_symbolic_quadratic_model(1.0, -25.0, -9.0)
    # build global tool
    model = QuadraticGlobalTool({15: -159.0, 16: -153.0, 14: -163.0})
    # check mu to N conversion
    assert_almost_equal(model.convert_mu_to_n(5.0), 15., decimal=6)
    assert_almost_equal(model.convert_mu_to_n(5.004), 15.002, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(3.472), 14.236, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(8.962), 16.981, decimal=6)
    # check grand potential (as a function of mu)
    assert_almost_equal(model.grand_potential_mu(deriv(15)), grand(15), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(14)), grand(14), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(16)), grand(16), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(13.7)), grand(13.7), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(14.8)), grand(14.8), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(12.3)), grand(12.3), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(11.4)), grand(11.4), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(15.05)), grand(15.05), decimal=6)
    # check grand potential derivative (as a function of mu)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(15.05), 1), -15.05, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(16.34), 1), -16.34, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(15.61), 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(16.67), 2), -0.5, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(14.31), 3), 0.0, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(12.67), 4), 0.0, decimal=6)


def test_global_quadratic_nnp_grand_potential_reactivity():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    model = QuadraticGlobalTool({15: -159.0, 16: -153.0, 14: -163.0})
    # check hyper-softnesses
    assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_global_quadratic_pnp_energy():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    energy, deriv, _ = make_symbolic_quadratic_model(3.0, -6.0, 30.0)
    # build global tool
    model = QuadraticGlobalTool({5: 75.0, 6: 102.0, 4: 54.0})
    # check parameters
    assert_almost_equal(model.params[0], 30.0, decimal=6)
    assert_almost_equal(model.params[1], -6.0, decimal=6)
    assert_almost_equal(model.params[2], 3.0, decimal=6)
    assert_almost_equal(model.n0, 5, decimal=6)
    assert_almost_equal(model.n_max, 1, decimal=6)
    # check E(N)
    assert_almost_equal(model.energy(20), energy(20), decimal=6)
    assert_almost_equal(model.energy(10), energy(10), decimal=6)
    assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    assert_almost_equal(model.energy_derivative(20, 2), 6.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 2), 6.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 2), 6.0, decimal=6)
    # check d^nE(N) for n > 2
    assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)


def test_global_quadratic_pnp_energy_reactivity():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    energy, _, _ = make_symbolic_quadratic_model(3.0, -6.0, 30.0)
    # build global tool
    model = QuadraticGlobalTool({5: 75.0, 6: 102.0, 4: 54.0})
    # check ionization potential and electron affinity
    ip = energy(4) - energy(5)
    ea = energy(5) - energy(6)
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.eta, ip - ea, decimal=6)
    assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)


def test_global_quadratic_pnp_grand_potential_n():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    _, _, grand = make_symbolic_quadratic_model(3.0, -6.0, 30.0)
    # build global tool
    model = QuadraticGlobalTool({5: 75.0, 6: 102.0, 4: 54.0})
    # check grand potential (as a function of N)
    assert_almost_equal(model.grand_potential(5.), grand(5.), decimal=6)
    assert_almost_equal(model.grand_potential(5.75), grand(5.75), decimal=6)
    assert_almost_equal(model.grand_potential(6.3), grand(6.3), decimal=6)
    assert_almost_equal(model.grand_potential(4.21), grand(4.21), decimal=6)
    assert_almost_equal(model.grand_potential(10.), grand(10.), decimal=6)
    assert_almost_equal(model.grand_potential(15.6), grand(15.6), decimal=6)
    # check grand potential derivative (as a function of N)
    assert_almost_equal(model.grand_potential_derivative(4.), -4, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5.), -5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6., 1), -6, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(3.0123), -3.0123, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(7.2, 1), -7.2, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(8.1, 2), -1 / 6., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6.12, 2), -1 / 6., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5.1, 3), 0., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6.3, 4), 0., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(7.45, 5), 0., decimal=6)


def test_global_quadratic_pnp_grand_potential_mu():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    _, deriv, grand = make_symbolic_quadratic_model(3.0, -6.0, 30.0)
    # build global tool
    model = QuadraticGlobalTool({5: 75.0, 6: 102.0, 4: 54.0})
    # check mu to N conversion
    assert_almost_equal(model.convert_mu_to_n(24.0), 5., decimal=6)
    assert_almost_equal(model.convert_mu_to_n(30.06), 6.01, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(20.16), 4.36, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(24.15), 5.025, decimal=6)
    # check grand potential (as a function of mu)
    assert_almost_equal(model.grand_potential_mu(deriv(5.15)), grand(5.15), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(3.2)), grand(3.2), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(7.67)), grand(7.67), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(5.15)), grand(5.15), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(6.31)), grand(6.31), decimal=6)
    # check grand potential derivative (as a function of mu)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81), 1), -5.81, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.2), 1), -4.2, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81), 2), -1/6., decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.89), 2), -1/6., decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(6.79), 3), 0., decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(3.12), 4), 0., decimal=6)


def test_global_quadratic_pnp_grand_potential_reactivity():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    model = QuadraticGlobalTool({5: 75.0, 6: 102.0, 4: 54.0})
    # check hyper-softnesses
    assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_global_quadratic_n0p_energy():
    # E(N) = -100 + 5*N^2, N0=5
    energy, deriv, _ = make_symbolic_quadratic_model(5.0, 0., -100.)
    # build global tool
    model = QuadraticGlobalTool({5: 25.0, 6: 80.0, 4: -20.0})
    # check parameters
    assert_almost_equal(model.params[0], -100.0, decimal=6)
    assert_almost_equal(model.params[1], 0.0, decimal=6)
    assert_almost_equal(model.params[2], 5.0, decimal=6)
    assert_almost_equal(model.n0, 5, decimal=6)
    assert_almost_equal(model.n_max, 0.0, decimal=6)
    # check E(N)
    assert_almost_equal(model.energy(20), energy(20), decimal=6)
    assert_almost_equal(model.energy(10), energy(10), decimal=6)
    assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    assert_almost_equal(model.energy_derivative(20, 2), 10.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 2), 10.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 2), 10.0, decimal=6)
    # check d^nE(N) for n > 2
    assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)


def test_global_quadratic_n0p_energy_reactivity():
    # E(N) = -100 + 5*N^2, N0=5
    energy, _, _ = make_symbolic_quadratic_model(5.0, 0., -100.)
    ip = energy(4) - energy(5)
    ea = energy(5) - energy(6)
    # build global tool
    model = QuadraticGlobalTool({5: 25.0, 6: 80.0, 4: -20.0})
    # check ionization potential and electron affinity
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    assert_almost_equal(model.eta, ip - ea, decimal=6)
    assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)


def test_global_quadratic_n0p_grand_potential_n():
    # E(N) = -100 + 5*N^2, N0=5
    _, _, grand = make_symbolic_quadratic_model(5.0, 0., -100.)
    # build global tool
    model = QuadraticGlobalTool({5: 25.0, 6: 80.0, 4: -20.0})
    # check grand potential (as a function of N)
    assert_almost_equal(model.grand_potential(5.), grand(5.), decimal=6)
    assert_almost_equal(model.grand_potential(4.), grand(4.), decimal=6)
    assert_almost_equal(model.grand_potential(6.), grand(6.), decimal=6)
    assert_almost_equal(model.grand_potential(4.678), grand(4.678), decimal=6)
    assert_almost_equal(model.grand_potential(6.123), grand(6.123), decimal=6)
    # check grand potential derivative (as a function of N)
    assert_almost_equal(model.grand_potential_derivative(5.), -5, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(4., 1), -4, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6.), -6, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5.6001), -5.6001, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(4.145, 1), -4.145, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5., 2), -0.1, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(7.001, 2), -0.1, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5., 3), 0., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(1.25, 4), 0., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6.5, 5), 0., decimal=6)


def test_global_quadratic_n0p_grand_potential_mu():
    # E(N) = -100 + 5*N^2, N0=5
    _, deriv, grand = make_symbolic_quadratic_model(5.0, 0., -100.)
    # build global tool
    model = QuadraticGlobalTool({5: 25.0, 6: 80.0, 4: -20.0})
    # check mu to N conversion
    assert_almost_equal(model.convert_mu_to_n(50), 5., decimal=6)
    assert_almost_equal(model.convert_mu_to_n(34.6), 3.46, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(69.8), 6.98, decimal=6)
    assert_almost_equal(model.convert_mu_to_n(100.0), 10., decimal=6)
    # check grand potential (as a function of mu)
    assert_almost_equal(model.grand_potential_mu(deriv(5.)), grand(5.), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(5.641)), grand(5.641), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(4.56)), grand(4.56), decimal=6)
    assert_almost_equal(model.grand_potential_mu(deriv(3.905)), grand(3.905), decimal=6)
    # check grand potential derivative (as a function of mu)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81), 1), -5.81, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.341), 1), -4.341, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(6.452), 1), -6.452, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(3.678), 2), -0.1, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.341), 2), -0.1, decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(2.001), 3), 0., decimal=6)
    assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.456), 4), 0., decimal=6)


def test_global_quadratic_n0p_grand_potential_reactivity():
    # E(N) = -100 + 5*N^2, N0=5
    model = QuadraticGlobalTool({5: 25.0, 6: 80.0, 4: -20.0})
    # check hyper-softnesses
    assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_local_quadratic_raises():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # check value of N keys
    assert_raises(ValueError, QuadraticLocalTool, {0.45: d0, 1.0: dp, 0.0: dm})
    assert_raises(ValueError, QuadraticLocalTool, {1.0: d0, 1.5: dp, 0.5: dm})
    assert_raises(ValueError, QuadraticLocalTool, {0.45: d0, 1.0: dp, 0.0: dm})
    assert_raises(ValueError, QuadraticLocalTool, {3.0: d0, 5.0: dp, 1.0: dm})
    # check number of items
    assert_raises(ValueError, QuadraticLocalTool, {2.0: d0, 3.0: dp})
    assert_raises(ValueError, QuadraticLocalTool, {2.0: d0, 3.0: dp, 1.0: dm, 4.0: 2 * d0})
    # check array size
    assert_raises(ValueError, QuadraticLocalTool, {2.0: d0, 3.0: dp, 1.0: dm[:-1]})
    assert_raises(ValueError, QuadraticLocalTool, {2.0: d0, 3.0: dp[:-2], 1.0: dm[:-1]})
    assert_raises(ValueError, QuadraticLocalTool, {2.0: d0[:-1], 3.0: dp, 1.0: dm[:-1]})
    # check invalid Nmax
    assert_raises(ValueError, QuadraticLocalTool, {2.: d0, 3.: dp, 1.: dm}, -1.23)
    assert_raises(ValueError, QuadraticLocalTool, {2.: d0, 3.: dp, 1.: dm}, -0.91, 2.03)


def test_local_quadratic_first_order():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model
    model = QuadraticLocalTool({5: d0, 6: dp, 4: dm})
    # check density
    assert_almost_equal(model.density(5.), d0, decimal=6)
    assert_almost_equal(model.density(6.), dp, decimal=6)
    assert_almost_equal(model.density(4.), dm, decimal=6)
    # check density
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    assert_almost_equal(model.density_derivative(5., 1), expected, decimal=6)
    expected = np.array([-0.75, 2.75, 0.5, -1.5, -2.5])
    assert_almost_equal(model.density_derivative(6., 1), expected, decimal=6)
    expected = np.array([0.25, -2.25, 2.5, 0.5, -0.5])
    assert_almost_equal(model.density_derivative(4., 1), expected, decimal=6)
    expected = np.array([-0.3, 0.5, 1.4, -0.6, -1.6])
    assert_almost_equal(model.density_derivative(5.1, 1), expected, decimal=6)
    expected = np.array([-0.2, 0.0, 1.6, -0.4, -1.4])
    assert_almost_equal(model.density_derivative(4.9, 1), expected, decimal=6)


def test_local_quadratic_higher_order():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build quadratic local model
    model = QuadraticLocalTool({5: d0, 6: dp, 4: dm})
    # check dual Descriptor
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0])
    assert_almost_equal(model.dual_descriptor, expected, decimal=6)
    # expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 3.5
    # assert_almost_equal(model.softness(3.5), expected, decimal=6)
    # expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 0.5
    # assert_almost_equal(model.softness(0.5), expected, decimal=6)
    # expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (1.5 * 1.5)
    # assert_almost_equal(model.hyper_softness(1.5), expected, decimal=6)
    # expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (0.25 * 0.25)
    # assert_almost_equal(model.hyper_softness(0.25), expected, decimal=6)


def test_condensed_quadratic_raises():
    # fake population array
    pop_0 = np.array([-0.5, 1.0, -0.5, 0.5])
    pop_m = np.array([1.0, -2.0, -0.9, 0.3])
    pop_p = np.array([0.8, -3.1, -0.1, 3.6])
    # check dict_pop
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p, 4.: pop_p})
    assert_raises(ValueError, QuadraticCondensedTool, {5.: pop_m, 6: pop_0, 7.5: pop_p})
    assert_raises(ValueError, QuadraticCondensedTool, {-1.: pop_m, 0.: pop_0, 1.: pop_p})
    assert_raises(ValueError, QuadraticCondensedTool, {0.: pop_m, 1.: pop_0, 3.: pop_p})
    assert_raises(ValueError, QuadraticCondensedTool, {0.: pop_m, 0.5: pop_0, 1.: pop_p})
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p[:-1]})
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m[:-1], 2.: pop_0, 3.: pop_p[:-1]})
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m, 2.: pop_0[:-1], 3.: pop_p[:-1]})
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m[:-1], 2.: pop_0, 3.: pop_p})
    # check N_max
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p}, -0.5, 1.5)
    assert_raises(ValueError, QuadraticCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p}, -2., 1.5)
    # make a model
    model = QuadraticCondensedTool({1.: pop_m, 2.: pop_0, 3.: pop_p})
    assert_raises(ValueError, model.population, -1.5)
    assert_raises(ValueError, model.population, "2.0")
    assert_raises(ValueError, model.population_derivative, -0.5)
    assert_raises(ValueError, model.population_derivative, "1.6")
    assert_raises(ValueError, model.population_derivative, 1., order=0.5)
    assert_raises(ValueError, model.population_derivative, 1., order=-1)
    assert_raises(ValueError, model.population_derivative, 1., order=1.)


def test_condensed_linear_h2o_population():
    # ESP charges of H2O at ub3lyp/ccpvtz
    pop_09 = np.array([4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    pop_10 = np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    pop_11 = np.array([-5.81613550E-01, -2.09193820E-01, -2.09192630E-01])
    # compute fukui function & dual descriptor
    ff, df = 0.5 * (pop_11 - pop_09), pop_11 - 2 * pop_10 + pop_09
    # build model
    model = QuadraticCondensedTool({9: pop_09, 10: pop_10, 11: pop_11}, None, None)
    # check N0, Nmax & S
    assert_equal(model.n_ref, 10)
    assert model.n_max is None
    assert model.global_softness is None
    # check population
    assert_almost_equal(model.population(8.9), pop_10 - ff * 1.1 + 0.5 * df * 1.21, decimal=8)
    assert_almost_equal(model.population(9.0), pop_09, decimal=8)
    assert_almost_equal(model.population(9.7), pop_10 - ff * 0.3 + 0.5 * df * 0.09, decimal=8)
    assert_almost_equal(model.population(10.0), pop_10, decimal=8)
    assert_almost_equal(model.population(10.5), pop_10 + ff * 0.5 + 0.5 * df * 0.25, decimal=8)
    assert_almost_equal(model.population(11.), pop_11, decimal=8)
    assert_almost_equal(model.population(11.1), pop_10 + ff * 1.1 + 0.5 * df * 1.21, decimal=8)
    # check population derivative
    assert_almost_equal(model.population_derivative(7.1), ff - 2.9 * df, decimal=8)
    assert_almost_equal(model.population_derivative(9.0), ff - 1.0 * df, decimal=8)
    assert_almost_equal(model.population_derivative(9.2), ff - 0.8 * df, decimal=8)
    assert_almost_equal(model.population_derivative(10.8), ff + 0.8 * df, decimal=8)
    assert_almost_equal(model.population_derivative(11.), ff + 1.0 * df, decimal=8)
    assert_almost_equal(model.population_derivative(11.7), ff + 1.7 * df, decimal=8)


def test_condensed_linear_h2o_reactivity():
    # ESP charges of H2O at ub3lyp/ccpvtz
    pop_09 = np.array([4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    pop_10 = np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    pop_11 = np.array([-5.81613550E-01, -2.09193820E-01, -2.09192630E-01])
    # build model
    model = QuadraticCondensedTool({9: pop_09, 10: pop_10, 11: pop_11}, 10., 2.3)
    # check N0, Nmax & S
    assert_equal(model.n_ref, 10)
    assert_equal(model.n_max, 10)
    assert_equal(model.global_softness, 2.3)
    # check reactivity indicators
    assert_almost_equal(model.fukui_function, 0.5 * (pop_11 - pop_09), decimal=8)
    assert_almost_equal(model.dual_descriptor, pop_11 - 2 * pop_10 + pop_09, decimal=8)
    assert_almost_equal(model.softness, 2.3 * 0.5 * (pop_11 - pop_09), decimal=8)
