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
"""Test chemtools.conceptual.quadratic Module."""

import numpy as np
from chemtools.conceptual.quadratic import QuadraticGlobalTool, QuadraticLocalTool


def test_global_quadratic1():
    # E(N) = -9.0 + (-25.0)*N + N^2, N0=15
    model = QuadraticGlobalTool(-159.0, -153.0, -163.0, 15)
    # check parameters
    np.testing.assert_almost_equal(model.params[0], -9.0, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -25.0, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 1.0, decimal=6)
    np.testing.assert_almost_equal(model.n0, 15, decimal=6)
    np.testing.assert_almost_equal(model.n_max, 12.5, decimal=6)
    # check E(N)
    energy = lambda n: -9.0 - 25.0 * n + n * n
    np.testing.assert_almost_equal(model.energy(20), energy(20), decimal=6)
    np.testing.assert_almost_equal(model.energy(10), energy(10), decimal=6)
    np.testing.assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    deriv = lambda n: -25.0 + 2 * n
    np.testing.assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    np.testing.assert_almost_equal(model.energy_derivative(20, 2), 2.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 2), 2.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 2), 2.0, decimal=6)
    # check d^nE(N) for n > 2
    np.testing.assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)
    # check ionization potential and electron affinity
    ip = energy(14) - energy(15)
    ea = energy(15) - energy(16)
    np.testing.assert_almost_equal(model.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, ea, decimal=6)
    np.testing.assert_almost_equal(model.ip, ip, decimal=6)
    np.testing.assert_almost_equal(model.ea, ea, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    np.testing.assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    np.testing.assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.eta, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)
    # check grand potential (as a function of N)
    grand = lambda n: energy(n) - deriv(n) * n
    np.testing.assert_almost_equal(model.grand_potential(15), grand(15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(14), grand(14), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(16.), grand(16.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(15.2), grand(15.2), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(14.62), grand(14.62), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(11.5), grand(11.5), decimal=6)
    # check grand potential derivative (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential_derivative(14.), -14, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(15.), -15, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(16.), -16, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(15.001), -15.001, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(14.67), -14.67, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(16.91), -16.91, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(model.n0, 1), -15, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(16., 1), -16., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(17.125, 1), -17.125, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(model.n0, 2), -0.5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(15.89, 2), -0.5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(14.03, 2), -0.5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(16.51, 2), -0.5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(model.n0, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(15, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(17.5, 5), 0.0, decimal=6)
    # check mu to N conversion
    np.testing.assert_almost_equal(model.convert_mu_to_n(5.0), 15., decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(5.004), 15.002, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(3.472), 14.236, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(8.962), 16.981, decimal=6)
    # check grand potential (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(15)), grand(15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(14)), grand(14), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(16)), grand(16), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(13.7)), grand(13.7), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(14.8)), grand(14.8), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(12.3)), grand(12.3), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(11.4)), grand(11.4), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(15.05)), grand(15.05), decimal=6)
    # check grand potential derivative (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(15.05)), -15.05,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(16.34), 1), -16.34,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(15.61), 2), -0.5,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(16.67), 2), -0.5,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(14.31), 3), 0.0,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(12.67), 4), 0.0,
                                   decimal=6)
    # check hyper-softnesses
    np.testing.assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_global_quadratic2():
    # E(N) = 30.0 + (-6.0)*N + 3*N^2, N0=10
    model = QuadraticGlobalTool(75.0, 102.0, 54.0, 5)
    # check parameters
    np.testing.assert_almost_equal(model.params[0], 30.0, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -6.0, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 3.0, decimal=6)
    np.testing.assert_almost_equal(model.n0, 5, decimal=6)
    np.testing.assert_almost_equal(model.n_max, 1, decimal=6)
    # check E(N)
    energy = lambda n: 30.0 - 6.0 * n + 3 * n * n
    np.testing.assert_almost_equal(model.energy(20), energy(20), decimal=6)
    np.testing.assert_almost_equal(model.energy(10), energy(10), decimal=6)
    np.testing.assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    deriv = lambda n: -6.0 + 6 * n
    np.testing.assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    np.testing.assert_almost_equal(model.energy_derivative(20, 2), 6.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 2), 6.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 2), 6.0, decimal=6)
    # check d^nE(N) for n > 2
    np.testing.assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)
    # check ionization potential and electron affinity
    ip = energy(4) - energy(5)
    ea = energy(5) - energy(6)
    np.testing.assert_almost_equal(model.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, ea, decimal=6)
    np.testing.assert_almost_equal(model.ip, ip, decimal=6)
    np.testing.assert_almost_equal(model.ea, ea, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    np.testing.assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    np.testing.assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.eta, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)
    # check grand potential (as a function of N)
    grand = lambda n: energy(n) - deriv(n) * n
    np.testing.assert_almost_equal(model.grand_potential(5.), grand(5.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(5.75), grand(5.75), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(6.3), grand(6.3), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(4.21), grand(4.21), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(10.), grand(10.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(15.6), grand(15.6), decimal=6)
    # check grand potential derivative (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.), -4, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5.), -5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6., 1), -6, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(3.0123), -3.0123, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(7.2, 1), -7.2, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(8.1, 2), -1 / 6., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.12, 2), -1 / 6., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5.1, 3), 0., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.3, 4), 0., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(7.45, 5), 0., decimal=6)
    # check mu to N conversion
    np.testing.assert_almost_equal(model.convert_mu_to_n(24.0), 5., decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(30.06), 6.01, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(20.16), 4.36, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(24.15), 5.025, decimal=6)
    # check grand potential (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(5.15)), grand(5.15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(3.2)), grand(3.2), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(7.67)), grand(7.67), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(5.15)), grand(5.15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(6.31)), grand(6.31), decimal=6)
    # check grand potential derivative (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81)), -5.81,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.2), 1), -4.2,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81), 2), -1/6.,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.89), 2), -1/6.,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(6.79), 3), 0.,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(3.12), 4), 0.,
                                   decimal=6)
    # check hyper-softnesses
    np.testing.assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_global_quadratic3():
    # E(N) = -100 + 5*N^2, N0=5
    model = QuadraticGlobalTool(25.0, 80.0, -20.0, 5)
    # check parameters
    np.testing.assert_almost_equal(model.params[0], -100.0, decimal=6)
    np.testing.assert_almost_equal(model.params[1], 0.0, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 5.0, decimal=6)
    np.testing.assert_almost_equal(model.n0, 5, decimal=6)
    np.testing.assert_almost_equal(model.n_max, 0.0, decimal=6)
    # check E(N)
    energy = lambda n: -100.0 + 5 * n * n
    np.testing.assert_almost_equal(model.energy(20), energy(20), decimal=6)
    np.testing.assert_almost_equal(model.energy(10), energy(10), decimal=6)
    np.testing.assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check dE(N)
    deriv = lambda n: 10 * n
    np.testing.assert_almost_equal(model.energy_derivative(20), deriv(20), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10), deriv(10), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5), deriv(16.5), decimal=6)
    # check d2E(N)
    np.testing.assert_almost_equal(model.energy_derivative(20, 2), 10.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 2), 10.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 2), 10.0, decimal=6)
    # check d^nE(N) for n > 2
    np.testing.assert_almost_equal(model.energy_derivative(20, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 5), 0.0, decimal=6)
    # check ionization potential and electron affinity
    ip = energy(4) - energy(5)
    ea = energy(5) - energy(6)
    np.testing.assert_almost_equal(model.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, ea, decimal=6)
    np.testing.assert_almost_equal(model.ip, ip, decimal=6)
    np.testing.assert_almost_equal(model.ea, ea, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    electrophil = (-0.5 * (ip + ea))**2 / (2 * (ip - ea))
    np.testing.assert_almost_equal(model.electrophilicity, -electrophil, decimal=6)
    nucleofugal = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.nucleofugality, nucleofugal, decimal=6)
    electrofugal = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(model.electrofugality, -electrofugal, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    np.testing.assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.eta, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)
    # check grand potential (as a function of N)
    grand = lambda n: energy(n) - deriv(n) * n
    np.testing.assert_almost_equal(model.grand_potential(5.), grand(5.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(4.), grand(4.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(6.), grand(6.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(4.678), grand(4.678), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(6.123), grand(6.123), decimal=6)
    # check grand potential derivative (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5.), -5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4., 1), -4, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.), -6, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5.6001), -5.6001, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.145, 1), -4.145, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5., 2), -0.1, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(7.001, 2), -0.1, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5., 3), 0., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(1.25, 4), 0., decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5, 5), 0., decimal=6)
    # check mu to N conversion
    np.testing.assert_almost_equal(model.convert_mu_to_n(50), 5., decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(34.6), 3.46, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(69.8), 6.98, decimal=6)
    np.testing.assert_almost_equal(model.convert_mu_to_n(100.0), 10., decimal=6)
    # check grand potential (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(5.)), grand(5.), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(5.641)), grand(5.641), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(4.56)), grand(4.56), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu(deriv(3.905)), grand(3.905), decimal=6)
    # check grand potential derivative (as a function of mu)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.81)), -5.81,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.341)), -4.341,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(6.452), 1), -6.452,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(3.678), 2), -0.1,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(4.341), 2), -0.1,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(2.001), 3), 0.,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(deriv(5.456), 4), 0.,
                                   decimal=6)
    # check hyper-softnesses
    np.testing.assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)


def test_local_quadratic():
    # Fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # Build a linear local model
    model = QuadraticLocalTool(d0, dp, dm, 5)
    # Check density
    np.testing.assert_almost_equal(model.density_zero, d0, decimal=6)
    np.testing.assert_almost_equal(model.density_plus, dp, decimal=6)
    np.testing.assert_almost_equal(model.density_minus, dm, decimal=6)
    np.testing.assert_almost_equal(model.density(5.), d0, decimal=6)
    np.testing.assert_almost_equal(model.density(6.), dp, decimal=6)
    np.testing.assert_almost_equal(model.density(4.), dm, decimal=6)
    np.testing.assert_almost_equal(model.density(None), d0, decimal=6)
    # Density
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    np.testing.assert_almost_equal(model.fukui_function(None), expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(5.), expected, decimal=6)
    expected = np.array([-0.75, 2.75, 0.5, -1.5, -2.5])
    np.testing.assert_almost_equal(model.fukui_function(6.), expected, decimal=6)
    expected = np.array([0.25, -2.25, 2.5, 0.5, -0.5])
    np.testing.assert_almost_equal(model.fukui_function(4.), expected, decimal=6)
    expected = np.array([-0.3, 0.5, 1.4, -0.6, -1.6])
    np.testing.assert_almost_equal(model.fukui_function(5.1), expected, decimal=6)
    expected = np.array([-0.2, 0.0, 1.6, -0.4, -1.4])
    np.testing.assert_almost_equal(model.fukui_function(4.9), expected, decimal=6)
    # Dual Descriptor
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0])
    np.testing.assert_almost_equal(model.dual_descriptor(), expected, decimal=6)
    # expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 3.5
    # np.testing.assert_almost_equal(model.softness(3.5), expected, decimal=6)
    # expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 0.5
    # np.testing.assert_almost_equal(model.softness(0.5), expected, decimal=6)
    # expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (1.5 * 1.5)
    # np.testing.assert_almost_equal(model.hyper_softness(1.5), expected, decimal=6)
    # expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (0.25 * 0.25)
    # np.testing.assert_almost_equal(model.hyper_softness(0.25), expected, decimal=6)
