# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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
"""Test chemtools.conceptual.general Module."""

import math
import numpy as np
import sympy as sp
from chemtools.conceptual.general import GeneralGlobalTool


def test_global_general_energy_quadratic():
    # E(N) = 31.0 - 28.0 * N + 4.0 * N^2
    n, n0, a, b, c = sp.symbols('n, n0, a, b, c')
    expr = a + b * n + c * (n**2)
    model = GeneralGlobalTool(expr, 3.45, {2.1: -10.16, 2.5: -14.0, 4.3: -15.44}, n, n0)
    np.testing.assert_almost_equal(model.params[a], 31.0, decimal=6)
    np.testing.assert_almost_equal(model.params[b], -28.0, decimal=6)
    np.testing.assert_almost_equal(model.params[c], 4.0, decimal=6)
    np.testing.assert_almost_equal(model.n0, 3.45, decimal=6)
    np.testing.assert_almost_equal(model.n_max, 28. / (2 * 4.0), decimal=6)
    # check energy
    energy = lambda n: 31.0 - 28.0 * n + 4.0 * (n**2)
    np.testing.assert_almost_equal(model.energy(5.23), energy(5.23), decimal=6)
    np.testing.assert_almost_equal(model.energy(3.45), energy(3.45), decimal=6)
    np.testing.assert_almost_equal(model.energy(3.00), energy(3.0), decimal=6)
    # check energy derivatives
    deriv = lambda n: -28.0 + 8.0 * n
    np.testing.assert_almost_equal(model.energy_derivative(6.30), deriv(6.30), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3.45), deriv(3.45), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.70), deriv(1.70), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.55, 2), 8.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3.45, 2), 8.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(16.5, 2), 8.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(9.20, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3.45, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.00, 5), 0.0, decimal=6)
    # check ionization potential and electron affinity
    ip = energy(2.45) - energy(3.45)
    ea = energy(3.45) - energy(4.45)
    np.testing.assert_almost_equal(model.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, ea, decimal=6)
    np.testing.assert_almost_equal(model.ip, ip, decimal=6)
    np.testing.assert_almost_equal(model.ea, ea, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.chemical_potential, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.5 * (ip + ea), decimal=6)
    np.testing.assert_almost_equal(model.eta, ip - ea, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0 / (ip - ea), decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, (ip + ea)**2 / (8*(ip - ea)), decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, (ip - 3*ea)**2 / (8*(ip - ea)), decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, (3*ip - ea)**2 / (8*(ip - ea)), decimal=6)
    # check grand potential
    grand = lambda n: energy(n) - deriv(n) * n
    np.testing.assert_almost_equal(model.grand_potential(15), grand(15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(10), grand(10), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(20), grand(20), decimal=6)


def test_global_general_energy_exponential():
    # E(N) = 6.91 * exp(-0.25 * (N - 7.0)) + 2.74
    n, n0, a, b, gamma = sp.symbols('n, n0, A, B, gamma')
    expr = a * sp.exp(- gamma * (n - 7.0)) + b
    n_energies = {7.5: 8.838053596859556, 1.25: 31.832186639954763, 3.6: 18.906959746808596}
    model = GeneralGlobalTool(expr, 7.0, n_energies, n, n0)
    np.testing.assert_almost_equal(model.params[a], 6.91, decimal=6)
    np.testing.assert_almost_equal(model.params[b], 2.74, decimal=6)
    np.testing.assert_almost_equal(model.params[gamma], 0.25, decimal=6)
    np.testing.assert_almost_equal(model.n0, 7, decimal=6)
    np.testing.assert_almost_equal(model.n_max, float('inf'), decimal=6)
    # check energy
    energy = lambda n: 6.91 * math.exp(-0.25 * (n - 7.0)) + 2.74
    np.testing.assert_almost_equal(model.energy(20), energy(20), decimal=6)
    np.testing.assert_almost_equal(model.energy(10), energy(10), decimal=6)
    np.testing.assert_almost_equal(model.energy(8), energy(8), decimal=6)
    np.testing.assert_almost_equal(model.energy(16.5), energy(16.5), decimal=6)
    # check energy derivatives
    dE = lambda n, r: 6.91 * math.pow(-0.25, r) * math.exp(-0.25 * (n - 7))
    np.testing.assert_almost_equal(model.energy_derivative(18.1), dE(18.1, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10), dE(10, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(8.5), dE(8.5, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(12.25), dE(12.25, 1), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(17.3, 2), dE(17.3, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 2), dE(10, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(9.1, 2), dE(9.1, 2), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(13.7, 3), dE(13.7, 3), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10, 4), dE(10, 4), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4.5, 5), dE(4.5, 5), decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(6.5, 10), dE(6.5, 10), decimal=6)
    # check ionization potential and electron affinity
    ip = energy(6) - energy(7)
    ea = energy(7) - energy(8)
    np.testing.assert_almost_equal(model.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, ea, decimal=6)
    np.testing.assert_almost_equal(model.ip, ip, decimal=6)
    np.testing.assert_almost_equal(model.ea, ea, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, -dE(7, 1), decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, 6.91, decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, -(energy(8.0) - 2.74), decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, energy(6.0) - 2.74, decimal=6)
    # check chemical potential, chemical hardness, and related tools
    np.testing.assert_almost_equal(model.chemical_potential, dE(7, 1), decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, dE(7, 2), decimal=6)
    np.testing.assert_almost_equal(model.mu, dE(7, 1), decimal=6)
    np.testing.assert_almost_equal(model.eta, dE(7, 2), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), dE(7, 3), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), dE(7, 4), decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), dE(7, 5), decimal=6)
    np.testing.assert_almost_equal(model.softness, 1.0/dE(7, 2), decimal=6)
    np.testing.assert_almost_equal(model.hyper_softness(2), -dE(7., 3)/dE(7, 2)**3, decimal=6)
    # check grand potential
    grand = lambda n: energy(n) - dE(n, 1) * n
    np.testing.assert_almost_equal(model.grand_potential(15), grand(15), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(10), grand(10), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(20), grand(20), decimal=6)


# def test_global_general_morse():
    # # E(N) = 4.01 * exp(-0.17 * (N - 17.0) + 0.32) + 6.95
    # n, n0, a, b, gamma, delta = sp.symbols('n, n0, A, B, gamma, delta')
    # expr = a * sp.exp(- gamma * (n - 4.0) + delta) + b
    # # energies = {1.0: 16.146208148459372, 4.0:12.472282334987188, 5.9:10.947988026968526,
    # #             7.7:9.894064887805316}
    # energies = {0.5:98.21717932504218, 4.2:55.606783873132294, 5.9:43.39452499661585,
    #             7.7:33.787260559941295}
    # model = GeneralGlobalTool(expr, 17.0, energies, n, n0)
    # print model._params
    # np.testing.assert_almost_equal(model._params[a], 4.01, decimal=6)
    # np.testing.assert_almost_equal(model._params[b], 6.95, decimal=6)
    # np.testing.assert_almost_equal(model._params[gamma], 0.17, decimal=6)
    # np.testing.assert_almost_equal(model._params[delta], 0.32, decimal=6)
    # np.testing.assert_almost_equal(model.n0, 7, decimal=6)
    # assert 5 == 6
    # pass
