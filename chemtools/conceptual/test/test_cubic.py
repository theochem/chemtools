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
# pragma pylint: disable=invalid-name
"""Test chemtools.conceptual.cubic Module."""


import numpy as np

from numpy.testing import assert_raises, assert_almost_equal

from chemtools.conceptual.cubic import CubicGlobalTool


def test_global_cubic_raises():
    # check invalid N0
    assert_raises(ValueError, CubicGlobalTool, {5: -5.5, 4: -6.0, 6: -7.0, 7: -7.5})
    assert_raises(ValueError, CubicGlobalTool, {0: -5.5, 1: -6.0, -1: -7.0})
    assert_raises(ValueError, CubicGlobalTool, {0.1: -5.5, 1.1: -6.0, -0.9: -7.0})
    assert_raises(ValueError, CubicGlobalTool, {0.99: -5.5, 1.99: -6.0, -0.1: -7.0})
    assert_raises(ValueError, CubicGlobalTool, {-1.: -5.5, 0.: -6.0, -2.: -7.0})
    assert_raises(ValueError, CubicGlobalTool, {-2: -5.5, -1: -6.0, -3: -7.0})
    assert_raises(ValueError, CubicGlobalTool, {0.0: -5.5, 0.5: -6.0, 1.0: -7.0})
    # check invalid N
    model = CubicGlobalTool({5.0: 5.0, 6.0: 10.0, 4.0: 8.0})
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


def test_global_cubic_omega_half():
    # E(N) = 100. + 12.6 * (N - 10) - 62.1 * (N - 10)**2 + 0. * (N - 10)**3, N0=10
    dict_energy = {9: 25.3, 10: 100., 11: 50.5}
    model = CubicGlobalTool(dict_energy, omega=0.5)
    # check n0, n_max, omega, ip & ea
    assert model.n_max is None
    assert_almost_equal(model.n0, 10, decimal=10)
    assert_almost_equal(model.omega, 0.5, decimal=10)
    assert_almost_equal(model.ip, 25.3 - 100., decimal=10)
    assert_almost_equal(model.ionization_potential, 25.3 - 100., decimal=10)
    assert_almost_equal(model.ea, 100. - 50.5, decimal=10)
    assert_almost_equal(model.electron_affinity, 100. - 50.5, decimal=10)
    # check parameters
    assert_almost_equal(model.params, np.array([100., 12.6, -62.1, 0.0]), decimal=10)
    # check energy values
    assert_almost_equal(model.energy(7.8), -228.284, decimal=10)
    assert_almost_equal(model.energy(9), 25.3, decimal=10)
    assert_almost_equal(model.energy(9.5), 78.175, decimal=10)
    assert_almost_equal(model.energy(10), 100., decimal=10)
    assert_almost_equal(model.energy(10.5), 90.775, decimal=10)
    assert_almost_equal(model.energy(11.), 50.5, decimal=10)
    assert_almost_equal(model.energy(11.9), -100.241, decimal=10)
    # check energy derivatives
    assert_almost_equal(model.energy_derivative(9., 1), 136.8, decimal=10)
    assert_almost_equal(model.energy_derivative(9.5, 1), 74.7, decimal=10)
    assert_almost_equal(model.energy_derivative(10., 1), 12.6, decimal=10)
    assert_almost_equal(model.energy_derivative(10.5, 1), -49.5, decimal=10)
    assert_almost_equal(model.energy_derivative(11., 1), -111.6, decimal=10)
    assert_almost_equal(model.energy_derivative(9., 2), -124.2, decimal=10)
    assert_almost_equal(model.energy_derivative(10., 2), -124.2, decimal=10)
    assert_almost_equal(model.energy_derivative(11.2, 2), -124.2, decimal=10)
    assert_almost_equal(model.energy_derivative(9., 3), 0., decimal=10)
    assert_almost_equal(model.energy_derivative(10., 4), 0., decimal=10)
    assert_almost_equal(model.energy_derivative(11., 5), 0., decimal=10)


def test_global_cubic_omega_half_reactivity():
    # E(N) = 100. + 12.6 * (N - 10) - 62.1 * (N - 10)**2 + 0. * (N - 10)**3, N0=10
    dict_energy = {9: 25.3, 10: 100., 11: 50.5}
    model = CubicGlobalTool(dict_energy, omega=0.5)
    # check reactivity indicators
    assert_almost_equal(model.mu, 12.6, decimal=10)
    assert_almost_equal(model.chemical_potential, 12.6, decimal=10)
    assert_almost_equal(model.eta, -124.2, decimal=10)
    assert_almost_equal(model.chemical_hardness, -124.2, decimal=10)
    assert_almost_equal(model.softness, -1. / 124.2, decimal=10)
    assert_almost_equal(model.hyper_hardness(), 0., decimal=10)
    assert_almost_equal(model.hyper_hardness(2), 0., decimal=10)
    assert_almost_equal(model.hyper_hardness(3), 0., decimal=10)


def test_global_cubic_omega_one():
    # E(N) = -100. - 4.48 * (N - 5) - 13.12 * (N - 5)**2 - 13.12 * (N - 5)**3, N0=5
    dict_energy = {4: -95.52, 5: -100., 6: -130.72}
    model = CubicGlobalTool(dict_energy, omega=1.0)
    # check n0, n_max, omega, ip & ea
    assert model.n_max is None
    assert_almost_equal(model.n0, 5, decimal=10)
    assert_almost_equal(model.omega, 1.0, decimal=10)
    assert_almost_equal(model.ip, -95.52 + 100., decimal=10)
    assert_almost_equal(model.ionization_potential, -95.52 + 100., decimal=10)
    assert_almost_equal(model.ea, -100. + 130.72, decimal=10)
    assert_almost_equal(model.electron_affinity, -100. + 130.72, decimal=10)
    # check parameters
    assert_almost_equal(model.params, np.array([-100., -4.48, -13.12, -13.12]), decimal=10)
    # check energy values
    assert_almost_equal(model.energy(4), -95.52, decimal=10)
    assert_almost_equal(model.energy(4.6), -99.46752, decimal=10)
    assert_almost_equal(model.energy(5), -100., decimal=10)
    assert_almost_equal(model.energy(5.2), -101.52576, decimal=10)
    assert_almost_equal(model.energy(6), -130.72, decimal=10)
    # check energy derivatives
    assert_almost_equal(model.energy_derivative(4.0, 1), -17.6, decimal=10)
    assert_almost_equal(model.energy_derivative(6.2, 1), -92.6464, decimal=10)
    assert_almost_equal(model.energy_derivative(5.0, 2), -26.24, decimal=10)
    assert_almost_equal(model.energy_derivative(5.7, 2), -81.344, decimal=10)
    assert_almost_equal(model.energy_derivative(6.0, 3), -78.72, decimal=10)
    assert_almost_equal(model.energy_derivative(7.3, 3), -78.72, decimal=10)
    assert_almost_equal(model.energy_derivative(6.5, 4), 0., decimal=10)
    assert_almost_equal(model.energy_derivative(7.0, 5), 0., decimal=10)
    assert_almost_equal(model.energy_derivative(7.4, 6), 0., decimal=10)


def test_global_cubic_omega_one_reactivity():
    # E(N) = -100. - 4.48 * (N - 5) - 13.12 * (N - 5)**2 - 13.12 * (N - 5)**3, N0=5
    dict_energy = {4: -95.52, 5: -100., 6: -130.72}
    model = CubicGlobalTool(dict_energy, omega=1.0)
    # check reactivity indicators
    assert_almost_equal(model.mu, -4.48, decimal=10)
    assert_almost_equal(model.chemical_potential, -4.48, decimal=10)
    assert_almost_equal(model.eta, -26.24, decimal=10)
    assert_almost_equal(model.chemical_hardness, -26.24, decimal=10)
    assert_almost_equal(model.softness, -1. / 26.24, decimal=10)
    assert_almost_equal(model.hyper_hardness(), -78.72, decimal=10)
    assert_almost_equal(model.hyper_hardness(2), -78.72, decimal=10)
    assert_almost_equal(model.hyper_hardness(3), 0., decimal=10)
    assert_almost_equal(model.hyper_hardness(4), 0., decimal=10)
