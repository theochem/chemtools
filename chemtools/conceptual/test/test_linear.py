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
"""Test chemtools.conceptual.linear Module."""

import numpy as np
from numpy.testing import assert_raises, assert_equal, assert_almost_equal
from chemtools.conceptual.linear import LinearGlobalTool, LinearLocalTool, LinearCondensedTool


def test_global_linear_raises():
    # check invalid N0
    assert_raises(ValueError, LinearGlobalTool, {0.: -10.0, 1.: -15.0, -1.: -16.0})
    assert_raises(ValueError, LinearGlobalTool, {0.5: -10.0, 1.5: -15.0, -0.5: -16.0})
    assert_raises(ValueError, LinearGlobalTool, {0.9: -10.0, 1.9: -15.0, -0.1: -16.0})
    assert_raises(ValueError, LinearGlobalTool, {-5.0: -10.0, -4.0: -15.0, -6.0: -16.0})
    assert_raises(ValueError, LinearGlobalTool, {0.9: 10.0, 0.5: 15.0, 1.2: 16.0})
    assert_raises(ValueError, LinearGlobalTool, {1.0: 10.0, 0.5: 15.0, 1.5: 16.0})
    # check invalid N
    model = LinearGlobalTool({5.: -10.5, 6.: -11.7, 4.: -12.3})
    assert_raises(ValueError, model.energy, -1.0)
    assert_raises(ValueError, model.energy, -2.3)
    assert_raises(ValueError, model.energy_derivative, -0.5, 1)
    assert_raises(ValueError, model.energy_derivative, -1.9, 2)
    # check invalid derivative order
    assert_raises(ValueError, model.energy_derivative, 0.5, 1.)
    assert_raises(ValueError, model.energy_derivative, 0.5, 0.2)
    assert_raises(ValueError, model.energy_derivative, 0.5, -1)
    assert_raises(ValueError, model.energy_derivative, 0.5, -3)
    assert_raises(ValueError, model.energy_derivative, 0.5, '1')
    assert_raises(ValueError, model.energy_derivative, 0.5, [1])


def test_global_linear_np_energy():
    # E(N) = -1.0 - 0.5 * N, N <= 10
    # E(N) = -7.0 + 0.1 * N, N >= 10
    # build global tool
    model = LinearGlobalTool({10: -6.0, 11: -5.9, 9: -5.5})
    # check parameters
    assert_equal(model.n0, 10)
    assert_equal(model.n_max, 10)
    # check E(N)
    assert_almost_equal(model.energy(10.), -6.0, decimal=6)
    assert_almost_equal(model.energy(11), -5.9, decimal=6)
    assert_almost_equal(model.energy(9.), -5.5, decimal=6)
    assert_almost_equal(model.energy(8.5), -5.25, decimal=6)
    assert_almost_equal(model.energy(11.2), -5.88, decimal=6)
    assert_almost_equal(model.energy(13.56), -5.644, decimal=6)
    assert_almost_equal(model.energy(5.01), -3.505, decimal=6)
    assert_almost_equal(model.energy(0), -1.0, decimal=6)
    assert_almost_equal(model.energy(model.n_max), -6.0, decimal=6)
    # check dE(N)
    assert_equal(model.energy_derivative(10.0, 1), None)
    assert_almost_equal(model.energy_derivative(9.0, 1), -0.5, decimal=6)
    assert_almost_equal(model.energy_derivative(7.25, 1), -0.5, decimal=6)
    assert_almost_equal(model.energy_derivative(11., 1), 0.1, decimal=6)
    assert_almost_equal(model.energy_derivative(10.45, 1), 0.1, decimal=6)
    # check d^nE(N) for n > 1
    assert_equal(model.energy_derivative(10.0, 2), None)
    assert_equal(model.energy_derivative(9.5, 2), 0.)
    assert_equal(model.energy_derivative(10.9, 2), 0.)
    assert_equal(model.energy_derivative(10.0, 3), None)
    assert_equal(model.energy_derivative(8.9, 2), 0.)
    assert_equal(model.energy_derivative(11.3, 2), 0.)


def test_global_linear_np_reactivity():
    # E(N) = -1.0 - 0.5 * N, N <= 10
    # E(N) = -7.0 + 0.1 * N, N >= 10
    # build global tool
    model = LinearGlobalTool({10: -6.0, 11: -5.9, 9: -5.5})
    # check ionization potential and electron affinity
    assert_almost_equal(model.ip, 0.5)
    assert_almost_equal(model.ea, -0.1)
    assert_almost_equal(model.ionization_potential, 0.5)
    assert_almost_equal(model.electron_affinity, -0.1)
    # check fundamental descriptors
    assert_almost_equal(model.mu_minus, -0.5)
    assert_almost_equal(model.mu_plus, 0.1)
    assert_almost_equal(model.mu_zero, -0.2)
    assert_equal(model.mu, None)
    assert_equal(model.chemical_potential, None)
    assert_equal(model.eta, None)
    assert_equal(model.chemical_hardness, None)
    assert_equal(model.electronegativity, None)
    # check derived descriptors
    assert_almost_equal(model.electrofugality, 0.5, decimal=6)
    assert_almost_equal(model.electrophilicity, 0.0, decimal=6)
    assert_almost_equal(model.nucleofugality, 0.1, decimal=6)


def test_global_linear_np_grand_potential_n():
    # E(N) = -1.0 - 0.5 * N, N <= 10
    # E(N) = -7.0 + 0.1 * N, N >= 10
    # build global tool
    model = LinearGlobalTool({10: -6.0, 11: -5.9, 9: -5.5})
    # check grand potential (as a function of N)
    assert_equal(model.grand_potential(None), None)
    assert_equal(model.grand_potential(10.), None)
    assert_almost_equal(model.grand_potential(9.), -1.0, decimal=6)
    assert_almost_equal(model.grand_potential(6.3), -1.0, decimal=6)
    assert_almost_equal(model.grand_potential(11.), -7.0, decimal=6)
    assert_almost_equal(model.grand_potential(10.1), -7.0, decimal=6)
    # check grand potential derivative (as a function of N)
    assert_almost_equal(model.grand_potential_derivative(9.3, 1), -9.3, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.0, 1), -10., decimal=6)
    assert_almost_equal(model.grand_potential_derivative(11.2, 1), -11.2, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(10.7, 1), -10.7, decimal=6)
    assert_equal(model.grand_potential_derivative(10., 2), None)
    assert_equal(model.grand_potential_derivative(10., 3), None)
    assert_equal(model.grand_potential_derivative(10., 4), None)
    assert_equal(model.grand_potential_derivative(9.3, 2), None)
    assert_equal(model.grand_potential_derivative(8.2, 3), None)
    assert_equal(model.grand_potential_derivative(6.1, 4), None)
    assert_equal(model.grand_potential_derivative(9.8, 5), None)
    assert_equal(model.grand_potential_derivative(10.1, 10), None)
    assert_equal(model.grand_potential_derivative(11.5, 2), None)
    assert_equal(model.grand_potential_derivative(12.9, 3), None)
    assert_equal(model.grand_potential_derivative(0.0, 4), None)
    assert_equal(model.grand_potential_derivative(None, 10), None)


# def test_global_linear_np_grand_potential_mu():
#     # E(N) = -1.0 - 0.5 * N, N <= 10
#     # E(N) = -7.0 + 0.1 * N, N >= 10
#     # build global tool
#     model = LinearGlobalTool(-6.0, -5.9, -5.5, 10)
#     # check mu to N conversion


def test_global_linear_np_grand_potential_reactivity():
    # E(N) = -1.0 - 0.5 * N, N <= 10
    # E(N) = -7.0 + 0.1 * N, N >= 10
    # build global tool
    model = LinearGlobalTool({10: -6.0, 11: -5.9, 9: -5.5})
    # check fundamental descriptors
    assert_equal(model.softness, None)
    assert_equal(model.hyper_softness(2), None)
    assert_equal(model.hyper_softness(3), None)
    assert_equal(model.hyper_softness(4), None)
    assert_equal(model.hyper_softness(10), None)


def test_global_linear_nn_energy():
    # E(N) =  0.0 - 0.7 * N, N <= 5
    # E(N) = -1.5 - 0.4 * N, N >= 5
    # build global tool
    model = LinearGlobalTool({5: -3.5, 6: -3.9, 4: -2.8})
    # check parameters
    assert_equal(model.n0, 5.)
    assert_equal(model.n_max, None)
    # check E(N)
    assert_almost_equal(model.energy(5), -3.5, decimal=6)
    assert_almost_equal(model.energy(6), -3.9, decimal=6)
    assert_almost_equal(model.energy(4), -2.8, decimal=6)
    assert_almost_equal(model.energy(0.), 0.0, decimal=6)
    assert_almost_equal(model.energy(3.5), -2.45, decimal=6)
    assert_almost_equal(model.energy(5.7), -3.78, decimal=6)
    # check dE(N)
    assert_equal(model.energy_derivative(5., 1), None)
    assert_almost_equal(model.energy_derivative(4., 1), -0.7, decimal=6)
    assert_almost_equal(model.energy_derivative(6., 1), -0.4, decimal=6)
    assert_almost_equal(model.energy_derivative(4.99, 1), -0.7, decimal=6)
    assert_almost_equal(model.energy_derivative(5.01, 1), -0.4, decimal=6)
    # check d^nE(N) for n > 1
    assert_equal(model.energy_derivative(5.0, 2), None)
    assert_equal(model.energy_derivative(4.5, 2), 0.)
    assert_equal(model.energy_derivative(3.9, 2), 0.)
    assert_equal(model.energy_derivative(5.0, 3), None)
    assert_equal(model.energy_derivative(6.1, 2), 0.)
    assert_equal(model.energy_derivative(7.0, 2), 0.)


def test_global_linear_nn_reactivity():
    # E(N) =  0.0 - 0.7 * N, N <= 5
    # E(N) = -1.5 - 0.4 * N, N >= 5
    # build global tool
    model = LinearGlobalTool({5: -3.5, 6: -3.9, 4: -2.8})
    # check ionization potential and electron affinity
    assert_almost_equal(model.ip, 0.7)
    assert_almost_equal(model.ea, 0.4)
    assert_almost_equal(model.ionization_potential, 0.7)
    assert_almost_equal(model.electron_affinity, 0.4)
    # check fundamental descriptors
    assert_almost_equal(model.mu_minus, -0.7)
    assert_almost_equal(model.mu_plus, -0.4)
    assert_almost_equal(model.mu_zero, -0.55)
    assert_equal(model.mu, None)
    assert_equal(model.chemical_potential, None)
    assert_equal(model.eta, None)
    assert_equal(model.chemical_hardness, None)
    assert_equal(model.hyper_hardness(2), None)
    assert_equal(model.hyper_hardness(3), None)
    assert_equal(model.hyper_hardness(4), None)
    # check derived descriptors
    assert_equal(model.electronegativity, None)
    assert_equal(model.electrofugality, None)
    assert_equal(model.electrophilicity, None)
    assert_equal(model.nucleofugality, None)


def test_global_linear_nn_grand_potential_n():
    # E(N) =  0.0 - 0.7 * N, N <= 5
    # E(N) = -1.5 - 0.4 * N, N >= 5
    # build global tool
    model = LinearGlobalTool({5: -3.5, 6: -3.9, 4: -2.8})
    # check grand potential (as a function of N)
    assert_equal(model.grand_potential(None), None)
    assert_equal(model.grand_potential(5.), None)
    assert_almost_equal(model.grand_potential(5.1), -1.5, decimal=6)
    assert_almost_equal(model.grand_potential(6.0), -1.5, decimal=6)
    assert_almost_equal(model.grand_potential(4.1), 0.0, decimal=6)
    assert_almost_equal(model.grand_potential(3.9), 0.0, decimal=6)
    # check derivative of grand potential
    assert_almost_equal(model.grand_potential_derivative(5.0, 1), -5.0, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(6.0, 1), -6.0, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(4.0, 1), -4.0, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(5.6, 1), -5.6, decimal=6)
    assert_almost_equal(model.grand_potential_derivative(3.3, 1), -3.3, decimal=6)
    assert_equal(model.grand_potential_derivative(5., 2), None)
    assert_equal(model.grand_potential_derivative(4., 3), None)
    assert_equal(model.grand_potential_derivative(5., 4), None)
    assert_equal(model.grand_potential_derivative(5.7, 5), None)
    assert_equal(model.grand_potential_derivative(4.2, 3), None)
    assert_equal(model.grand_potential_derivative(6.3, 2), None)
    assert_equal(model.grand_potential_derivative(3.6, 3), None)
    assert_equal(model.grand_potential_derivative(0.0, 4), None)
    assert_equal(model.grand_potential_derivative(None, 10), None)


# def test_global_linear_nn_grand_potential_mu():
#     # E(N) =  0.0 - 0.7 * N, N <= 5
#     # E(N) = -1.5 - 0.4 * N, N >= 5
#     # build global tool
#     model = LinearGlobalTool(-3.5, -3.9, -2.8, 5)


def test_global_linear_nn_grand_potential_reactivity():
    # E(N) =  0.0 - 0.7 * N, N <= 5
    # E(N) = -1.5 - 0.4 * N, N >= 5
    # build global tool
    model = LinearGlobalTool({5: -3.5, 6: -3.9, 4: -2.8})
    # check fundamental descriptors
    assert_equal(model.softness, None)
    assert_equal(model.hyper_softness(2), None)
    assert_equal(model.hyper_softness(3), None)
    assert_equal(model.hyper_softness(4), None)
    assert_equal(model.hyper_softness(8), None)


def test_local_linear_raises():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # check dictionary
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 11.: dp})
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 11.: dp, 9: dm, 8: d0})
    # check N
    assert_raises(ValueError, LinearLocalTool, {-10.: d0, 11.: dp, 9: dm})
    assert_raises(ValueError, LinearLocalTool, {10.: d0, -11.: dp, -9: dm})
    assert_raises(ValueError, LinearLocalTool, {-10.: d0, -11.: dp, -9: dm})
    assert_raises(ValueError, LinearLocalTool, {0.9: d0, 1.4: dp, 0.4: dm})
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 11.5: dp, 9: dm})
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 12.: dp, 8: dm})
    # check Nmax
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 11.: dp, 9: dm}, -1.5)
    assert_raises(ValueError, LinearLocalTool, {10.: d0, 11.: dp, 9: dm}, -0.2, 1.5)
    # check array size
    assert_raises(ValueError, LinearLocalTool, {2.0: d0, 3.0: dp, 1.0: dm[:-1]})
    assert_raises(ValueError, LinearLocalTool, {2.0: d0, 3.0: dp[:-2], 1.0: dm[:-1]})
    assert_raises(ValueError, LinearLocalTool, {2.0: d0[:-1], 3.0: dp, 1.0: dm[:-1]})
    # build a linear local model
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm})
    # check invalid N
    assert_raises(ValueError, model.density, '10.0')
    assert_raises(ValueError, model.density, -1.)
    assert_raises(ValueError, model.density_derivative, '9.5')
    assert_raises(ValueError, model.density_derivative, -0.1)
    assert_raises(ValueError, model.density_derivative, 2.0, 0.5)
    assert_raises(ValueError, model.density_derivative, 2.0, -1)
    assert_raises(ValueError, model.density_derivative, 2.0, 1.)


def test_local_linear_fake_density():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm})
    # check density
    assert_equal(model.n0, 10.)
    assert_almost_equal(model.density(10.), d0, decimal=6)
    assert_almost_equal(model.density(11.), dp, decimal=6)
    assert_almost_equal(model.density(9.), dm, decimal=6)
    assert_almost_equal(model.density(10.50), 0.5 * d0 + 0.5 * dp, decimal=6)
    assert_almost_equal(model.density(10.20), 0.8 * d0 + 0.2 * dp, decimal=6)
    assert_almost_equal(model.density(10.75), 0.25 * d0 + 0.75 * dp, decimal=6)
    assert_almost_equal(model.density(9.50), 0.5 * dm + 0.5 * d0, decimal=6)
    assert_almost_equal(model.density(9.32), 0.68 * dm + 0.32 * d0, decimal=6)
    assert_almost_equal(model.density(9.61), 0.39 * dm + 0.61 * d0, decimal=6)


def test_local_linear_fake_density_derivative():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm})
    # check first order derivatives
    expected = np.array([-0.5, 1.5, 1.0, -1.0, -2.0])
    assert_almost_equal(model.ff_plus, expected, decimal=6)
    assert_almost_equal(model.density_derivative(10.10, 1), expected, decimal=6)
    assert_almost_equal(model.density_derivative(10.73, 1), expected, decimal=6)
    expected = np.array([0.0, -1.0, 2.0, 0.0, -1.0])
    assert_almost_equal(model.ff_minus, expected, decimal=6)
    assert_almost_equal(model.density_derivative(9.40, 1), expected, decimal=6)
    assert_almost_equal(model.density_derivative(9.95, 1), expected, decimal=6)
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    assert_almost_equal(model.ff_zero, expected, decimal=6)
    assert_almost_equal(model.fukui_function, expected, decimal=6)
    assert_almost_equal(model.density_derivative(10., 1), expected, decimal=6)
    # check higher order derivatives
    assert model.density_derivative(10, 2) is None
    assert model.density_derivative(10., 2) is None
    assert model.density_derivative(10, 3) is None
    assert model.density_derivative(10, 4) is None
    assert_almost_equal(model.density_derivative(9., 2), 0., decimal=6)
    assert_almost_equal(model.density_derivative(10.1, 2), 0., decimal=6)
    assert_almost_equal(model.density_derivative(11., 2), 0., decimal=6)
    assert_almost_equal(model.density_derivative(10.5, 3), 0., decimal=6)
    assert_almost_equal(model.density_derivative(9.86, 4), 0., decimal=6)


def test_local_linear_fake_softness():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model & check softness
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm}, global_softness=0.5)
    expected = 0.5 * np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    assert_almost_equal(model.softness, expected, decimal=6)
    # check softness evaluated at other number of electrons
    expected = 1.5 * np.array([0.0, -1.0, 2.0, 0.0, -1.0])
    assert_almost_equal(1.5 * model.density_derivative(9.34, 1), expected, decimal=6)
    assert_almost_equal(1.5 * model.density_derivative(8.51, 1), expected, decimal=6)
    expected = 0.8 * np.array([-0.5, 1.5, 1.0, -1.0, -2.0])
    assert_almost_equal(0.8 * model.density_derivative(10.42, 1), expected, decimal=6)
    assert_almost_equal(0.8 * model.density_derivative(11.20, 1), expected, decimal=6)
    # build a linear local model & check softness
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm}, global_softness=1.63)
    expected = 1.63 * np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    assert_almost_equal(model.softness, expected, decimal=6)
    # check hyper softness
    assert model.hyper_softness is None
    model = LinearLocalTool({10: d0, 11.: dp, 9: dm}, n_max=10., global_softness=1.63)
    assert model.hyper_softness is None


def test_condensed_linear_raises():
    # fake population array
    pop_0 = np.array([-0.5, 1.0, -0.5, 0.5])
    pop_m = np.array([1.0, -2.0, -0.9, 0.3])
    pop_p = np.array([0.8, -3.1, -0.1, 3.6])
    # check dict_pop
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p, 4.: pop_p})
    assert_raises(ValueError, LinearCondensedTool, {5.: pop_m, 6: pop_0, 7.5: pop_p})
    assert_raises(ValueError, LinearCondensedTool, {-1.: pop_m, 0.: pop_0, 1.: pop_p})
    assert_raises(ValueError, LinearCondensedTool, {0.: pop_m, 1.: pop_0, 3.: pop_p})
    assert_raises(ValueError, LinearCondensedTool, {0.: pop_m, 0.5: pop_0, 1.: pop_p})
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p[:-1]})
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m[:-1], 2.: pop_0, 3.: pop_p[:-1]})
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m, 2.: pop_0[:-1], 3.: pop_p[:-1]})
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m[:-1], 2.: pop_0, 3.: pop_p})
    # check N_max
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p}, -0.5, 1.5)
    assert_raises(ValueError, LinearCondensedTool, {1.: pop_m, 2.: pop_0, 3.: pop_p}, -2., 1.5)
    # make a model
    model = LinearCondensedTool({1.: pop_m, 2.: pop_0, 3.: pop_p})
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
    # build model
    model = LinearCondensedTool({9: pop_09, 10: pop_10, 11: pop_11}, None, None)
    # check N0, Nmax & S
    assert_equal(model.n_ref, 10)
    assert model.n_max is None
    assert model.global_softness is None
    # check population
    assert_almost_equal(model.population(8.5), 1.5 * pop_09 - 0.5 * pop_10, decimal=8)
    assert_almost_equal(model.population(9), pop_09, decimal=8)
    assert_almost_equal(model.population(9.2), 0.8 * pop_09 + 0.2 * pop_10, decimal=8)
    assert_almost_equal(model.population(10), pop_10, decimal=8)
    assert_almost_equal(model.population(10.7), 0.3 * pop_10 + 0.7 * pop_11, decimal=8)
    assert_almost_equal(model.population(11), pop_11, decimal=8)
    assert_almost_equal(model.population(11.3), 1.3 * pop_11 - 0.3 * pop_10, decimal=8)
    # check population derivative
    assert_almost_equal(model.population_derivative(7.1), pop_10 - pop_09, decimal=8)
    assert_almost_equal(model.population_derivative(9.), pop_10 - pop_09, decimal=8)
    assert_almost_equal(model.population_derivative(9.6), pop_10 - pop_09, decimal=8)
    assert_almost_equal(model.population_derivative(10.), 0.5 * (pop_11 - pop_09), decimal=8)
    assert_almost_equal(model.population_derivative(10.9), pop_11 - pop_10, decimal=8)
    assert_almost_equal(model.population_derivative(11), pop_11 - pop_10, decimal=8)
    # check higher order derivatives
    assert model.population_derivative(10, 2) is None
    assert model.population_derivative(10, 3) is None
    assert_almost_equal(model.population_derivative(9., 2), 0., decimal=8)
    assert_almost_equal(model.population_derivative(11., 2), 0., decimal=8)
    assert_almost_equal(model.population_derivative(10.5, 3), 0., decimal=8)
    assert_almost_equal(model.population_derivative(9.75, 4), 0., decimal=8)


def test_condensed_linear_h2o_reactivity():
    # ESP charges of H2O at ub3lyp/ccpvtz
    pop_09 = np.array([4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    pop_10 = np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    pop_11 = np.array([-5.81613550E-01, -2.09193820E-01, -2.09192630E-01])
    # build model
    model = LinearCondensedTool({9: pop_09, 10: pop_10, 11: pop_11}, 10., 2.3)
    # check N0, Nmax & S
    assert_equal(model.n_ref, 10)
    assert_equal(model.n_max, 10)
    assert_equal(model.global_softness, 2.3)
    # check fukui functions
    assert_almost_equal(model.ff_minus, pop_10 - pop_09, decimal=8)
    assert_almost_equal(model.ff_zero, 0.5 * (pop_11 - pop_09), decimal=8)
    assert_almost_equal(model.ff_plus, pop_11 - pop_10, decimal=8)
    assert_almost_equal(model.fukui_function, 0.5 * (pop_11 - pop_09), decimal=8)
    # check other indicators
    assert_almost_equal(model.softness, 2.3 * 0.5 * (pop_11 - pop_09), decimal=8)
    assert model.hyper_softness is None
    assert model.dual_descriptor is None
