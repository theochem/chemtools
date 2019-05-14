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
# pragma pylint: disable=protected-access
"""Test chemtools.conceptual.squareroot module."""

import numpy as np
from numpy.testing import assert_almost_equal, assert_raises, assert_equal
import sympy as sp
import scipy

from chemtools.conceptual.squareroot import SquareRootGlobalTool


def make_symbolic_square_root_model(energy_vals):
    r"""Create symbolic square root energy model and it's derivatives, used for testing."""
    n_ref = 5
    dict_energy = {n_ref - 1.: energy_vals[0], n_ref: energy_vals[1], n_ref + 1.: energy_vals[2]}
    n_elec = sp.symbols('n_elec')

    a1_function = dict_energy[n_ref + 1] - 2. * dict_energy[n_ref] + dict_energy[n_ref - 1]
    a1_function /= (sp.sqrt(n_ref + 1.) - 2. * sp.sqrt(n_ref) + sp.sqrt(n_ref - 1.))

    a2_function = (dict_energy[n_ref + 1] - dict_energy[n_ref - 1]) / 2.
    a2_function -= (a1_function * (sp.sqrt(n_ref + 1.) - sp.sqrt(n_ref - 1.)) / 2.)

    a0_function = dict_energy[n_ref] - a1_function * sp.sqrt(n_ref) - a2_function * n_ref
    energy_function = a0_function + a1_function * sp.sqrt(n_elec) + a2_function * n_elec

    parameters = [a0_function, a1_function, a2_function]
    # Take the derivatives.
    first_deriv = sp.diff(energy_function, n_elec)
    sec_deriv = sp.diff(first_deriv, n_elec)
    third_deriv = sp.diff(sec_deriv, n_elec)
    fourth_deriv = sp.diff(third_deriv, n_elec)

    expr = (energy_function, first_deriv, sec_deriv, third_deriv, fourth_deriv)
    return dict_energy, expr, parameters, n_ref


def test_square_root_raises():
    # check invalid N0
    assert_raises(ValueError, SquareRootGlobalTool, {0: -15.0, 1: -14.4, -1: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {0.3: -15.0, 1.3: -14.4, -0.7: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {0.98: -15.0, 1.98: -14.4, -0.02: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {-1.: -15.0, 0.: -14.9, -2.: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {-2: -15.0, -1: -14.9, -3: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {-2: -15.0, -1: -14.9, -3: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {0.0: -15.0, 0.5: -14.9, 1.0: -14.0})
    assert_raises(ValueError, SquareRootGlobalTool, {2.: -15.0, 3.5: -14.9, 4.0: -14.0})


def test_parameters():
    # Test parameters for the square root model.
    for energy_minus in np.arange(-1., 1000., 200):
        for energy0 in np.arange(energy_minus - 2, 1000., 200):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, _, params, _ = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)
                assert_almost_equal(sqrt_root.params[0], params[0], decimal=5)
                assert_almost_equal(sqrt_root.params[1], params[1], decimal=5)
                assert_almost_equal(sqrt_root.params[2], params[2], decimal=5)


def test_energy():
    # Test energy values for the square root model.
    for energy_minus in np.arange(-1., 1000., 300):
        for energy0 in np.arange(energy_minus - 2, 1000., 400):
            for energy1 in np.arange(energy0 + 1, 1000., 410):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, _, _ = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)
                assert_almost_equal(sqrt_root.energy(4), expr[0].subs("n_elec", 4).evalf(),
                                    decimal=5)
                assert_almost_equal(sqrt_root.energy(5), expr[0].subs("n_elec", 5).evalf(),
                                    decimal=5)
                assert_almost_equal(sqrt_root.energy(4.5),
                                    expr[0].subs("n_elec", 4.5).evalf(), decimal=5)
                assert_raises(ValueError, sqrt_root.energy, -5)

                # Test at infinity
                actual = sqrt_root.energy(np.inf)
                desired = sp.limit(expr[0], "n_elec", sp.oo)
                if desired == sp.oo:
                    desired = np.inf
                elif desired == -sp.oo:
                    desired = -np.inf
                assert_equal(actual, desired)


def test_energy_derivative():
    # Test derivative of energy for the square root model.
    for energy_minus in np.arange(-1., 1000., 400):
        for energy0 in np.arange(energy_minus - 2, 1000., 400):
            for energy1 in np.arange(energy0 + 1, 1000., 410):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, _, _ = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)

                # Go through points ranging from 4 to 8.
                for n_elec in range(4, 8):
                    assert_almost_equal(sqrt_root.energy_derivative(n_elec, 1),
                                        expr[1].subs("n_elec", n_elec).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n_elec, 2),
                                        expr[2].subs("n_elec", n_elec).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n_elec, 3),
                                        expr[3].subs("n_elec", n_elec).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n_elec, 4),
                                        expr[4].subs("n_elec", n_elec).evalf())

                # Test at infinity
                assert_almost_equal(sqrt_root.energy_derivative(np.inf, 1),
                                    sp.limit(expr[1], "n_elec", sp.oo).evalf())
                assert_almost_equal(sqrt_root.energy_derivative(np.inf, 2),
                                    sp.limit(expr[2], "n_elec", sp.oo).evalf())

    # Test assertions
    assert_raises(ValueError, sqrt_root.energy_derivative, 4, -2)
    assert_raises(ValueError, sqrt_root.energy_derivative, 4, 1.5)


def test_chemical_concepts():
    # Test chemical concepts for the square root model.
    for energy_minus in np.arange(-1., 1000., 300):
        for energy0 in np.arange(energy_minus - 2, 1000., 300):
            for energy1 in np.arange(energy0 + 1, 1000., 310):
                energy, expr, _, n_ref = make_symbolic_square_root_model([energy_minus,
                                                                          energy0,
                                                                          energy1])
                sqrt_root = SquareRootGlobalTool(energy)

                # Test Electronegativity
                desired = -1. * expr[1].subs("n_elec", n_ref)
                assert_almost_equal(sqrt_root.electronegativity, desired)

                # Test Chemical Potential
                assert_almost_equal(sqrt_root.chemical_potential, expr[1].subs("n_elec", n_ref))

                # Test Chemical Hardness
                assert_almost_equal(sqrt_root.chemical_hardness, expr[2].subs("n_elec", n_ref))

                # Test Hyper-Hardness
                assert_almost_equal(sqrt_root.hyper_hardness(2), expr[3].subs("n_elec", n_ref))
                assert_almost_equal(sqrt_root.hyper_hardness(3), expr[4].subs("n_elec", n_ref))
                assert_almost_equal(sqrt_root.hyper_hardness(4),
                                    sp.diff(expr[4], "n_elec").subs("n_elec", n_ref))

                # Test Softness
                assert_almost_equal(sqrt_root.softness, 1. / expr[2].subs("n_elec", n_ref))

                # Test Hyper Softness
                desired = -1. * expr[3].subs("n_elec", n_ref) / expr[2].subs("n_elec", n_ref) ** 3.
                assert_almost_equal(sqrt_root.hyper_softness(2), desired.evalf(), decimal=5)

                # Test Grand Potential
                n_elec = sp.symbols("n_elec")
                grand_function = expr[0] - expr[1] * n_elec
                assert_almost_equal(sqrt_root.grand_potential(6),
                                    grand_function.subs(n_elec, 6).evalf())
                assert_almost_equal(sqrt_root.grand_potential(5),
                                    grand_function.subs(n_elec, 5).evalf())

                # Test Nucleofugality
                nucleofugality = expr[0].subs(n_elec, n_elec + 1)
                sign = 1
                if np.sign(n_ref + 1 - sqrt_root.n_max) == -1:
                    sign = -1
                # Test Nucleofugality at infinite values.
                if sqrt_root.n_max == np.inf:
                    limit = sp.limit(nucleofugality, "n_elec", sp.oo)
                    if limit == -sp.oo:
                        desired = sign * np.inf
                else:
                    nucleofugality -= expr[0].subs(n_elec, sqrt_root.n_max)
                    desired = nucleofugality.subs(n_elec, n_ref).evalf()
                assert_almost_equal(sqrt_root.nucleofugality, desired)

                # Test Electrofugality
                electrofugality_f = expr[0].subs(n_elec, n_elec - 1)
                electrofugality_f -= expr[0].subs(n_elec, sqrt_root.n_max)
                sign = 1
                if np.sign(sqrt_root.n_max - n_ref + 1.) == -1.:
                    sign = -1
                desired = sign * electrofugality_f.subs(n_elec, n_ref).evalf()
                if sqrt_root.n_max == np.inf:
                    desired = sign * np.inf
                assert_almost_equal(sqrt_root.electrofugality, desired)

                # Test Electrophilicity
                n_elec = sp.symbols("n_elec")
                assert_almost_equal(sqrt_root.electronegativity, (-1) * expr[1].subs(n_elec, n_ref))

                # Test ionization
                ionization = expr[0].subs(n_elec, n_ref - 1) - expr[0].subs(n_elec, n_ref)
                assert_almost_equal(sqrt_root.ionization_potential, ionization)
                assert_almost_equal(sqrt_root.ip, ionization)

                # Test Electron affinity
                electron_aff = expr[0].subs(n_elec, n_ref) - expr[0].subs(n_elec, n_ref + 1)
                assert_almost_equal(sqrt_root.electron_affinity, electron_aff)
                assert_almost_equal(sqrt_root.ea, electron_aff)


def test_nmax_using_scipy():
    # Test nmax using scipy.optimize for the square root model.
    for energy_minus in np.arange(-1., 1000., 102):
        for energy0 in np.arange(energy_minus - 2, 1000., 101):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, _, _ = make_symbolic_square_root_model(energy_vals)
                sqrt = SquareRootGlobalTool(energy)

                # Test using scipy minimize.
                energy_minimize = sp.utilities.lambdify(sp.symbols("n_elec"), expr[0])
                scipy_solution = scipy.optimize.minimize_scalar(energy_minimize, 5.,
                                                                bounds=(0, 100.),
                                                                method='bounded')

                # Minima occurs at the upper bound. indicating the energy going to -infinity.
                if np.abs(scipy_solution.x - 100.) < 1e-4:
                    scipy_solution.x = np.inf
                else:
                    # Test first derivative is zero at minima.
                    assert_almost_equal(expr[1].subs("n_elec", sqrt.n_max).evalf(), 0)
                assert_almost_equal(sqrt.n_max, scipy_solution.x, decimal=5)


def test_nmax_using_fixed_examples():
    # Test nmax using fixed examples for the square root model.
    energy_vals = [1., 0.5, 3.]
    energy, _, _, _ = make_symbolic_square_root_model(energy_vals)
    sqrt = SquareRootGlobalTool(energy)

    sqrt._params = [1., 1., 1.]
    desired = 0.
    assert_equal(sqrt._compute_nmax(), desired)

    sqrt._params = [1, -1., 1.]
    desired = 0.25
    assert_equal(sqrt._compute_nmax(), desired)

    for case in [0., -1.]:
        sqrt._params = [1., -1., case]
        desired = np.inf
        assert_equal(sqrt._compute_nmax(), desired)

    sqrt._params = [1., 1., -1]
    desired = np.inf
    assert_equal(sqrt._compute_nmax(), desired)

    sqrt._params = [1., 0., 5.]
    assert_raises(ValueError, sqrt._compute_nmax)

    sqrt._params = [1., 1., 0.]
    desired = 0.
    assert_equal(sqrt._compute_nmax(), desired)
