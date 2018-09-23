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
"""Test chemtools.conceptual.squareroot module."""

import numpy as np
import sympy as sp
import scipy

from numpy.testing import TestCase, assert_almost_equal, assert_raises, assert_equal
from chemtools.conceptual.squareroot import SquareRootGlobalTool


def make_symbolic_square_root_model(energy_vals):
    n0 = 5
    dict_energy = {n0 - 1.: energy_vals[0], n0: energy_vals[1], n0 + 1.: energy_vals[2]}
    a0, a1, a2, n = sp.symbols('a0 a1 a2 n')

    a1_function = dict_energy[n0 + 1] - 2. * dict_energy[n0] + dict_energy[n0 - 1]
    a1_function /= (sp.sqrt(n0 + 1.) - 2. * sp.sqrt(n0) + sp.sqrt(n0 - 1.))

    a2_function = (dict_energy[n0 + 1] - dict_energy[n0 - 1]) / 2.
    a2_function -= (a1_function * (sp.sqrt(n0 + 1.) - sp.sqrt(n0 - 1.)) / 2.)

    a0_function = dict_energy[n0] - a1_function * sp.sqrt(n0) - a2_function * n0
    energy_function = a0_function + a1_function * sp.sqrt(n) + a2_function * n

    parameters = [a0_function, a1_function, a2_function]
    # Take the derivatives.
    first_deriv = sp.diff(energy_function, n)
    sec_deriv = sp.diff(first_deriv, n)
    third_deriv = sp.diff(sec_deriv, n)
    fourth_deriv = sp.diff(third_deriv, n)

    expr = (energy_function, first_deriv, sec_deriv, third_deriv, fourth_deriv)
    return dict_energy, expr, parameters, n0


def test_parameters():
    # Test parameters for the square root model.
    for energy_minus in np.arange(-1., 1000., 200):
        for energy0 in np.arange(energy_minus - 2, 1000., 200):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)
                assert_almost_equal(sqrt_root.params[0], params[0], decimal=5)
                assert_almost_equal(sqrt_root.params[1], params[1], decimal=5)
                assert_almost_equal(sqrt_root.params[2], params[2], decimal=5)


def test_energy():
    # Test energy values for the square root model.
    for energy_minus in np.arange(-1., 1000., 300):
        for energy0 in np.arange(energy_minus - 2, 1000., 300):
            for energy1 in np.arange(energy0 + 1, 1000., 310):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)
                assert_almost_equal(sqrt_root.energy(4), expr[0].subs('n', 4).evalf(),
                                    decimal=5)
                assert_almost_equal(sqrt_root.energy(5), expr[0].subs('n', 5).evalf(),
                                    decimal=5)
                assert_almost_equal(sqrt_root.energy(4.5),
                                    expr[0].subs('n', 4.5).evalf(), decimal=5)
                assert_raises(ValueError, sqrt_root.energy, -5)

                # Test at infinity
                actual = sqrt_root.energy(np.inf)
                if np.isinf(actual):
                    if actual > 0.:
                        actual = sp.oo
                    else:
                        actual = -sp.oo
                assert_equal(actual, sp.limit(expr[0], "n", sp.oo))


def test_energy_derivative():
    # Test derivative of energy for the square root model.
    for energy_minus in np.arange(-1., 1000., 200):
        for energy0 in np.arange(energy_minus - 2, 1000., 200):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)

                # Go through points ranging from 4 to 8.
                for n in range(4, 8):
                    assert_almost_equal(sqrt_root.energy_derivative(n, 1),
                                        expr[1].subs("n", n).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n, 2),
                                        expr[2].subs("n", n).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n, 3),
                                        expr[3].subs("n", n).evalf())
                    assert_almost_equal(sqrt_root.energy_derivative(n, 4),
                                        expr[4].subs("n", n).evalf())


def test_chemical_concepts():
    # Test chemical concepts for the square root model.
    for energy_minus in np.arange(-1., 1000., 200):
        for energy0 in np.arange(energy_minus - 2, 1000., 200):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)

                sqrt_root = SquareRootGlobalTool(energy)

                # Test Electronegativity
                assert_almost_equal(sqrt_root.electronegativity, -expr[1].subs("n", n0))

                # Test Chemical Potential
                assert_almost_equal(sqrt_root.chemical_potential, expr[1].subs("n", n0))

                # Test Chemical Hardness
                assert_almost_equal(sqrt_root.chemical_hardness, expr[2].subs("n", n0))

                # Test Hyper-Hardness
                assert_almost_equal(sqrt_root.hyper_hardness(2), expr[3].subs("n", n0))
                assert_almost_equal(sqrt_root.hyper_hardness(3), expr[4].subs("n", n0))
                assert_almost_equal(sqrt_root.hyper_hardness(4),
                                    sp.diff(expr[4], "n").subs("n", n0))

                # Test Softness
                assert_almost_equal(sqrt_root.softness, 1. / expr[2].subs("n", n0))

                # Test Hyper Softness
                desired = -expr[3].subs("n", n0) / expr[2].subs("n", n0) ** 3.
                assert_almost_equal(sqrt_root.hyper_softness(2), desired.evalf(), decimal=5)

                # Test Grand Potential
                n = sp.symbols("n")
                grand_function = expr[0] - expr[1] * n
                assert_almost_equal(sqrt_root.grand_potential(6), grand_function.subs(n, 6).evalf())
                assert_almost_equal(sqrt_root.grand_potential(5), grand_function.subs(n, 5).evalf())

                # Test Nucleofugality
                nucleofugality = expr[0].subs(n, n + 1)
                sign = 1
                if np.sign(n0 + 1 - sqrt_root.n_max) == -1:
                    sign = -1
                # Test Nucleofugality at infinite values.
                if sqrt_root.n_max == np.inf:
                    limit = sp.limit(nucleofugality, "n", sp.oo)
                    if limit == sp.oo:
                        desired = sign * (-np.inf)
                    elif limit == -sp.oo:
                        desired = sign * (np.inf)
                    else:
                        desired = sign * (nucleofugality - limit)
                else:
                    nucleofugality -= expr[0].subs(n, sqrt_root.n_max)
                    desired = nucleofugality.subs(n, n0).evalf()
                assert_almost_equal(sqrt_root.nucleofugality, desired)

                # Test Electrofugality
                electrofugality_f = expr[0].subs(n, n - 1)
                electrofugality_f -= expr[0].subs(n, sqrt_root.n_max)
                sign = 1
                if np.sign(sqrt_root.n_max - n0 + 1.) == -1.:
                    sign = -1
                desired = sign * electrofugality_f.subs(n, n0).evalf()
                if sqrt_root.n_max == np.inf:
                    desired = sign * np.inf
                assert_almost_equal(sqrt_root.electrofugality, desired)

                # Test Electrophilicity
                n = sp.symbols('n')
                assert_almost_equal(sqrt_root.electronegativity, (-1) * expr[1].subs(n, n0))

                # Test ionization
                ionization = expr[0].subs(n, n0 - 1) - expr[0].subs(n, n0)
                assert_almost_equal(sqrt_root.ionization_potential, ionization)
                assert_almost_equal(sqrt_root.ip, ionization)

                # Test Electron affinity
                ea = expr[0].subs(n, n0) - expr[0].subs(n,  n0 + 1)
                assert_almost_equal(sqrt_root.electron_affinity, ea)
                assert_almost_equal(sqrt_root.ea, ea)


def test_nmax_using_scipy():
    # Test nmax using scipy.optimize for the square root model.
    for energy_minus in np.arange(-1., 1000., 102):
        for energy0 in np.arange(energy_minus - 2, 1000., 101):
            for energy1 in np.arange(energy0 + 1, 1000., 210):
                energy_vals = [energy_minus, energy0, energy1]
                energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)
                sqrt = SquareRootGlobalTool(energy)

                # Test using scipy minimize.
                energy_minimize = sp.utilities.lambdify(sp.symbols('n'), expr[0])
                scipy_solution = scipy.optimize.minimize_scalar(energy_minimize, 5.,
                                                                bounds=(0, 100.),
                                                                method='bounded')

                # Minima occurs at the upper bound. indicating the energy going to -infinity.
                if np.abs(scipy_solution.x - 100.) < 1e-4:
                    scipy_solution.x = np.inf
                else:
                    # Test first derivative is zero at minima.
                    assert_almost_equal(expr[1].subs('n', sqrt.n_max).evalf(), 0)
                assert_almost_equal(sqrt.n_max, scipy_solution.x, decimal=5)


def test_nmax_using_fixed_examples():
    # Test nmax using fixed examples for the square root model.
    energy_vals = [1., 0.5, 3.]
    energy, expr, params, n0 = make_symbolic_square_root_model(energy_vals)
    sqrt = SquareRootGlobalTool(energy)

    sqrt._params = [1., 1., 1.]
    desired = 0.
    assert_equal(sqrt._compute_nmax(), desired)

    sqrt._params = [1, -1., 1.]
    desired = 0.25
    assert_equal(sqrt._compute_nmax(), desired)

    for x in [0., -1.]:
        sqrt._params = [1., -1., x]
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
