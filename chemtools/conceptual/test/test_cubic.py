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


import numpy as np
import sympy as sp

from numpy.testing import TestCase

from chemtools.conceptual.cubic import CubicGlobalTool


class TestCubicRootModel(TestCase):
    """Cubic Objects are instantiated with different omega values."""

    def setUp(self):
        # Set Up Different Cubic Objects
        self.energy_zero1 = 100.
        self.energy_minus1 = 25.3
        self.energy_plus1 = 50.5
        dict_energy = {9: self.energy_minus1, 10: self.energy_zero1, 11: self.energy_plus1}
        self.n0 = 10
        self.omega1 = 1. / 2.
        self.omega2 = 1. / 3.
        self.omega3 = 1.
        self.cubic_half = CubicGlobalTool(dict_energy, self.omega1)
        self.cubic_third = CubicGlobalTool(dict_energy, self.omega2)
        self.cubic_one = CubicGlobalTool(dict_energy, self.omega3)

        # Set Up Parameters for substitution
        self.parameters_half = [('energy_plus', self.energy_plus1),
                                ('energy_minus', self.energy_minus1),
                                ('energy_zero', self.energy_zero1),
                                ('omega', self.omega1)]
        self.parameters_third = [('energy_plus', self.energy_plus1),
                                 ('energy_minus', self.energy_minus1),
                                 ('energy_zero', self.energy_zero1),
                                 ('omega', self.omega2)]
        self.parameters_one = [('energy_plus', self.energy_plus1),
                               ('energy_minus', self.energy_minus1),
                               ('energy_zero', self.energy_zero1),
                               ('omega', self.omega3)]
        # Set Up Sympy Symbolic Functions for Testing
        N = sp.symbols('N')
        energy_plus, energy_zero, energy_minus, omega = sp.symbols(
            'energy_plus, energy_zero, energy_minus, omega')
        self.I = energy_minus - energy_zero
        self.A = energy_zero - energy_plus

        self.a0 = energy_zero
        self.a1 = - (omega * self.I + (1 - omega) * self.A)
        self.a2 = (self.I - self.A) / 2.
        self.a3 = (2. * omega - 1) * (self.I - self.A) / 2.

        # Set Up Energy Symbolic Function for Testing
        self.energy_function = self.a0 + self.a1 * N + self.a2 * N ** 2 + self.a3 * N ** 3
        self.first_deriv = sp.diff(self.energy_function, N)
        self.sec_deriv = sp.diff(self.first_deriv, N)
        self.third_deriv = sp.diff(self.sec_deriv, N)

    def test_omega(self):
        np.testing.assert_equal(self.cubic_half.omega, .5)
        np.testing.assert_equal(self.cubic_third.omega, 1./3.)
        np.testing.assert_equal(self.cubic_one.omega, 1.)

    def test_parameters_with_omega_half(self):
        a0_sympy = self.a0.subs(self.parameters_half).evalf()
        a1_sympy_half = self.a1.subs(self.parameters_half).evalf()
        a2_sympy_half = self.a2.subs(self.parameters_half).evalf()
        a3_sympy_half = self.a3.subs(self.parameters_half).evalf()
        np.testing.assert_array_almost_equal(self.cubic_half.params,
                                             [a0_sympy, a1_sympy_half, a2_sympy_half,
                                              a3_sympy_half],
                                             decimal=5)

    def test_parameters_with_omega_third(self):
        a0_sympy = self.a0.subs(self.parameters_third).evalf()
        a1_sympy_third = self.a1.subs(self.parameters_third).evalf()
        a2_sympy_third = self.a2.subs(self.parameters_third).evalf()
        a3_sympy_third = self.a3.subs(self.parameters_third).evalf()
        np.testing.assert_array_almost_equal(self.cubic_third.params,
                                             [a0_sympy, a1_sympy_third, a2_sympy_third,
                                              a3_sympy_third],
                                             decimal=5)

    def test_parameters_with_omega_one(self):
        a0_sympy = self.a0.subs(self.parameters_one).evalf()
        a1_sympy_one = self.a1.subs(self.parameters_one).evalf()
        a2_sympy_one = self.a2.subs(self.parameters_one).evalf()
        a3_sympy_one = self.a3.subs(self.parameters_one).evalf()
        np.testing.assert_array_almost_equal(self.cubic_one.params,
                                             [a0_sympy, a1_sympy_one, a2_sympy_one, a3_sympy_one],
                                             decimal=5)

    def test_energy_model_with_omega_half(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        energy_function_half = self.energy_function.subs(self.parameters_half)
        energy_values = [energy_function_half.subs('N', n).evalf() for n in n_values]
        actual_energy_values = [self.cubic_half.energy(n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_energy_values, energy_values)

    def test_energy_model_with_omega_third(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        energy_function_third = self.energy_function.subs(self.parameters_third)
        energy_values = [energy_function_third.subs('N', n).evalf() for n in n_values]
        actual_energy_values = [self.cubic_third.energy(n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_energy_values, energy_values)

    def test_energy_model_with_omega_one(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        energy_function_one = self.energy_function.subs(self.parameters_one)
        energy_values = [energy_function_one.subs('N', n).evalf() for n in n_values]
        actual_energy_values = [self.cubic_one.energy(n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_energy_values, energy_values)

    def test_energy_first_derivative_with_omega_half(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        first_deriv_half = self.first_deriv.subs(self.parameters_half)
        actual_values = [self.cubic_half.energy_derivative(n) for n in n_values]
        desired_values = [first_deriv_half.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_first_derivative_with_omega_third(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        first_deriv_third = self.first_deriv.subs(self.parameters_third)
        actual_values = [self.cubic_third.energy_derivative(n) for n in n_values]
        desired_values = [first_deriv_third.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_first_derivative_with_omega_one(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        first_deriv_one = self.first_deriv.subs(self.parameters_one)
        actual_values = [self.cubic_one.energy_derivative(n) for n in n_values]
        desired_values = [first_deriv_one.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_sec_derivative_with_omega_half(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        sec_deriv_half = self.sec_deriv.subs(self.parameters_half)
        actual_values = [self.cubic_half.energy_derivative(n, order=2) for n in n_values]
        desired_values = [sec_deriv_half.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_sec_derivative_with_omega_third(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        sec_deriv_third = self.sec_deriv.subs(self.parameters_third)
        actual_values = [self.cubic_third.energy_derivative(n, order=2) for n in n_values]
        desired_values = [sec_deriv_third.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_sec_derivative_with_omega_one(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        sec_deriv_one = self.sec_deriv.subs(self.parameters_one)
        actual_values = [self.cubic_one.energy_derivative(n, order=2) for n in n_values]
        desired_values = [sec_deriv_one.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_third_derivative_with_omega_half(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        sec_deriv_half = self.third_deriv.subs(self.parameters_half)
        actual_values = [self.cubic_half.energy_derivative(n, order=3) for n in n_values]
        desired_values = [sec_deriv_half.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_third_derivative_with_omega_third(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        third_deriv_third = self.third_deriv.subs(self.parameters_third)
        actual_values = [self.cubic_third.energy_derivative(n, order=3) for n in n_values]
        desired_values = [third_deriv_third.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_energy_third_derivative_with_omega_one(self):
        n_values = [1., 1.5, 2., 3.5, 4., 4.5, 5., 8., 9.5, 10., 11.]
        third_deriv_one = self.third_deriv.subs(self.parameters_one)
        actual_values = [self.cubic_one.energy_derivative(n, order=3) for n in n_values]
        desired_values = [third_deriv_one.subs('N', n) for n in n_values]
        np.testing.assert_array_almost_equal(actual_values, desired_values)

    def test_higher_energy_derivative_with_omega_half(self):
        order_terms = [4., 5., 6., 7., 8.]
        actual_energy = [self.cubic_half.energy_derivative(5., order=o) for o in order_terms]
        desired_values = np.zeros(len(actual_energy))
        np.testing.assert_array_equal(actual_energy, desired_values)

    def test_higher_energy_derivative_with_omega_third(self):
        order_terms = [4., 5., 6., 7., 8.]
        actual_energy = [self.cubic_third.energy_derivative(5., order=o) for o in order_terms]
        desired_values = np.zeros(len(actual_energy))
        np.testing.assert_array_equal(actual_energy, desired_values)

    def test_higher_energy_derivative_with_omega_one(self):
        order_terms = [4., 5., 6., 7., 8.]
        actual_energy = [self.cubic_one.energy_derivative(5., order=o) for o in order_terms]
        desired_values = np.zeros(len(actual_energy))
        np.testing.assert_array_equal(actual_energy, desired_values)

    def test_ionization_with_omega_half(self):
        ionization = self.I.subs([('energy_minus', self.energy_minus1),
                                  ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_half.ionization_potential, ionization)

    def test_ionization_with_omega_third(self):
        ionization = self.I.subs([('energy_minus', self.energy_minus1),
                                  ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_third.ionization_potential, ionization)

    def test_ionization_with_omega_one(self):
        ionization = self.I.subs([('energy_minus', self.energy_minus1),
                                  ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_one.ionization_potential, ionization)

    def test_electron_affinity_with_omega_half(self):
        electron_affin = self.A.subs([('energy_plus', self.energy_plus1),
                                      ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_half.electron_affinity, electron_affin)

    def test_electron_affinity_with_omega_third(self):
        electron_affin = self.A.subs([('energy_plus', self.energy_plus1),
                                      ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_third.electron_affinity, electron_affin)

    def test_electron_affinity_with_omega_one(self):
        electron_affin = self.A.subs([('energy_plus', self.energy_plus1),
                                      ('energy_zero', self.energy_zero1)]).evalf()
        np.testing.assert_almost_equal(self.cubic_one.electron_affinity, electron_affin)

    def test_chemical_potential_with_omega_half(self):
        first_deriv_half = self.first_deriv.subs(self.parameters_half)
        np.testing.assert_almost_equal(self.cubic_half.chemical_potential,
                                       first_deriv_half.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_half.mu, first_deriv_half.subs('N', self.n0))

    def test_chemical_potential_with_omega_third(self):
        first_deriv_third = self.first_deriv.subs(self.parameters_third)
        np.testing.assert_almost_equal(self.cubic_third.chemical_potential,
                                       first_deriv_third.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_third.mu, first_deriv_third.subs('N', self.n0))

    def test_chemical_potential_with_omega_one(self):
        first_deriv_one = self.first_deriv.subs(self.parameters_one)
        np.testing.assert_almost_equal(self.cubic_one.chemical_potential,
                                       first_deriv_one.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_one.mu, first_deriv_one.subs('N', self.n0))

    def test_chemical_hardness_with_omega_half(self):
        sec_deriv_half = self.sec_deriv.subs(self.parameters_half)
        np.testing.assert_almost_equal(self.cubic_half.chemical_hardness,
                                       sec_deriv_half.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_half.eta, sec_deriv_half.subs('N', self.n0))

    def test_chemical_hardness_with_omega_third(self):
        sec_deriv_third = self.sec_deriv.subs(self.parameters_third)
        np.testing.assert_almost_equal(self.cubic_third.chemical_hardness,
                                       sec_deriv_third.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_third.eta, sec_deriv_third.subs('N', self.n0))

    def test_chemical_hardness_with_omega_one(self):
        sec_deriv_one = self.sec_deriv.subs(self.parameters_one)
        np.testing.assert_almost_equal(self.cubic_one.chemical_hardness,
                                       sec_deriv_one.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_one.eta, sec_deriv_one.subs('N', self.n0))

    def test_chemical_softness_with_omega_half(self):
        sec_deriv_half = self.sec_deriv.subs(self.parameters_half)
        np.testing.assert_almost_equal(self.cubic_half.softness,
                                       1. / sec_deriv_half.subs('N', self.n0))

    def test_chemical_softness_with_omega_third(self):
        sec_deriv_third = self.sec_deriv.subs(self.parameters_third)
        np.testing.assert_almost_equal(self.cubic_third.softness,
                                       1. / sec_deriv_third.subs('N', self.n0))

    def test_chemical_softness_with_omega_one(self):
        sec_deriv_one = self.sec_deriv.subs(self.parameters_one)
        np.testing.assert_almost_equal(self.cubic_one.softness,
                                       1. / sec_deriv_one.subs('N', self.n0))

    def test_chemical_hyper_softness_with_omega_half(self):
        sec_deriv_half = self.sec_deriv.subs(self.parameters_half)
        third_deriv_half = self.third_deriv.subs(self.parameters_half)
        np.testing.assert_almost_equal(self.cubic_half.hyper_softness(2),
                                       -third_deriv_half.subs('N', self.n0) / sec_deriv_half.subs(
                                           'N', self.n0)**3)

    def test_chemical_hyper_softness_with_omega_third(self):
        sec_deriv_third = self.sec_deriv.subs(self.parameters_third)
        third_deriv_third = self.third_deriv.subs(self.parameters_third)
        np.testing.assert_almost_equal(self.cubic_third.hyper_softness(2),
                                       -third_deriv_third.subs('N', self.n0) / sec_deriv_third.subs(
                                           'N', self.n0)**3)

    def test_chemical_hyper_softness_with_omega_one(self):
        sec_deriv_one = self.sec_deriv.subs(self.parameters_one)
        third_deriv_one = self.third_deriv.subs(self.parameters_one)
        np.testing.assert_almost_equal(self.cubic_one.hyper_softness(2),
                                       -third_deriv_one.subs('N', self.n0) / sec_deriv_one.subs(
                                           'N', self.n0)**3)

    def test_chemical_hyper_hardness_with_omega_half(self):
        third_deriv_half = self.third_deriv.subs(self.parameters_half)
        np.testing.assert_almost_equal(self.cubic_half.hyper_hardness(2),
                                       third_deriv_half.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_half.hyper_hardness(3), 0)
        np.testing.assert_almost_equal(self.cubic_half.hyper_hardness(4), 0)

    def test_chemical_hyper_hardness_with_omega_third(self):
        third_deriv_third = self.third_deriv.subs(self.parameters_third)
        np.testing.assert_almost_equal(self.cubic_third.hyper_hardness(2),
                                       third_deriv_third.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_third.hyper_hardness(3), 0)
        np.testing.assert_almost_equal(self.cubic_third.hyper_hardness(4), 0)

    def test_chemical_hyper_hardness_with_omega_one(self):
        third_deriv_one = self.third_deriv.subs(self.parameters_one)
        np.testing.assert_almost_equal(self.cubic_one.hyper_hardness(2),
                                       third_deriv_one.subs('N', self.n0))
        np.testing.assert_almost_equal(self.cubic_one.hyper_hardness(3), 0)
        np.testing.assert_almost_equal(self.cubic_one.hyper_hardness(4), 0)
