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
from chemtools.conceptual.rational import RationalGlobalTool


def test_global_rational1():
    # E(N) = (0.5 - 2.2 N) / (1 + 0.7 N)
    model = RationalGlobalTool(-1.6250, -1.96774193, -1.0, 2.0)
    # check parameters
    np.testing.assert_almost_equal(model.n0, 2.0, decimal=6)
    np.testing.assert_almost_equal(model.params[0], 0.5, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -2.2, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 0.7, decimal=6)
    # check energy values (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy(0), 0.5, decimal=6)
    np.testing.assert_almost_equal(model.energy(1), -1.0, decimal=6)
    np.testing.assert_almost_equal(model.energy(2), -1.6250, decimal=6)
    np.testing.assert_almost_equal(model.energy(3), -1.96774193, decimal=6)
    np.testing.assert_almost_equal(model.energy(4), -2.18421052, decimal=6)
    np.testing.assert_almost_equal(model.energy(5), -2.33333333, decimal=6)
    np.testing.assert_almost_equal(model.energy(6), -2.44230769, decimal=6)
    np.testing.assert_almost_equal(model.energy(1.5), -1.36585365, decimal=6)
    np.testing.assert_almost_equal(model.energy(0.8), -0.80769230, decimal=6)
    # check energy derivatives (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy_derivative(0), -2.55, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1), -0.88235294, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(2), -0.44270833, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3), -0.26534859, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4), -0.17659279, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(5), -0.12592592, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(6), -0.09430473, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.5), -0.60678167, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.8), -1.04783037, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.5, 2), 0.41438748, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.8, 2), 0.94036059, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(1.1, 3), -0.7638260, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(2.5, 4), 0.13346951, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(0.65, 5), -7.74347011, decimal=5)
    np.testing.assert_almost_equal(model.energy_derivative(1.90, 6), 0.827697092, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(3.20, 3), -0.06803109, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4.05, 7), -0.03231737, decimal=6)
    # check global descriptors (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.ip, 0.625, decimal=6)
    np.testing.assert_almost_equal(model.ea, 0.34274193, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.44270833, decimal=6)
    np.testing.assert_almost_equal(model.eta, 0.258246527, decimal=6)
    np.testing.assert_almost_equal(model.ionization_potential, 0.625, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, 0.34274193, decimal=6)
    np.testing.assert_almost_equal(model.chemical_potential, -0.44270833, decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, 0.258246527, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.44270833, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), -0.22596571, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.263626663, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), -0.38445555, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(5), 0.672797214, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(6), -1.37362764, decimal=6)
    np.testing.assert_almost_equal(model.softness, 3.87226890, decimal=6)
    # check n_max and related descriptors (expected values are computed symbolically)
    np.testing.assert_equal(model.n_max, float('inf'))
    np.testing.assert_almost_equal(model.energy(model.n_max), -3.14285714, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, 1.51785714, decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, -1.17511520, decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, 2.14285714, decimal=6)
    # check grand potential (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential(1.), -0.11764705, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(2.), -0.73958333, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(3.), -1.17169614, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(2.78), -1.08950656, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(5.2), -1.74186236, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(0.), 0.5, decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential(model.n_max), , decimal=6)
    # check grand potential derivative (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.), -2.0, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(1.4, 1), -1.4, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.86), -2.86, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2., 2), -3.87226890, decimal=6)
    # expected values based on derived formulas
    n0, a0, a1, b1 = 2.0, 0.5, -2.2, 0.7
    dE = lambda n, r: (-b1)**(r-1) * (a1 - a0*b1) * math.factorial(r) / (1 + b1*n)**(r+1)
    d2omega = lambda n: -1. / dE(n, 2)
    d3omega = lambda n: dE(n, 3) / dE(n, 2)**3
    # d4omega = lambda n: dE(n, 4) / dE(n, 2)**4 - (3 * dE(n, 3)**2) / dE(n, 2)**5
    # d5omega = lambda n: (dE(n, 5)/dE(n, 2)**5 - (10 * dE(n, 3) * dE(n, 4))/dE(n, 2)**6 +
    #                      (15 * dE(n, 3)**3)/dE(n, 2)**7)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.1, 2), d2omega(4.1),
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(3.5, 2), d2omega(3.5),
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.9, 3), d3omega(2.9),
                                   decimal=5)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.67, 1), -4.67, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.67, 2), d2omega(4.67),
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(4.67, 3), d3omega(4.67),
                                   decimal=4)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(1.6, 4), d4omega(1.6),
    #                                decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(2.92, 4), d4omega(2.92),
    #                                decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(5.01, 5), d5omega(5.01),
    #                                decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(4.101, 5), d5omega(4.101),
    #                                decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(model.n_max, 4), ,
    #                                decimal=6)
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
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(5.81, 1)), -5.81,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(4.67, 1), 2),
                                   d2omega(4.67), decimal=5)
    # np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(6.45, 1), 3),
    #                                d3omega(6.45), decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(5.12, 1), 4),
    #                                d4omega(5.12), decimal=6)
    # check hyper-softnesses
    expected = 3.0 * (1 + b1 * n0)**5 / (4 * b1 * (a1 - a0 * b1)**2)
    np.testing.assert_almost_equal(model.hyper_softness(2), expected, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.0, 3), -expected, decimal=6)
    expected = -15 * (1 + b1 * n0)**7 / (8 * b1 * (a1 - a0 * b1)**3)
    np.testing.assert_almost_equal(model.hyper_softness(3), expected, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(2.0, 4), -expected, decimal=6)


def test_global_rational2():
    # E(N) = (-0.15 - 4.2 N) / (1 + 0.45 N)
    model = RationalGlobalTool(-6.99363057, -7.23428571, -6.69064748, 6.5)
    # check parameters
    np.testing.assert_almost_equal(model.n0, 6.5, decimal=6)
    np.testing.assert_almost_equal(model.params[0], -0.15, decimal=6)
    np.testing.assert_almost_equal(model.params[1], -4.2, decimal=6)
    np.testing.assert_almost_equal(model.params[2], 0.45, decimal=6)
    # check energy values (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy(6.5), -6.99363057, decimal=6)
    np.testing.assert_almost_equal(model.energy(7.5), -7.23428571, decimal=6)
    np.testing.assert_almost_equal(model.energy(5.5), -6.69064748, decimal=6)
    np.testing.assert_almost_equal(model.energy(5.0), -6.507692307, decimal=6)
    np.testing.assert_almost_equal(model.energy(8.0), -7.336956521, decimal=6)
    # check energy derivatives (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.energy_derivative(6.5), -0.26824617, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(7.5), -0.21590204, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(5.5), -0.34221831, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(4, 2), 0.16942647, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(10., 3), -0.00548704, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(9.5, 4), 0.002212836, decimal=6)
    # check global descriptors (expected values are computed symbolically)
    np.testing.assert_almost_equal(model.ip, 0.30298309, decimal=6)
    np.testing.assert_almost_equal(model.ea, 0.24065514, decimal=6)
    np.testing.assert_almost_equal(model.mu, -0.26824617, decimal=6)
    np.testing.assert_almost_equal(model.eta, 0.06150867, decimal=6)
    np.testing.assert_almost_equal(model.ionization_potential, 0.30298309, decimal=6)
    np.testing.assert_almost_equal(model.electron_affinity, 0.24065514, decimal=6)
    np.testing.assert_almost_equal(model.chemical_potential, -0.26824617, decimal=6)
    np.testing.assert_almost_equal(model.chemical_hardness, 0.06150867, decimal=6)
    np.testing.assert_almost_equal(model.electronegativity, 0.26824617, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(2), -0.0211558, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(3), 0.00970204, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(4), -0.0055616, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(5), 0.00382587, decimal=6)
    np.testing.assert_almost_equal(model.hyper_hardness(6), -0.0030704, decimal=6)
    np.testing.assert_almost_equal(model.softness, 16.25786868, decimal=6)
    # check n_max and related descriptors (expected values are computed symbolically)
    np.testing.assert_equal(model.n_max, float('inf'))
    np.testing.assert_almost_equal(model.energy(model.n_max), -9.33333333, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 2), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.energy_derivative(model.n_max, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(model.electrophilicity, 2.33970276, decimal=6)
    np.testing.assert_almost_equal(model.nucleofugality, -2.099047619, decimal=6)
    np.testing.assert_almost_equal(model.electrofugality, 2.64268585, decimal=6)
    # check grand potential (as a function of N)
    np.testing.assert_almost_equal(model.grand_potential(6.5), -5.2500304, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(7.91), -5.7468530, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential(0.), -0.15, decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential(model.n_max), , decimal=6)
    # check grand potential derivative (as a function of N)
    # expected values based on derived formulas
    n0, a0, a1, b1 = 6.5, -0.15, -4.2, 0.45
    dE = lambda n, r: (-b1)**(r-1) * (a1 - a0*b1) * math.factorial(r) / (1 + b1*n)**(r+1)
    d2omega = lambda n: -1. / dE(n, 2)
    d3omega = lambda n: dE(n, 3) / dE(n, 2)**3
    # d4omega = lambda n: dE(n, 4) / dE(n, 2)**4 - (3 * dE(n, 3)**2) / dE(n, 2)**5
    # d5omega = lambda n: (dE(n, 5)/dE(n, 2)**5 - (10 * dE(n, 3) * dE(n, 4))/dE(n, 2)**6 +
    #                      (15 * dE(n, 3)**3)/dE(n, 2)**7)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5), -6.5, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(7.1, 1), -7.1, decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(5.8, 2), d2omega(5.8),
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_derivative(0.0, 3), d3omega(0.0),
                                   decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(8.01, 4), d4omega(8.01),
    #                                decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_derivative(6.901, 5), d5omega(6.901),
    #                                decimal=6)
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
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(6.301, 1)), -6.301,
                                   decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(5.55, 1), 2),
                                   d2omega(5.55), decimal=6)
    np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(6.99, 1), 3),
                                   d3omega(6.99), decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(7.1, 1), 4),
    #                                d4omega(7.1), decimal=6)
    # np.testing.assert_almost_equal(model.grand_potential_mu_derivative(dE(7.6, 1), 5),
    #                                d5omega(7.6), decimal=6)
    # check hyper-softnesses
    expected = 3.0 * (1 + b1 * n0)**5 / (4 * b1 * (a1 - a0 * b1)**2)
    np.testing.assert_almost_equal(model.hyper_softness(2), expected, decimal=5)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5, 3), -expected, decimal=5)
    expected = -15 * (1 + b1 * n0)**7 / (8 * b1 * (a1 - a0 * b1)**3)
    np.testing.assert_almost_equal(model.hyper_softness(3), expected, decimal=4)
    np.testing.assert_almost_equal(model.grand_potential_derivative(6.5, 4), -expected, decimal=4)
