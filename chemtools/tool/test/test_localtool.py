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
#pylint: skip-file


from chemtools import *


def test_local_linear():
    # Fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # Build a linear local model
    model = LinearLocalTool(d0, dp, dm)
    # Check density
    np.testing.assert_almost_equal(model.density_zero, d0, decimal=8)
    np.testing.assert_almost_equal(model.density_plus, dp, decimal=8)
    np.testing.assert_almost_equal(model.density_minus, dm, decimal=8)
    # Check descriptos
    expected = np.array([-0.5, 1.5, 1.0, -1.0, -2.0])
    np.testing.assert_almost_equal(model.ff_plus, expected, decimal=8)
    expected = np.array([0.0, -1.0, 2.0, 0.0, -1.0])
    np.testing.assert_almost_equal(model.ff_minus, expected, decimal=8)
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    np.testing.assert_almost_equal(model.ff_zero, expected, decimal=8)
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0])
    np.testing.assert_almost_equal(model.dual_descriptor, expected, decimal=8)


def test_local_quadratic():
    # Fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # Build a linear local model
    model = QuadraticLocalTool(d0, dp, dm)
    # Check density
    np.testing.assert_almost_equal(model.density_zero, d0, decimal=8)
    np.testing.assert_almost_equal(model.density_plus, dp, decimal=8)
    np.testing.assert_almost_equal(model.density_minus, dm, decimal=8)
    # Check descriptos
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    np.testing.assert_almost_equal(model.fukui_function, expected, decimal=8)
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0])
    np.testing.assert_almost_equal(model.dual_descriptor, expected, decimal=8)
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 3.5
    np.testing.assert_almost_equal(model.softness(3.5), expected, decimal=8)
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5]) / 0.5
    np.testing.assert_almost_equal(model.softness(0.5), expected, decimal=8)
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (1.5 * 1.5)
    np.testing.assert_almost_equal(model.hyper_softness(1.5), expected, decimal=8)
    expected = np.array([-0.5, 2.5, -1.0, -1.0, -1.0]) / (0.25 * 0.25)
    np.testing.assert_almost_equal(model.hyper_softness(0.25), expected, decimal=8)
