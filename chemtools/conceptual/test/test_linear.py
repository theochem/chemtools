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
"""Test chemtools.conceptual.linear Module."""

import numpy as np
from chemtools.conceptual.linear import LinearLocalTool


def test_local_linear_fake_density():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model
    model = LinearLocalTool(d0, dp, dm, 10)
    # check density
    np.testing.assert_equal(model.n0, 10.)
    np.testing.assert_almost_equal(model.density_zero, d0, decimal=6)
    np.testing.assert_almost_equal(model.density_plus, dp, decimal=6)
    np.testing.assert_almost_equal(model.density_minus, dm, decimal=6)
    np.testing.assert_almost_equal(model.density(10.), d0, decimal=6)
    np.testing.assert_almost_equal(model.density(11.), dp, decimal=6)
    np.testing.assert_almost_equal(model.density(9.), dm, decimal=6)
    np.testing.assert_almost_equal(model.density(None), d0, decimal=6)
    np.testing.assert_almost_equal(model.density(10.50), 0.5 * d0 + 0.5 * dp, decimal=6)
    np.testing.assert_almost_equal(model.density(10.20), 0.8 * d0 + 0.2 * dp, decimal=6)
    np.testing.assert_almost_equal(model.density(10.75), 0.25 * d0 + 0.75 * dp, decimal=6)
    np.testing.assert_almost_equal(model.density(9.50), 0.5 * dm + 0.5 * d0, decimal=6)
    np.testing.assert_almost_equal(model.density(9.32), 0.68 * dm + 0.32 * d0, decimal=6)
    np.testing.assert_almost_equal(model.density(9.61), 0.39 * dm + 0.61 * d0, decimal=6)


def test_local_linear_fake_fukui_function():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # build a linear local model
    model = LinearLocalTool(d0, dp, dm, 10)
    # check descriptor
    expected = np.array([-0.5, 1.5, 1.0, -1.0, -2.0])
    np.testing.assert_almost_equal(model.ff_plus, expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(10.10), expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(10.73), expected, decimal=6)
    expected = np.array([0.0, -1.0, 2.0, 0.0, -1.0])
    np.testing.assert_almost_equal(model.ff_minus, expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(9.40), expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(9.95), expected, decimal=6)
    expected = np.array([-0.25, 0.25, 1.5, -0.5, -1.5])
    np.testing.assert_almost_equal(model.ff_zero, expected, decimal=6)
    np.testing.assert_almost_equal(model.fukui_function(None), expected, decimal=6)
    np.testing.assert_equal(model.fukui_function(10.), expected)
