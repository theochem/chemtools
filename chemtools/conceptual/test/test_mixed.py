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

from numpy.testing import assert_almost_equal

from chemtools.conceptual.mixed import MixedGlobalTool


def test_mixed_global_mu_gcv_h2o():
    # energy values of H2O cation, neutral & anion at ub3lyp/ccpvtz
    e09, e10, e11 = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea = e09 - e10, e10 - e11
    model = MixedGlobalTool({9: e09, 10: e10, 11: e11})
    # check mu GCV
    expected = - np.array([ip + 3 * ea, 3 * ip + ea]) / 4.
    assert_almost_equal(model.chemical_potential_gcv, expected, decimal=8)
    # check omega GCV
    expected = np.array([ip + 3 * ea, 3 * ip + ea])**2 / (32. * (ip - ea))
    assert_almost_equal(model.electron_transfer_power_gcv, expected, decimal=8)


def test_mixed_global_mu_ma_h2o():
    # energy values of H2O cation, neutral & anion at ub3lyp/ccpvtz
    e09, e10, e11 = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea = e09 - e10, e10 - e11
    model = MixedGlobalTool({9: e09, 10: e10, 11: e11})
    expected = - np.array([ea, ip]) / 1.
    assert_almost_equal(model.chemical_potential_ma(0.), expected, decimal=8)
    expected = - np.array([ip + ea, ip + ea]) / 2.
    assert_almost_equal(model.chemical_potential_ma(1.), expected, decimal=8)
    expected = - np.array([1.5 * ip + ea, ip + 1.5 * ea]) / 2.5
    assert_almost_equal(model.chemical_potential_ma(1.5), expected, decimal=8)
    expected = - np.array([2.0 * ip + ea, ip + 2.0 * ea]) / 3.0
    assert_almost_equal(model.chemical_potential_ma(2.), expected, decimal=8)
    expected = - np.array([3.0 * ip + ea, ip + 3.0 * ea]) / 4.
    assert_almost_equal(model.chemical_potential_ma(3.0), expected, decimal=8)
    expected = - np.array([3.61 * ip + ea, ip + 3.61 * ea]) / 4.61
    assert_almost_equal(model.chemical_potential_ma(3.61), expected, decimal=8)
    expected = - np.array([4.25 * ip + ea, ip + 4.25 * ea]) / 5.25
    assert_almost_equal(model.chemical_potential_ma(4.25), expected, decimal=8)
