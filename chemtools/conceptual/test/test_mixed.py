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
"""Test chemtools.conceptual.mixed Module."""


import numpy as np

from numpy.testing import assert_almost_equal, assert_raises

from chemtools.conceptual.mixed import MixedGlobalTool, MixedLocalTool, MixedCondensedTool


def test_mixed_raises():
    # check global raises
    assert_raises(ValueError, MixedGlobalTool, {0.0: -0.5, 0.5: -1.0, 1.0: -0.75})
    assert_raises(ValueError, MixedGlobalTool, {0.0: -0.5, 2.0: -1.0, -1.0: -0.75})
    # check local raises
    dict_d = {1.0: np.array([0.1, 0.3]), 2.0: np.array([0.2, 0.5]), 3.0: np.array([0.2, 0.6])}
    dict_e = {1.0: -0.5, 2.0: -1.0, 4.0: -0.75}
    assert_raises(ValueError, MixedLocalTool, dict_e, dict_d)
    dict_d = {1.0: np.array([0.1, 0.3]), 2.0: np.array([0.2, 0.5]), 4.0: np.array([0.2, 0.6])}
    assert_raises(ValueError, MixedLocalTool, dict_e, dict_d)
    # check condensed raises
    dict_p = {1.0: np.array([0.1, 0.3]), 2.0: np.array([0.2, 0.5]), 3.0: np.array([0.2, 0.6])}
    dict_e = {1.0: -0.5, 2.0: -1.0, 4.0: -0.75}
    assert_raises(ValueError, MixedCondensedTool, dict_e, dict_p)
    dict_p = {1.0: np.array([0.1, 0.3]), 2.0: np.array([0.2, 0.5]), 4.0: np.array([0.2, 0.6])}
    assert_raises(ValueError, MixedCondensedTool, dict_e, dict_p)


def test_mixed_global_mu_gcv_h2o():
    # energy values of H2O cation, neutral & anion at ub3lyp/ccpvtz
    e09, e10, e11 = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea = e09 - e10, e10 - e11
    model = MixedGlobalTool({9: e09, 10: e10, 11: e11})
    # check mu GCV
    expected = np.array([ip + 3 * ea, 3 * ip + ea]) / -4.
    assert_almost_equal(model.chemical_potential_gcv, expected, decimal=8)
    # check omega GCV
    expected = np.array([ip + 3 * ea, 3 * ip + ea])**2 / (32. * (ip - ea))
    assert_almost_equal(model.electron_transfer_power_gcv, expected, decimal=8)


def test_mixed_global_mu_ma_h2o():
    # energy values of H2O cation, neutral & anion at ub3lyp/ccpvtz
    e09, e10, e11 = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea = e09 - e10, e10 - e11
    model = MixedGlobalTool({9: e09, 10: e10, 11: e11})
    expected = np.array([ea, ip]) / -1.
    assert_almost_equal(model.chemical_potential_ma(0.), expected, decimal=8)
    expected = np.array([ip + ea, ip + ea]) / -2.
    assert_almost_equal(model.chemical_potential_ma(1.), expected, decimal=8)
    expected = np.array([1.5 * ip + ea, ip + 1.5 * ea]) / -2.5
    assert_almost_equal(model.chemical_potential_ma(1.5), expected, decimal=8)
    expected = np.array([2.0 * ip + ea, ip + 2.0 * ea]) / -3.0
    assert_almost_equal(model.chemical_potential_ma(2.), expected, decimal=8)
    expected = np.array([3.0 * ip + ea, ip + 3.0 * ea]) / -4.
    assert_almost_equal(model.chemical_potential_ma(3.0), expected, decimal=8)
    expected = np.array([3.61 * ip + ea, ip + 3.61 * ea]) / -4.61
    assert_almost_equal(model.chemical_potential_ma(3.61), expected, decimal=8)
    expected = np.array([4.25 * ip + ea, ip + 4.25 * ea]) / -5.25
    assert_almost_equal(model.chemical_potential_ma(4.25), expected, decimal=8)


def test_mixed_local_fake():
    # fake density & energy values
    dict_d = {1: np.array([0, 1.2, 0.7]), 2: np.array([3.1, 0.4, 0]), 3: np.array([5.6, 0.2, 1.4])}
    dict_e = {1: -0.5, 2: -1.0, 3: -0.75}
    ip, ea = -0.5 - (-1.0), -1.0 - (-0.75)
    ffp, ff0, ffm = dict_d[3] - dict_d[2], 0.5 * (dict_d[3] - dict_d[1]), dict_d[2] - dict_d[1]
    f2 = dict_d[3] - 2 * dict_d[2] + dict_d[1]
    # build mixed local model
    model = MixedLocalTool(dict_e, dict_d)
    # check softness
    assert_almost_equal(model.softness_yp[0], ffp / (ip - ea), decimal=8)
    assert_almost_equal(model.softness_yp[1], ff0 / (ip - ea), decimal=8)
    assert_almost_equal(model.softness_yp[2], ffm / (ip - ea), decimal=8)
    # check philicity CMS
    assert_almost_equal(model.philicity_cms[0], ffp * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    assert_almost_equal(model.philicity_cms[1], ff0 * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    assert_almost_equal(model.philicity_cms[2], ffm * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    # check philicity MGV
    pp = ea * ffp / (ip - ea) + 0.5 * ea**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[0], pp, decimal=8)
    p0 = 0.5 * (ip + ea) * ff0 / (ip - ea) + 0.5**3 * (ip + ea)**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[1], p0, decimal=8)
    pm = -ip * ffm / (ip - ea) + 0.5 * ip**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[2], pm, decimal=8)


def test_mixed_condensed_eap_h2o():
    # energy values of H2O cation, neutral & anion at ub3lyp/ccpvtz
    e09, e10, e11 = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea = e09 - e10, e10 - e11
    # ESP population values of H2O cation, nautral & anion at ub3lyp/ccpvtz
    pop09 = np.array([8, 1, 1]) - np.array([4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    pop10 = np.array([8, 1, 1]) - np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    pop11 = np.array([8, 1, 1]) - np.array([-5.81613550E-01, -2.09193820E-01, -2.09192630E-01])
    ffp, ff0, ffm = pop11 - pop10, 0.5 * (pop11 - pop09), pop10 - pop09
    f2 = pop11 - 2 * pop10 + pop09
    # build mixed condensed model
    model = MixedCondensedTool({9: e09, 10: e10, 11: e11}, {9: pop09, 10: pop10, 11: pop11})
    # check softness
    assert_almost_equal(model.softness_yp[0], ffp / (ip - ea), decimal=8)
    assert_almost_equal(model.softness_yp[1], ff0 / (ip - ea), decimal=8)
    assert_almost_equal(model.softness_yp[2], ffm / (ip - ea), decimal=8)
    # check philicity CMS
    assert_almost_equal(model.philicity_cms[0], ffp * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    assert_almost_equal(model.philicity_cms[1], ff0 * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    assert_almost_equal(model.philicity_cms[2], ffm * (ip + ea)**2 / (8 * (ip - ea)), decimal=8)
    # check philicity MGV
    pp = ea * ffp / (ip - ea) + 0.5 * ea**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[0], pp, decimal=8)
    p0 = 0.5 * (ip + ea) * ff0 / (ip - ea) + 0.5**3 * (ip + ea)**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[1], p0, decimal=8)
    pm = -ip * ffm / (ip - ea) + 0.5 * ip**2 * f2 / (ip - ea)**2
    assert_almost_equal(model.philicity_mgvgc[2], pm, decimal=8)
