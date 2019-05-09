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
"""Test chemtools.topology.eigenvalues module."""


import warnings
import numpy as np
from numpy.testing import assert_almost_equal, assert_raises, assert_equal

from chemtools.topology.eigenvalues import EigenDescriptor


def test_raises():
    eigenvalues = [[]]
    assert_raises(TypeError, EigenDescriptor, eigenvalues)
    assert_raises(TypeError, EigenDescriptor, np.array([5.]))


def test_ellipticity():
    eigenvalues = np.array([[10., 20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.],
                            [30., 10., 20.]])
    result = EigenDescriptor(eigenvalues).ellipticity
    assert_almost_equal(result, [0.5, 0., -0.6, np.inf, 0.5])


def test_bond_descriptor():
    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.]])
    result = EigenDescriptor(eigenvalues).bond_descriptor
    # TODO: Check [30. / (-30. / 2.), np.nan, 0., 1. / -10.])
    assert_equal(result, [30. / (-30. / 2.), np.nan, np.nan, 1. / -10.])


def test_eccentricity():
    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.]])
    result = EigenDescriptor(eigenvalues).eccentricity
    assert_almost_equal(result, [np.nan, 2.**0.5, (-2. / -10.)**0.5, np.nan])


def test_index_critical_pt():
    eigenvalues = np.array([[-10., -20., 3.], [1., 2., 3.],
                            [10., -20., -1e-10], [10., -20., -1e-20],
                            [10., -20., 0.]])
    result = EigenDescriptor(eigenvalues).index
    assert_equal(result, [2, 0, 2, 1, 1])


def test_morse_critical_pt():
    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.],
                            [-1e-20, 1e-20, 1e-20],
                            [1e-5, 1e-5, 1e-5]])
    result = EigenDescriptor(eigenvalues).morse
    # Catch warning about zero eigenvalue.
    # with warnings.catch_warnings(record=True) as warn_msg:
    #     warnings.simplefilter("always")
    #     assert_equal(eigen_obj.morse_critical_pt(4), (0., 0.))
    #     assert len(warn_msg) > 0.
    assert_equal(result[0], (3, -1))
    assert_equal(result[1], (3, 3))
    assert_equal(result[2], (3, -3))
    assert_equal(result[5], (3., 3.))
