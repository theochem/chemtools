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

r"""Test chemtools.utils.eigenvalues module."""

import warnings
import numpy as np
from numpy.testing import assert_almost_equal, assert_raises, assert_equal

from chemtools.utils.eigenvalues import EigenDescriptor


def test_raises():
    eigenvalues = [[]]
    assert_raises(TypeError, EigenDescriptor, eigenvalues)
    assert_raises(TypeError, EigenDescriptor, np.array([5.]))


def test_ellipticity():
    eigenvalues = np.array([[10., 20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.],
                            [30., 10., 20.]])
    eigen_obj = EigenDescriptor(eigenvalues)
    actual = eigen_obj.ellipticity(0)
    assert_equal(actual, (eigenvalues[0, 2] / eigenvalues[0, 1]) - 1.)

    actual = eigen_obj.ellipticity(1)
    assert_equal(actual, (eigenvalues[1, 0] / eigenvalues[1, 2]) - 1.)

    actual = eigen_obj.ellipticity(2)
    assert_equal(actual, (eigenvalues[2, 1] / eigenvalues[2, 2]) - 1.)

    actual = eigen_obj.ellipticity(3)
    assert_equal(actual, np.inf)

    assert_equal(eigen_obj.ellipticity(4), (30. / 20.) - 1.)


def test_bond_descriptor():
    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.]])
    eigen_obj = EigenDescriptor(eigenvalues)

    actual = eigen_obj.bond_descriptor(0)
    assert_equal(actual, 30. / (-30. / 2.))

    actual = eigen_obj.bond_descriptor(1)
    assert_equal(actual, None)

    actual = eigen_obj.bond_descriptor(2)
    assert_almost_equal(actual, 0.)

    actual = eigen_obj.bond_descriptor(3)
    assert_equal(actual, 1. / -10.)


def test_eccentricity():
    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.]])
    eigen_obj = EigenDescriptor(eigenvalues)

    assert_equal(eigen_obj.eccentricity(0), None)
    assert_almost_equal(eigen_obj.eccentricity(1), 2.**0.5)
    assert_equal(eigen_obj.eccentricity(2), (-2. / -10.)**0.5)
    assert_equal(eigen_obj.eccentricity(3), None)


def test_index_critical_pt():
    eigenvalues = np.array([[-10., -20., 3.], [1., 2., 3.],
                            [10., -20., -1e-10], [10., -20., -1e-20],
                            [10., -20., 0.]])
    eigen_obj = EigenDescriptor(eigenvalues)
    assert_equal(eigen_obj.index_critical_pt(0), 2)
    assert_equal(eigen_obj.index_critical_pt(1), 0)
    assert_equal(eigen_obj.index_critical_pt(2), 2)
    assert_equal(eigen_obj.index_critical_pt(3), 1)
    assert_equal(eigen_obj.index_critical_pt(4), 1)


def test_morse_critical_pt():

    eigenvalues = np.array([[-10., -20., 30.], [20., 10., 20.],
                            [-10., -2., -5.], [-10., 0., 1.],
                            [-1e-20, 1e-20, 1e-20],
                            [1e-5, 1e-5, 1e-5]])
    eigen_obj = EigenDescriptor(eigenvalues)

    # Catch warning about zero eigenvalue.
    with warnings.catch_warnings(record=True) as warn_msg:
        warnings.simplefilter("always")
        assert_equal(eigen_obj.morse_critical_pt(4), (0., 0.))
        assert len(warn_msg) > 0.
    assert_equal(eigen_obj.morse_critical_pt(0), (3, -1))
    assert_equal(eigen_obj.morse_critical_pt(1), (3, 3))
    assert_equal(eigen_obj.morse_critical_pt(2), (3, -3))
    assert_equal(eigen_obj.morse_critical_pt(5), (3., 3.))
