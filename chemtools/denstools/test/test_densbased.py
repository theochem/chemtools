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
# pragma pylint: disable=invalid-name,bad-whitespace
"""Test chemtools.denstools.densbased."""


from numpy.testing import assert_raises, assert_almost_equal
import numpy as np
from chemtools.denstools.densbased import DensTool, DensGradTool
from chemtools.denstools.densbased import DensGradLapTool, DensGradLapKedTool
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_densbased_raises():
    # fake dens & grad arrays
    d = np.array([1.00, 3.00, 5.00, 2.00, 7.00])
    g = np.array([[ 0.50,  0.50,  0.50],
                  [ 0.35, -0.35,  0.40],
                  [-0.30, -0.50, -0.50],
                  [ 0.40,  0.40,  0.60],
                  [ 0.25, -0.10, -0.50]])
    # check ValueError
    assert_raises(ValueError, DensGradLapTool, np.array([[0.], [0.]]), g, None)
    assert_raises(ValueError, DensGradLapTool, d, np.array([0.]), None)
    assert_raises(ValueError, DensGradLapTool, d, g, lap=np.array([0.]))


def test_dens_based_fake():
    # fake density, gradient and laplacian arrays
    d = np.array([1.00, 3.00, 5.00, 2.00, 7.00])
    # build a model
    model = DensTool(d)
    # check attributes
    np.testing.assert_almost_equal(model.density, d, decimal=6)
    # check Shannon information
    expected = np.array([0.00000000, 3.29583687, 8.04718956, 1.38629436, 13.62137104])
    np.testing.assert_almost_equal(model.shannon_information, expected, decimal=6)
    # check TF kinetic energy density
    expected = np.array([2.871234, 17.91722219, 41.97769574, 9.115599745, 73.5470608])
    np.testing.assert_almost_equal(model.ked_thomas_fermi, expected, decimal=6)


def test_dens_grad_based_fake():
    # fake density, gradient and laplacian arrays
    d = np.array([1.00, 3.00, 5.00, 2.00, 7.00])
    g = np.array([[ 0.50,  0.50,  0.50],
                  [ 0.35, -0.35,  0.40],
                  [-0.30, -0.50, -0.50],
                  [ 0.40,  0.40,  0.60],
                  [ 0.25, -0.10, -0.50]])
    # build a model
    model = DensGradTool(d, g)
    # check attributes
    np.testing.assert_almost_equal(model.density, d, decimal=6)
    np.testing.assert_almost_equal(model.gradient, g, decimal=6)
    # check Shannon information
    expected = np.array([0.00000000, 3.29583687, 8.04718956, 1.38629436, 13.62137104])
    np.testing.assert_almost_equal(model.shannon_information, expected, decimal=6)
    # check gradient norm
    expected = np.array([0.86602540, 0.63639610, 0.76811457, 0.82462113, 0.56789083])
    np.testing.assert_almost_equal(model.gradient_norm, expected, decimal=6)
    # check reduced density gradient
    expected = np.array([0.13996742, 0.02377181, 0.01451986, 0.05289047, 0.00685431])
    np.testing.assert_almost_equal(model.reduced_density_gradient, expected, decimal=6)
    # check Weizsacker kinetic energy
    expected = np.array([0.09375000, 0.01687500, 0.01475000, 0.04250000, 0.00575893])
    np.testing.assert_almost_equal(model.ked_weizsacker, expected, decimal=6)
    # check TF kinetic energy density
    expected = np.array([2.871234, 17.91722219, 41.97769574, 9.115599745, 73.5470608])
    np.testing.assert_almost_equal(model.ked_thomas_fermi, expected, decimal=6)


def test_dens_grad_lap_based_fake():
    # fake density, gradient and laplacian arrays
    d = np.array([1.00, 3.00, 5.00, 2.00, 7.00])
    g = np.array([[ 0.50,  0.50,  0.50],
                  [ 0.35, -0.35,  0.40],
                  [-0.30, -0.50, -0.50],
                  [ 0.40,  0.40,  0.60],
                  [ 0.25, -0.10, -0.50]])
    l = np.array([1.5, 0.0, -0.4, 1.0, -1.75])
    # build a model
    model = DensGradLapTool(d, g, l)
    # check attributes
    np.testing.assert_almost_equal(model.density, d, decimal=6)
    np.testing.assert_almost_equal(model.gradient, g, decimal=6)
    np.testing.assert_almost_equal(model.laplacian, l, decimal=6)
    # check Shannon information
    expected = np.array([0.00000000, 3.29583687, 8.04718956, 1.38629436, 13.62137104])
    np.testing.assert_almost_equal(model.shannon_information, expected, decimal=6)
    # check gradient norm
    expected = np.array([0.86602540, 0.63639610, 0.76811457, 0.82462113, 0.56789083])
    np.testing.assert_almost_equal(model.gradient_norm, expected, decimal=6)
    # check reduced density gradient
    expected = np.array([0.13996742, 0.02377181, 0.01451986, 0.05289047, 0.00685431])
    np.testing.assert_almost_equal(model.reduced_density_gradient, expected, decimal=6)
    # check Weizsacker kinetic energy
    expected = np.array([0.09375000, 0.01687500, 0.01475000, 0.04250000, 0.00575893])
    np.testing.assert_almost_equal(model.ked_weizsacker, expected, decimal=6)
    # check TF kinetic energy density
    expected = np.array([2.871234, 17.91722219, 41.97769574, 9.115599745, 73.5470608])
    np.testing.assert_almost_equal(model.ked_thomas_fermi, expected, decimal=6)


def test_dens_grad_lap_ked_based_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    # check local properties at the position of nuclei
    args = (data['nuc_dens'], data['nuc_grad'], data['nuc_lap'], data['nuc_ked_pd'])
    model = DensGradLapKedTool(*args)
    assert_almost_equal(model.density, data['nuc_dens'], decimal=6)
    assert_almost_equal(model.gradient, data['nuc_grad'], decimal=6)
    assert_almost_equal(model.laplacian, data['nuc_lap'], decimal=6)
    assert_almost_equal(model.ked_positive_definite, data['nuc_ked_pd'], decimal=6)
    assert_almost_equal(model.gradient_norm, data['nuc_grad_norm'], decimal=6)
    # assert_almost_equal(model.reduced_density_gradient, data['nuc_rdg'], decimal=6)
    assert_almost_equal(model.ked_hamiltonian, data['nuc_ked_ham'], decimal=4)
    assert_almost_equal(model.ked_general(0.), data['nuc_ked_ham'], decimal=4)
