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
# pragma pylint: disable=invalid-name,bad-whitespace
"""Test chemtools.denstools.densbased."""


from numpy.testing import assert_raises
import numpy as np
from horton import BeckeMolGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensityBasedTool
from chemtools.utils.cube import CubeGen
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
    assert_raises(ValueError, DensityBasedTool, np.array([[0.], [0.]]), g)
    assert_raises(ValueError, DensityBasedTool, d, np.array([0.]))
    assert_raises(ValueError, DensityBasedTool, d, g, lap=np.array([0.]))


def test_density_local_tool_fake():
    # fake density, gradient and laplacian arrays
    d = np.array([1.00, 3.00, 5.00, 2.00, 7.00])
    g = np.array([[ 0.50,  0.50,  0.50],
                  [ 0.35, -0.35,  0.40],
                  [-0.30, -0.50, -0.50],
                  [ 0.40,  0.40,  0.60],
                  [ 0.25, -0.10, -0.50]])
    l = np.array([1.5, 0.0, -0.4, 1.0, -1.75])
    # build a model
    model = DensityBasedTool(d, g, l, kin=None)
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
    np.testing.assert_almost_equal(model.kinetic_energy_density_weizsacker, expected, decimal=6)
    expected = np.array([2.871234, 17.91722219, 41.97769574, 9.115599745, 73.5470608])
    np.testing.assert_almost_equal(model.kinetic_energy_density_thomas_fermi, expected, decimal=6)


def test_density_local_tool_electrostatic_potential():
    with path('chemtools.data', 'water_b3lyp_sto3g.fchk') as file_path:
        mol = Molecule.from_file(str(file_path))
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='coarse', random_rotate=False, mode='keep')

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # build a density local model
    model = DensityBasedTool(mol.compute_density(grid.points), mol.compute_gradient(grid.points))

    # mep results obtained from Fortran code:
    expected = np.array([-0.01239766, -0.02982537, -0.02201149,   -0.01787292, -0.05682143,
                         -0.02503563, -0.00405942, -0.00818772,   -0.00502268,  0.00321181,
                         -0.03320573, -0.02788605,  0.02741914, 1290.21135500, -0.03319778,
                          0.01428660,  0.10127092,  0.01518299,    0.01530548,  0.00197975,
                         -0.00894206,  0.04330806,  0.03441681,   -0.00203017,  0.02272626,
                          0.03730846,  0.01463959])

    test = model.compute_electrostatic_potential(mol.numbers, mol.coordinates, grid.weights,
                                                 grid.points, cube.points)

    np.testing.assert_almost_equal(test, expected, decimal=1)

    # check ValueError
    assert_raises(ValueError, model.compute_electrostatic_potential, np.array([0.]),
                  mol.coordinates, grid.weights, grid.points, cube.points)
    assert_raises(ValueError, model.compute_electrostatic_potential, mol.numbers, np.array([0.]),
                  grid.weights, grid.points, cube.points)
    assert_raises(ValueError, model.compute_electrostatic_potential, mol.numbers, mol.coordinates,
                  np.array([0.]), grid.points, cube.points)
    assert_raises(ValueError, model.compute_electrostatic_potential, mol.numbers, mol.coordinates,
                  grid.weights, np.array([0.]), cube.points)
    assert_raises(ValueError, model.compute_electrostatic_potential, mol.numbers, mol.coordinates,
                  grid.weights, grid.points, np.array([0.]))
    assert_raises(ValueError, model.compute_electrostatic_potential, mol.numbers, mol.coordinates,
                  grid.weights, grid.points, np.array([[0.]]))
