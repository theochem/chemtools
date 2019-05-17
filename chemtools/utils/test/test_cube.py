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
"""Test chemtools.utils.cube."""


import shutil
import tempfile
from contextlib import contextmanager
from numpy.testing import assert_raises
import numpy as np
from horton import IOData
from chemtools.toolbox.conceptual import LocalConceptualDFT
from chemtools.utils.cube import UniformGrid
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


@contextmanager
def tmpdir(name):
    """Create temporary directory that gets deleted after accessing it."""
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)


def test_cubegen_o2_uhf():
    with path('chemtools.data', 'o2_uhf.fchk') as path_file:
        mol = IOData.from_file(str(path_file))

    # create cube file from file:
    cube = UniformGrid.from_file(path_file, spacing=0.5, threshold=6.0, rotate=False)

    # test the cube gives the right result:
    origin_result = [-6.0, -6.0, -7.25]
    axes_result = [[0.5, 0.0, 0.0],
                   [0.0, 0.5, 0.0],
                   [0.0, 0.0, 0.5]]
    shape_result = [24, 24, 29]
    weight_result = np.full(cube.npoints, 0.125)

    np.testing.assert_array_almost_equal(cube.origin, origin_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.axes, axes_result, decimal=7)
    np.testing.assert_array_equal(cube.shape, shape_result)
    np.testing.assert_array_equal(cube.numbers, mol.numbers)
    np.testing.assert_array_almost_equal(cube.coordinates, mol.coordinates, decimal=10)
    np.testing.assert_array_almost_equal(cube.pseudo_numbers, mol.pseudo_numbers, decimal=10)
    np.testing.assert_array_almost_equal(cube.weights(), weight_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.weights(method='R0'), weight_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.weights(method='R'), weight_result, decimal=7)

    # create cube file from molecule:
    cube = UniformGrid.from_molecule(mol.numbers, mol.pseudo_numbers, mol.coordinates)

    # test the cube gives the right result:
    origin_result = [-5.0, -5.0, -6.1]
    axes_result = [[0.0, 0.0, 0.2],
                   [0.0, 0.2, 0.0],
                   [0.2, 0.0, 0.0]]
    shape_result = [61, 50, 50]
    weight_result = np.full(cube.npoints, 0.0080)

    np.testing.assert_array_almost_equal(cube.origin, origin_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.axes, axes_result, decimal=7)
    np.testing.assert_array_equal(cube.shape, shape_result)
    np.testing.assert_array_almost_equal(cube.weights(), weight_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.weights(method='R0'), weight_result, decimal=7)
    np.testing.assert_array_almost_equal(cube.weights(method='R'), weight_result, decimal=7)

    # test integration of Fukui functions:

    tool = LocalConceptualDFT.from_file(str(path_file), model='linear', points=cube.points)

    ffm_default = cube.integrate(tool.ff_minus)
    ffm_r = cube.integrate(tool.ff_minus, method='R')
    ffm_r0 = cube.integrate(tool.ff_minus, method='R0')

    np.testing.assert_almost_equal(ffm_default, ffm_r, decimal=7)
    np.testing.assert_almost_equal(ffm_r, ffm_r0, decimal=7)
    np.testing.assert_almost_equal(ffm_default, 1.000, decimal=2)
    np.testing.assert_almost_equal(ffm_default, 1.000, decimal=2)
    np.testing.assert_almost_equal(ffm_default, 1.000, decimal=2)

    o = np.array([-3.0, -3.0, -3.0])
    a = np.array([[2.0,  0.0,  0.0],
                  [0.0,  2.0,  0.0],
                  [0.0,  0.0,  2.0]])
    s = np.array([ 3,    3,    3])
    # check ValueError
    assert_raises(ValueError, UniformGrid, mol.numbers, mol.pseudo_numbers, mol.coordinates,
                  np.array([0.]), a, s)
    assert_raises(ValueError, UniformGrid, mol.numbers, mol.pseudo_numbers, mol.coordinates,
                  o, np.array([0.]), s)
    assert_raises(ValueError, UniformGrid, mol.numbers, mol.pseudo_numbers, mol.coordinates,
                  o, a, np.array([0.]))
    assert_raises(ValueError, UniformGrid.from_cube, 'test.wrong_end')
    assert_raises(ValueError, cube.dump_cube, 'test.wrong_end', tool.ff_minus)
    assert_raises(ValueError, cube.dump_cube, 'test.cube', np.array([0.]))
    assert_raises(ValueError, cube.weights, method='erroneous')
    assert_raises(ValueError, cube.integrate, np.array([0.]))


def test_cube_h2o_dimer():
    with path('chemtools.data', 'h2o_dimer_pbe_sto3g-dens.cube') as file_path:
    # Build the cube
    # Check against previous generated .cube files
        cube = UniformGrid.from_cube(file_path)
        mol1 = IOData.from_file(str(file_path))

    with tmpdir('chemtools.test.test_base.test_cube_h2o_dimer') as dn:
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        cube.dump_cube(cube2, mol1.cube_data)
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        np.testing.assert_array_almost_equal(mol1.coordinates, mol2.coordinates, decimal=6)
        np.testing.assert_equal(mol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = mol1.grid
        ugrid2 = mol2.grid
        np.testing.assert_array_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        np.testing.assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = mol1.cube_data / mol1.cube_data
        data2 = mol2.cube_data / mol1.cube_data
        np.testing.assert_array_almost_equal(data1, data2, decimal=4)
        np.testing.assert_equal(mol1.pseudo_numbers, mol2.pseudo_numbers)
