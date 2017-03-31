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


import tempfile, shutil
from chemtools import *
from contextlib import contextmanager


@contextmanager
def tmpdir(name):
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)


def test_analyze_nci_h2o_dimer_wfn():
    file_path = context.get_fn('test/h2o_dimer_pbe_sto3g.wfn')
    # Check against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    dens_cube1_path = context.get_fn('test/h2o_dimer_pbe_sto3g-dens.cube')
    cube = CubeGen.from_cube(dens_cube1_path)
    # Build the NCI tool
    desp = NCI.from_file(file_path, cube)
    # Check against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    grad_cube1_path = context.get_fn('test/h2o_dimer_pbe_sto3g-grad.cube')
    dmol1 = IOData.from_file(dens_cube1_path)
    gmol1 = IOData.from_file(grad_cube1_path)

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_h2o_dimer_fchk') as dn:
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g')
        desp.dump_files(cube2)
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        np.testing.assert_array_almost_equal(dmol1.coordinates, mol2.coordinates, decimal=6)
        np.testing.assert_equal(dmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        np.testing.assert_array_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        np.testing.assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = dmol1.cube_data / dmol1.cube_data
        data2 = mol2.cube_data / dmol1.cube_data
        np.testing.assert_array_almost_equal(data1, data2, decimal=4)
        np.testing.assert_equal(dmol1.pseudo_numbers, mol2.pseudo_numbers)

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        np.testing.assert_array_almost_equal(gmol1.coordinates,  mol2.coordinates, decimal=6)
        np.testing.assert_equal(gmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        np.testing.assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        np.testing.assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = gmol1.cube_data / gmol1.cube_data
        data2 = mol2.cube_data / gmol1.cube_data
        np.testing.assert_array_almost_equal(data1, data2, decimal=4)
        np.testing.assert_equal(gmol1.pseudo_numbers, mol2.pseudo_numbers)


def test_analyze_nci_h2o_dimer_fchk():
    file_path = context.get_fn('test/h2o_dimer_pbe_sto3g.fchk')
    mol = IOData.from_file(file_path)
    # Check against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    dens_cube1_path = context.get_fn('test/h2o_dimer_pbe_sto3g-dens.cube')
    cube = CubeGen.from_cube(dens_cube1_path)
    # Build the NCI tool
    desp = NCI.from_iodata(mol, cube)
    # Check against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    grad_cube1_path = context.get_fn('test/h2o_dimer_pbe_sto3g-grad.cube')
    dmol1 = IOData.from_file(dens_cube1_path)
    gmol1 = IOData.from_file(grad_cube1_path)

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_h2o_dimer_fchk') as dn:
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g')
        desp.dump_files(cube2)
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        np.testing.assert_array_almost_equal(dmol1.coordinates, mol2.coordinates, decimal=6)
        np.testing.assert_equal(dmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        np.testing.assert_array_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        np.testing.assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = dmol1.cube_data / dmol1.cube_data
        data2 = mol2.cube_data / dmol1.cube_data
        np.testing.assert_array_almost_equal(data1, data2, decimal=4)
        np.testing.assert_equal(dmol1.pseudo_numbers, mol2.pseudo_numbers)

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        np.testing.assert_array_almost_equal(gmol1.coordinates,  mol2.coordinates, decimal=6)
        np.testing.assert_equal(gmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        np.testing.assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        np.testing.assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = gmol1.cube_data / gmol1.cube_data
        data2 = mol2.cube_data / gmol1.cube_data
        np.testing.assert_array_almost_equal(data1, data2, decimal=4)
        np.testing.assert_equal(gmol1.pseudo_numbers, mol2.pseudo_numbers)
