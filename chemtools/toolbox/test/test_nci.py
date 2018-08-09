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
"""Test chemtools.analysis.nci."""


import os
import shutil
import tempfile
from contextlib import contextmanager

import numpy as np
from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import IOData
from chemtools import context
from chemtools.utils import CubeGen
from chemtools.toolbox.nci import NCI
from chemtools.toolbox.molecule import make_molecule


@contextmanager
def tmpdir(name):
    """Create temporary directory that gets deleted after accessing it."""
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
        assert_almost_equal(dmol1.coordinates, mol2.coordinates, decimal=6)
        assert_equal(dmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = dmol1.cube_data / dmol1.cube_data
        data2 = mol2.cube_data / dmol1.cube_data
        assert_almost_equal(data1, data2, decimal=4)
        assert_equal(dmol1.pseudo_numbers, mol2.pseudo_numbers)

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        assert_almost_equal(gmol1.coordinates, mol2.coordinates, decimal=6)
        assert_equal(gmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = gmol1.cube_data / gmol1.cube_data
        data2 = mol2.cube_data / gmol1.cube_data
        assert_almost_equal(data1, data2, decimal=4)
        assert_equal(gmol1.pseudo_numbers, mol2.pseudo_numbers)


def test_analyze_nci_h2o_dimer_fchk():
    file_path = context.get_fn('test/h2o_dimer_pbe_sto3g.fchk')
    mol = make_molecule(file_path)
    # Check against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    dens_cube1_path = context.get_fn('test/h2o_dimer_pbe_sto3g-dens.cube')
    cube = CubeGen.from_cube(dens_cube1_path)
    # Build the NCI tool
    desp = NCI.from_molecule(mol, cube)
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
        assert_almost_equal(dmol1.coordinates, mol2.coordinates, decimal=6)
        assert_equal(dmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = dmol1.cube_data / dmol1.cube_data
        data2 = mol2.cube_data / dmol1.cube_data
        assert_almost_equal(data1, data2, decimal=4)
        assert_equal(dmol1.pseudo_numbers, mol2.pseudo_numbers)

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)
        # Check coordinates
        assert_almost_equal(gmol1.coordinates, mol2.coordinates, decimal=6)
        assert_equal(gmol1.numbers, mol2.numbers)
        # Check grid data
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        assert_almost_equal(ugrid1.grid_rvecs, ugrid2.grid_rvecs, decimal=6)
        assert_equal(ugrid1.shape, ugrid2.shape)
        data1 = gmol1.cube_data / gmol1.cube_data
        data2 = mol2.cube_data / gmol1.cube_data
        assert_almost_equal(data1, data2, decimal=4)
        assert_equal(gmol1.pseudo_numbers, mol2.pseudo_numbers)


def test_analyze_nci_assert_errors():
    file_path = context.get_fn('test/h2o_dimer_pbe_sto3g.fchk')
    mol = make_molecule(file_path)
    cube = CubeGen.from_file(file_path, spacing=2., threshold=0.0)

    dens = np.array([2.10160232e-04, 1.11307672e-05, 3.01244062e-04, 2.31768360e-05,
                     6.56282686e-03, 2.62815892e-04, 2.46559574e-02, 1.82760928e-03,
                     1.89299475e-02, 1.39689069e-03, 1.10257641e+00, 3.63942662e-02,
                     8.01150391e-03, 2.79542971e-04, 1.98278511e-02, 1.89336116e-03])

    rdg = np.array([6.434055294420, 19.330749118092, 5.927230664196, 13.593077700571,
                    2.411123457672,  5.203993648485, 1.482717816902,  3.136335572700,
                    1.779261001623,  3.395839461280, 0.226436405120,  0.912678557191,
                    2.223911275208,  4.990189542067, 1.676113282597,  3.171756800841])

    assert_raises(ValueError, NCI, np.array([0.]), rdg, cube)
    assert_raises(ValueError, NCI, dens, np.array([0.]), cube)
    assert_raises(ValueError, NCI, dens, rdg, cube, hessian=np.array([0.]))
    assert_raises(ValueError, NCI.from_file, [file_path, file_path])
    assert_raises(ValueError, NCI.from_file, file_path, cube=1)

    desp = NCI(dens, rdg, cube)
    assert desp.signed_density is None
    assert desp.eigvalues is None

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_assert_errors') as dn:
        test = '%s/%s' % (dn, 'test')
        desp.dump_files(test)
        test = '%s/%s' % (dn, 'test-dens.cube')
        assert os.path.isfile(test) and os.access(test, os.R_OK)
        test = '%s/%s' % (dn, 'test-grad.cube')
        assert os.path.isfile(test) and os.access(test, os.R_OK)
        test = '%s/%s' % (dn, 'test.vmd')
        assert os.path.isfile(test) and os.access(test, os.R_OK)

    desp = NCI.from_molecule(mol)
    assert desp.signed_density.shape == desp._density.shape

    desp = NCI.from_file(file_path, cube)

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_assert_errors') as dn:
        test = '%s/%s' % (dn, 'test.png')
        desp.plot(test)
        assert os.path.isfile(test) and os.access(test, os.R_OK)
