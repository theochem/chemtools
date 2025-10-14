# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2024 The ChemTools Development Team
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

"""Test chemtools.toolbox.iqa."""

import glob
import pytest
import numpy as np
from numpy.testing import assert_raises, assert_almost_equal
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal, assert_allclose

from grid.onedgrid import UniformInteger, GaussChebyshev
from grid.rtransform import ExpRTransform, PowerRTransform, BeckeRTransform
from grid.atomgrid import AtomGrid
from iodata import load_one
from gbasis.wrappers import from_iodata
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart
from chemtools.toolbox.iqa import IQA


try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_iqa_raises():

    with path('chemtools.data', 'h2o_rhf_sto3g.fchk')as fname:
        with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname2:
            mol_1 = Molecule.from_file(str(fname))
            mol_2 = Molecule.from_file(str(fname2))
            basis_gbasis = from_iodata(mol_1._iodata)
            one_rdm_gbasis = mol_1._iodata.one_rdms.get("post_scf", mol_1._iodata.one_rdms.get("scf"))
            # wrong_grid = AtomicGrid(6, 6, np.array([0., 0., 0]))
            onedg = UniformInteger(100)
            rgrid = ExpRTransform(1e-5, 2e1).transform_1d_grid(onedg)
            wrong_grid = AtomGrid(rgrid, center=mol_1.coordinates[0])
            wrong_grid2 = MolecularGrid.from_file(fname2)
            grid = MolecularGrid.from_molecule(mol_1, specs="insane", k=3,
                                                   rotate=False)
            part = DensPart.from_molecule(mol_1, grid=grid, scheme="h", local=False)
            # Check invalid grid
            # assert_raises(ValueError, IQA.from_file, str('h2o_rhf_sto3g.fchk'), 'Atomic')
            assert_raises(TypeError, IQA.from_molecule, mol_1,  wrong_grid)
            assert_raises(TypeError, IQA, mol_1, basis_gbasis[0], one_rdm_gbasis, wrong_grid, part)
            assert_raises(TypeError, IQA, mol_1, basis_gbasis[0], one_rdm_gbasis, wrong_grid2, part)
            # Check wrong molecule
            assert_raises(ValueError, IQA.from_molecule, 'wrong_mol_iodata', grid)
            # Check part/scheme
            assert_raises(TypeError, IQA.from_molecule, mol_1, grid, part='wrong_part')
            assert_raises(TypeError, IQA, mol_1, basis_gbasis, one_rdm_gbasis, grid, part='wrong_part')
            # Check basis
            assert_raises(TypeError, IQA, mol_1, ['wrong_basis'], one_rdm_gbasis, grid, part='wrong_part')
            assert_raises(TypeError, IQA, mol_1, basis_gbasis, np.array([[1, 2], [1, 2]], dtype=bool), grid, part='wrong_part')
            assert_raises(TypeError, IQA, mol_1, basis_gbasis, np.array([1.0, 2.0, 3.0]), grid, part='wrong_part')
            assert_raises(ValueError, IQA, mol_1, basis_gbasis, np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]), grid, part='wrong_part')
            assert_raises(ValueError, IQA, mol_1, basis_gbasis, np.array([[1.0, 2.0], [3.0, 4.0]]), grid, part='wrong_part')


def test_h2o_rhf_sto3g():
    # check total values of decomposition against in h2o_rhf_sto3g.log
    with path('chemtools.data', 'h2o_rhf_sto3g.fchk') as fname:
        mol_iqa = IQA.from_file(fname, scheme='H', ee_interatomic=False, threshold=1e-1)
        results_iqa = mol_iqa.run_atomic()
        assert_allclose(results_iqa['nn_total'],  9.1559536481, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['en_atomic']),  -1.968747462907e+02, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['kin_atomic']),  7.457892107183e+01, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['coul_atomic']),  47.276218, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['x_atomic']),  -9.100148, rtol=1.e-4, atol=0.)
        # Check sum of atomic components equals total values
        assert_allclose(np.sum(results_iqa['en_atomic']), -1.968747462907e+02, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['kin_atomic']), 7.457892107183e+01, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['coul_atomic']), 47.276218, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['x_atomic']), -9.100148, rtol=1.e-4, atol=0.)

@pytest.mark.skip(reason="ignore IQA DFT")
def test_h2o_rpbepbe_ccpvtz():
    # check total values of decomposition against in h2o_rpbepbe_sto3g.log
    with path('chemtools.data', 'h2o_rpbepbe_sto3g.fchk') as fname:
        mol_iqa = IQA.from_file(fname, scheme='H')
        results_iqa = mol_iqa.iqa(dft_exch="gga_x_pbe", dft_corr="gga_c_pbe")
        assert_allclose(results_iqa['nn_total'], 9.1559536481, rtol=1.e-4, atol=0.)
        assert_allclose(results_iqa['en_total'], -1.968736981462e+02    , rtol=1.e-4, atol=0.)
        assert_allclose(results_iqa['kin_total'], 7.459229919735e+01, rtol=1.e-4, atol=0.)
        assert_allclose(results_iqa['c_total'], -0.343334, rtol=1.e-4, atol=0.)
        assert_allclose(results_iqa['x_total'], -9.018861, rtol=1.e-4, atol=0.)
        # Check sum of atomic components equals total values
        assert_allclose(np.sum(results_iqa['en_atomic']), -1.968736981462e+02, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['kin_atomic']), 7.459229919735e+01, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['c_atomic']), -0.343334, rtol=1.e-4, atol=0.)
        assert_allclose(np.sum(results_iqa['x_total']), -9.018861, rtol=1.e-4, atol=0.)
