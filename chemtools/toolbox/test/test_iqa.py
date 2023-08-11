# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2023 The ChemTools Development Team
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
import numpy as np
from numpy.testing import assert_raises, assert_almost_equal
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal, assert_allclose

from horton import ProAtomDB, AtomicGrid
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
        with path('chemtools.data', 'ch3_utpsstpss_321g.fchk') as fname2:
            mol_iodata = load_one(str(fname))
            mol_chemtools = Molecule.from_file(str(fname))
            mol_chemtools2 = Molecule.from_file(str(fname2))
            basis_gbasis = from_iodata(mol_iodata)
            one_rdm_gbasis = mol_iodata.one_rdms.get("post_scf", mol_iodata.one_rdms.get("scf"))
            wrong_grid = AtomicGrid(6, 6, np.array([0., 0., 0]))
            wrong_grid2 = MolecularGrid.from_file('ch3_utpsstpss_321g.fchk')
            grid = MolecularGrid.from_molecule(mol_chemtools, specs="insane", k=3,
                                                   rotate=False)
            proatomdb = ProAtomDB.from_refatoms(mol_chemtools.numbers)
            part = DensPart.from_molecule(mol_chemtools, grid=grid, scheme="h",
                                          proatomdb=proatomdb, local=False)
            # Check invalid grid
            assert_raises(ValueError, IQA.from_file, str(fname), 'Atomic')
            assert_raises(TypeError, IQA.from_molecule, mol_iodata, mol_chemtools, wrong_grid)
            assert_raises(TypeError, IQA, mol_iodata, basis_gbasis[0], one_rdm_gbasis, wrong_grid, part, molecule_chemtools=mol_chemtools)
            assert_raises(TypeError, IQA, mol_iodata, basis_gbasis[0], one_rdm_gbasis, wrong_grid2, part, molecule_chemtools=mol_chemtools)
            # Check wrong molecule iodata
            assert_raises(ValueError, IQA.from_molecule, 'wrong_mol_iodata', mol_chemtools, grid)
            # Check part/scheme
            assert_raises(TypeError, IQA.from_molecule, mol_iodata, mol_chemtools, grid, part='wrong_part')
            assert_raises(TypeError, IQA, mol_iodata, basis_gbasis, one_rdm_gbasis, grid, part='wrong_part', molecule_chemtools=mol_chemtools)
            assert_raises(NotImplementedError, IQA.from_molecule, mol_iodata, mol_chemtools, grid, scheme='HI')
            # Check basis
            assert_raises(TypeError, IQA, mol_iodata, ['wrong_basis'], one_rdm_gbasis, grid, part='wrong_part', molecule_chemtools=mol_chemtools)
            assert_raises(TypeError, IQA, mol_iodata, basis_gbasis, np.array([[1, 2], [1, 2]], dtype=bool), grid, part='wrong_part', molecule_chemtools=mol_chemtools)
            assert_raises(TypeError, IQA, mol_iodata, basis_gbasis, np.array([1.0, 2.0, 3.0]), grid, part='wrong_part', molecule_chemtools=mol_chemtools)
            assert_raises(ValueError, IQA, mol_iodata, basis_gbasis, np.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]), grid, part='wrong_part', molecule_chemtools=mol_chemtools)
            assert_raises(ValueError, IQA, mol_iodata, basis_gbasis, np.array([[1.0, 2.0], [3.0, 4.0]]), grid, part='wrong_part', molecule_chemtools=mol_chemtools)


def test_h2o_rhf_sto3g():
    # check total values of decomposition against in h2o_rhf_sto3g.log
    mol_iqa = IQA.from_file('h2o_rhf_sto3g.fchk', grid_type='becke', scheme='H')
    results_iqa = mol_iqa.iqa()
    assert_allclose(results_iqa['nn_total'],  9.1559536481, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['en_total'],  -1.968747462907e+02, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['kin_total'],  7.457892107183e+01, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['col_total'],  47.276218, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['ex_total'],  -9.100148, rtol=1.e-4, atol=0.)
    # Check sum of atomic components equals total values
    assert_allclose(np.sum(results_iqa['en_atomic']), -1.968747462907e+02, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['kin_atomic']), 7.457892107183e+01, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['col_atomic']), 47.276218, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['ex_total']), -9.100148, rtol=1.e-4, atol=0.)

def test_h2o_rpbepbe_ccpvtz():
    # check total values of decomposition against in h2o_rpbepbe_sto3g.log
    mol_iqa = IQA.from_file('h2o_rpbepbe_sto3g.fchk', grid_type='becke', scheme='H')
    results_iqa = mol_iqa.iqa(dft_exch="gga_x_pbe", dft_corr="gga_c_pbe")
    assert_allclose(results_iqa['nn_total'], 9.1559536481, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['en_total'], -1.968736981462e+02    , rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['kin_total'], 7.459229919735e+01, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['col_total'], -0.343334, rtol=1.e-4, atol=0.)
    assert_allclose(results_iqa['ex_total'], -9.018861, rtol=1.e-4, atol=0.)
    # Check sum of atomic components equals total values
    assert_allclose(np.sum(results_iqa['en_atomic']), -1.968736981462e+02, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['kin_atomic']), 7.459229919735e+01, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['col_atomic']), -0.343334, rtol=1.e-4, atol=0.)
    assert_allclose(np.sum(results_iqa['ex_total']), -9.018861, rtol=1.e-4, atol=0.)