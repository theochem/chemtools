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
"""Test chemtools.toolbox.oxidation."""


import numpy as np

from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart
from chemtools.toolbox.oxidation import EOS
from horton import ProAtomDB

import glob as glob
from numpy.testing import assert_equal, assert_almost_equal

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_eos_h_h2o_3fragments():
    # test against APOST-3D (version 3.1)
    atoms = glob.glob('chemtools/data/atom_0*')
    proatomdb = ProAtomDB.from_files(atoms, "power:5e-8:20:40:146")

    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        # make proatom & pass it to Denspart
        grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                             specs='power:5e-8:20:40:146', rotate=False, k=4)
        part = DensPart.from_molecule(mol, scheme="h", grid=grid, local=False, proatomdb=proatomdb)
        eos = EOS.from_molecule(mol, part, grid)
    # test occupations default fragments
    occupation_o = np.array([0.9948, 0.8953, 0.7958, 0.5573, 0.5348])
    occupation_h = np.array([0.1858, 0.0194, 0.0126, 0.0030, 0.0])
    result = eos.compute_fragment_occupation(spin='a')
    assert_almost_equal(occupation_o, result[0], decimal=3)
    assert_almost_equal(occupation_h, result[1], decimal=2)
    assert_almost_equal(occupation_h, result[2], decimal=2)
    # test oxidation states default fragments
    assert_equal([-2.0, 1.0, 1.0],  eos.compute_oxidation_state())

    # test occupations given 3 fragemts
    result = eos.compute_fragment_occupation([[0], [1], [2]], spin='a')
    assert_almost_equal(occupation_o, result[0], decimal=2)
    assert_almost_equal(occupation_h, result[1], decimal=2)
    assert_almost_equal(occupation_h, result[2], decimal=2)
    # test oxidation states given 3 fragments
    assert_equal([-2.0, 1.0, 1.0],  eos.compute_oxidation_state([[0], [1], [2]]))


def test_eos_h_h2o_2fragments():
    # test against APOST-3D (version 3.1)
    atoms = glob.glob('chemtools/data/atom_0*')
    proatomdb = ProAtomDB.from_files(atoms, "power:5e-8:20:40:146")

    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        # make proatom & pass it to Denspart
        grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                             specs='power:5e-8:20:40:146', rotate=False, k=4)
        part = DensPart.from_molecule(mol, scheme="h", grid=grid, local=False, proatomdb=proatomdb)
        eos = EOS.from_molecule(mol, part, grid)
    # test occupations
    occs_f1 = np.array([0.9974,   0.9678,   0.9194,   0.8942,   0.5931])
    occs_f2 = np.array([0.1859, 0.0194, 0.0126, 0.0030, 0.000])
    result = eos.compute_fragment_occupation([[0, 1], [2]], spin='a')
    assert_almost_equal(occs_f1, result[0], decimal=3)
    assert_almost_equal(occs_f2, result[1], decimal=2)
    # test oxidation states
    assert_equal([-1.0, 1.0],  eos.compute_oxidation_state([[0, 1], [2]]))


def test_eos_h_h2o_1fragments():
    # test against APOST-3D (version 3.1)
    atoms = glob.glob('chemtools/data/atom_0*')
    proatomdb = ProAtomDB.from_files(atoms, "power:5e-8:20:40:146")

    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                             specs='power:5e-8:20:40:146', rotate=False, k=4)
        part = DensPart.from_molecule(mol, scheme="h", grid=grid, local=False,  proatomdb=proatomdb)
        eos = EOS.from_molecule(mol, part, grid)
    # test occupations
#    occs = np.array([1.0, 1.0, 1.0, 0.9, 0.9])
    occs = np.array([1.0000,   1.0000,   1.0000,   1.0000,   0.9999])
    result = eos.compute_fragment_occupation([[0, 1, 2]], spin='a')
    assert_almost_equal(occs, result[0], decimal=3)
    result = eos.compute_fragment_occupation([[0, 1, 2]], spin='b')
    assert_almost_equal(occs, result[0], decimal=3)
    # test oxidation states
    assert_equal([0.0], eos.compute_oxidation_state([[0, 1, 2]]))
