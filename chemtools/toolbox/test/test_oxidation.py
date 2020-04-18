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
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def _get_eos(filename, scheme):
    # build proatom database
    atoms = glob.glob('chemtools/data/atom_0*')
    proatomdb = ProAtomDB.from_files(atoms, "power:5e-8:20:40:146")
    # load molecule & make grid, denspart, and eos instances
    with path('chemtools.data', filename) as file_path:
        mol = Molecule.from_file(file_path)
    grid = MolecularGrid.from_molecule(mol, specs='power:5e-8:20:40:146', k=4, rotate=False)
    part = DensPart.from_molecule(mol, scheme=scheme, grid=grid, local=False, proatomdb=proatomdb)
    return EOS.from_molecule(mol, part, grid)


def test_eos_h_h2o_3fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    occs_o = np.array([0.9948, 0.8953, 0.7958, 0.5573, 0.5348])
    occs_h = np.array([0.1858, 0.0194, 0.0126, 0.0030, 0.0])
    # test occupations default fragments
    result = eos.compute_fragment_occupation(spin='a')
    assert_almost_equal(occs_o, result[0], decimal=3)
    assert_almost_equal(occs_h, result[1], decimal=2)
    assert_almost_equal(occs_h, result[2], decimal=2)
    # test oxidation states default fragments
    assert_equal([-2.0, 1.0, 1.0], eos.compute_oxidation_state())
    assert_approx_equal(84.895, eos.reliability, significant=3)

    # test occupations given 3 fragments
    result = eos.compute_fragment_occupation([[0], [1], [2]], spin='a')
    assert_almost_equal(occs_o, result[0], decimal=2)
    assert_almost_equal(occs_h, result[1], decimal=2)
    assert_almost_equal(occs_h, result[2], decimal=2)
    # test oxidation states given 3 fragments
    assert_equal([-2.0, 1.0, 1.0], eos.compute_oxidation_state([[0], [1], [2]]))
    assert_approx_equal(84.895, eos.reliability, significant=3)


def test_eos_h_h2o_2fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    occs_f1 = np.array([0.9974, 0.9678, 0.9194, 0.8942, 0.5931])
    occs_f2 = np.array([0.1859, 0.0194, 0.0126, 0.0030, 0.000])
    # test occupations
    result = eos.compute_fragment_occupation([[0, 1], [2]], spin='a')
    assert_almost_equal(occs_f1, result[0], decimal=3)
    assert_almost_equal(occs_f2, result[1], decimal=2)
    # test oxidation states
    assert_equal([-1.0, 1.0], eos.compute_oxidation_state([[0, 1], [2]]))
    assert_approx_equal(90.726, eos.reliability, significant=3)


def test_eos_h_h2o_1fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    occs = np.array([1.0000, 1.0000, 1.0000, 1.0000, 0.9999])
    # test occupations
    result = eos.compute_fragment_occupation([[0, 1, 2]], spin='a')
    assert_almost_equal(occs, result[0], decimal=3)
    result = eos.compute_fragment_occupation([[0, 1, 2]], spin='b')
    assert_almost_equal(occs, result[0], decimal=3)
    # test oxidation states
    assert_equal([0.0], eos.compute_oxidation_state([[0, 1, 2]]))
    assert_approx_equal(100.0, eos.reliability, significant=7)
