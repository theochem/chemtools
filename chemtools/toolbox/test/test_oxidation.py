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


import glob
import numpy as np
from numpy.testing import assert_equal, assert_almost_equal, assert_approx_equal

from horton import ProAtomDB
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart
from chemtools.toolbox.oxidation import EOS


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
    return EOS(mol, part)


def check_oxidation_states(eos, occs, fragments, spin, oxidations, reliability, decimal=3):
    # test occupations
    result = eos.compute_fragment_occupation(fragments, spin=spin)
    assert occs.shape == result.shape
    assert_almost_equal(occs, result, decimal=decimal)
    # test oxidation states & reliability
    assert_equal(oxidations, eos.compute_oxidation_state(fragments))
    assert_approx_equal(reliability, eos.reliability, significant=3)


def test_eos_h_h2o_3fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    charges = np.array([-0.302250, 0.151138, 0.151185])
    occs = np.array([[0.9948, 0.8953, 0.7958, 0.5573, 0.5348],
                     [0.1858, 0.0194, 0.0126, 0.0030, 0.0],
                     [0.1858, 0.0194, 0.0126, 0.0030, 0.0]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=3)
    # test occupations for alpha & beta orbitals using default fragments
    check_oxidation_states(eos, occs, None, 'a', [-2.0, 1.0, 1.0], 84.895, decimal=3)
    check_oxidation_states(eos, occs, None, 'b', [-2.0, 1.0, 1.0], 84.895, decimal=3)
    # test occupations for alpha & beta orbitals given 3 fragments
    check_oxidation_states(eos, occs, [[0], [1], [2]], 'a', [-2.0, 1.0, 1.0], 84.895, decimal=3)
    check_oxidation_states(eos, occs, [[0], [1], [2]], 'b', [-2.0, 1.0, 1.0], 84.895, decimal=3)


def test_eos_h_h2o_2fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    charges = np.array([-0.302250, 0.151138, 0.151185])
    occs = np.array([[0.9974, 0.9678, 0.9194, 0.8942, 0.5931],
                     [0.1859, 0.0194, 0.0126, 0.0030, 0.000]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=3)
    # test occupations & oxidation states for alpha & beta orbitals
    check_oxidation_states(eos, occs, [[0, 1], [2]], 'a', [-1.0, 1.0], 90.726, decimal=3)
    check_oxidation_states(eos, occs, [[0, 1], [2]], 'b', [-1.0, 1.0], 90.726, decimal=3)


def test_eos_h_h2o_1fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'h')
    charges = np.array([-0.302250, 0.151138, 0.151185])
    occs = np.array([[1.0000, 1.0000, 1.0000, 1.0000, 0.9999]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=3)
    # test occupations & oxidation states for alpha & beta orbitals
    check_oxidation_states(eos, occs, [[0, 1, 2]], 'a', [0.0], 100.0, decimal=3)
    check_oxidation_states(eos, occs, [[0, 1, 2]], 'b', [0.0], 100.0, decimal=3)


def test_eos_hi_h2o_3fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'hi')
    charges = np.array([-0.90460, 0.45225, 0.45228])
    occs = np.array([[0.9972, 0.9431, 0.8837, 0.6757, 0.6483],
                     [0.1023, 0.0, 0.0, 0.0, 0.0], [0.1023, 0.0, 0.0, 0.0, 0.0]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=2)
    # test occupations for alpha & beta orbitals using default fragments
    check_oxidation_states(eos, occs, None, 'a', [-2.0, 1.0, 1.0], 100.0, decimal=2)
    check_oxidation_states(eos, occs, None, 'b', [-2.0, 1.0, 1.0], 100.0, decimal=2)
    # test occupations for alpha & beta orbitals given 3 fragments
    check_oxidation_states(eos, occs, [[0], [1], [2]], 'a', [-2.0, 1.0, 1.0], 100.0, decimal=2)
    check_oxidation_states(eos, occs, [[0], [1], [2]], 'b', [-2.0, 1.0, 1.0], 100.0, decimal=2)


def test_eos_hi_h2o_2fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'hi')
    charges = np.array([-0.90460, 0.45225, 0.45228])
    occs = np.array([[0.9986, 0.9813, 0.9535, 0.9408, 0.6933], [0.1023, 0.0, 0.0, 0.0, 0.0]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=2)
    # test occupations & oxidation states for alpha & beta orbitals
    check_oxidation_states(eos, occs, [[0, 1], [2]], 'a', [-1.0, 1.0], 100.0, decimal=2)
    check_oxidation_states(eos, occs, [[0, 1], [2]], 'b', [-1.0, 1.0], 100.0, decimal=2)


def test_eos_hi_h2o_1fragments():
    # test against APOST-3D (version 3.1)
    eos = _get_eos('h2o_q+0_ub3lyp_ccpvtz.fchk', 'hi')
    charges = np.array([-0.90460, 0.45225, 0.45228])
    occs = np.array([[1.0000, 1.0000, 1.0000, 1.0000, 0.9999]])
    # test atomic charges
    assert_almost_equal(charges, eos.part.charges, decimal=2)
    # test occupations & oxidation states for alpha & beta orbitals
    check_oxidation_states(eos, occs, [[0, 1, 2]], 'a', [0.0], 100.0, decimal=3)
    check_oxidation_states(eos, occs, [[0, 1, 2]], 'b', [0.0], 100.0, decimal=3)
