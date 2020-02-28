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

from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart
from chemtools.toolbox.oxidation import EOS

from numpy.testing import assert_equal, assert_almost_equal

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_eos_h_h2o_3fragments():
    # test against APOST-3D (version 3.1)
    with path('chemtools.data.examples', 'h2o.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        # make proatom & pass it to Denspart
        part = DensPart.from_molecule(mol, scheme="h")
        eos = EOS.from_molecule(mol, part, part.grid)
    # test occupations
    occupation_o = np.array([0.9950, 0.8933, 0.7876, 0.5591, 0.5384])
    occupation_h = np.array([0.1860, 0.0219, 0.0135, 0.0, 0.0])
    result = eos.compute_fragment_occupation()
    # assert_almost_equal(occupation_o, result[0], decimal=3)
    # assert_almost_equal(occupation_h, result[1], decimal=3)
    # assert_almost_equal(occupation_h, result[2], decimal=3)
    # test oxidation states
    # assert_equal(eos.compute_oxidation_state(), [-2.0, 1.0, 1.0])


def test_eos_h_h2o_2fragments():
    # test against APOST-3D (version 3.1)
    with path('chemtools.data.examples', 'h2o.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        # make proatom & pass it to Denspart
        part = DensPart.from_molecule(mol, scheme="h")
        eos = EOS.from_molecule(mol, part, part.grid)
    # test occupations
    occs_f1 = np.array([0.9975, 0.9680, 0.9184, 0.8895, 0.5958])
    occs_f2 = np.array([0.1860, 0.0219, 0.0135])
    result = eos.compute_fragment_occupation()
    # assert_almost_equal(occs_f1, result[0], decimal=3)
    # assert_almost_equal(occs_f2, result[1], decimal=3)
    # test oxidation states
    # assert_equal(eos.compute_oxidation_state(), [-1.0, 1.0])


def test_eos_h_h2o_1fragments():
    # test against APOST-3D (version 3.1)
    with path('chemtools.data.examples', 'h2o.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        # make proatom & pass it to Denspart
        part = DensPart.from_molecule(mol, scheme="h")
        eos = EOS.from_molecule(mol, part, part.grid)
    # test occupations
    occs = np.array([1.0000, 1.0000, 1.0000, 1.0000, 1.0000])
    # assert_almost_equal(occs, eos.compute_fragment_occupation([[0, 1, 2]]), decimal=3)
    # test oxidation states
    # assert_equal(eos.compute_oxidation_state(), [0.0])
