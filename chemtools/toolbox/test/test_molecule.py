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
"""Test chemtools.toolbox.molecule."""

from numpy.testing import assert_raises
import numpy as np
from chemtools import context
from chemtools.toolbox import molecule
from chemtools.utils.molecule import BaseMolecule
from chemtools.utils.wrappers import HortonMolecule


def test_make_molecule():
    """Test chemtools.toolbox.molecule.make_molecule."""
    assert_raises(NotImplementedError, molecule.make_molecule, package_name='gibberish')
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    test = molecule.make_molecule(file_path)
    assert isinstance(test, BaseMolecule)
    assert isinstance(test, HortonMolecule)
    for key, val in HortonMolecule.from_file(file_path).__dict__.items():
        if key == '_iodata':
            continue
        try:
            assert getattr(test, key) == val
        except ValueError:
            assert np.allclose(getattr(test, key), val)
