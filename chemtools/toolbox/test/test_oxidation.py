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


from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.part import DensPart
from chemtools.toolbox.oxidation import EOS
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_eos_h2o():
    # test against Apost
    with path('chemtools.data.examples', 'h2o.fchk') as file_path:
        mol = Molecule.from_file(file_path)
        part = DensPart.from_molecule(mol, scheme="h")
        eos = EOS.from_molecule(mol, part, part.grid)

    # test oxidation state of each atom and then fragments
    # use eos.compute_oxidation_state(fragments)
