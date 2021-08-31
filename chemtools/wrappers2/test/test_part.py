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

import numpy as np

from chemtools.wrappers2.molecule import Molecule
from chemtools.wrappers2.part import DensPart
from chemtools.wrappers2.grid import MolecularGrid
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_condense_linear_from_file_fmr_h_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        part = DensPart.from_file(fname, scheme='h')
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    computed = part.numbers - part.charges
    assert np.all(abs(expected - computed) < 1.e-3)
    assert np.all(abs(part.condense_to_atoms(part.density) - computed) < 1.e-2)


def test_condense_quadratic_from_molecule_fmr_mbis_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        molecule = Molecule.from_file(fname)
    part = DensPart.from_molecule(molecule, scheme='mbis')
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    computed = part.numbers - part.charges
    assert np.all(abs(expected - computed) < 1.e-2)
    assert np.all(abs(part.condense_to_atoms(part.density) - computed) < 1.e-2)


def test_condense_quadratic_from_molecule_fmr_mbis_ch4():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        molecule = Molecule.from_file(fname)
        grid = MolecularGrid.from_molecule(molecule, 'fine')
    part = DensPart.from_molecule(molecule, scheme='mbis', grid=grid)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    computed = part.numbers - part.charges
    assert np.all(abs(expected - computed) < 1.e-2)
    assert np.all(abs(part.condense_to_atoms(part.density) - computed) < 1.e-2)
