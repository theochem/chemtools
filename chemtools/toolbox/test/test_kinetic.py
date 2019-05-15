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
# pragma pylint: disable=invalid-name,fixme
"""Test chemtools.toolbox.kinetic."""


import numpy as np

from numpy.testing import assert_array_almost_equal, assert_allclose

from chemtools.wrappers.molecule import Molecule
from chemtools.toolbox.kinetic import KED
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def check_kinetic_energy_density(tool, data):
    """Test KineticEnergyDensity against stored data arrays."""
    # check density & gradient
    assert_array_almost_equal(tool.density, data['density'], decimal=6)
    # check kinetic energy density
    assert_array_almost_equal(tool.ked_positive_definite, data['ke_positive_definite'], decimal=6)
    assert_array_almost_equal(tool.ked_weizsacker, data['ke_weizsacker'], decimal=6)
    # TODO: Why assert_array_almost_equal doesn't work for ke_thomas_fermi
    assert_allclose(tool.ked_thomas_fermi, data['ke_thomas_fermi'], rtol=1e-08, atol=1e-08)


def test_kinetic_from_file_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = KED.from_file(fname, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_molecule_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_molecule initialization with exp_beta & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['points'])
    print('tool = ', tool)
    check_kinetic_energy_density(tool, data)
    # test from_molecule initialization without exp_beta & check against Fortran code
    del molecule._exp_beta
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_file_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = KED.from_file(fname, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_molecule_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['points'])
    check_kinetic_energy_density(tool, data)
