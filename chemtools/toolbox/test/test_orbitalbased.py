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
# pragma pylint: disable=invalid-name
"""Test chemtools.toolbox.orbitalbased."""


import numpy as np

from numpy.testing import assert_raises, assert_array_almost_equal, assert_allclose

from chemtools.utils.utils import context
from chemtools.wrappers.molecule import Molecule
from chemtools.toolbox.orbitalbased import OrbitalLocalTool


def test_orbital_based_raises():
    # check file name
    assert_raises(ValueError, OrbitalLocalTool.from_file, "gibberish", np.array([0.0, 1.0]))
    filename = context.get_fn("test/h2o_dimer_pbe_sto3g.wf")
    assert_raises(ValueError, OrbitalLocalTool.from_file, filename, np.array([0.0, 1.0]))
    filename = context.get_fn("test/h2o_dimer_pbe_sto3g.wfn")
    assert_raises(ValueError, OrbitalLocalTool.from_file, filename, np.array([0.0, 1.0]))
    assert_raises(ValueError, OrbitalLocalTool.from_file, filename, np.array([0.0, 1.0, 0.0]))
    assert_raises(ValueError, OrbitalLocalTool.from_file, filename, np.array([[0, 1], [1, 0.]]))


def check_orbital_based(tool, data, index_homo, index_lumo):
    """Test OrbitalLocalTool against stored data arrays."""
    # check density & gradient
    assert_array_almost_equal(tool.density, data["density"], decimal=6)
    assert_array_almost_equal(tool.gradient, data["gradient"], decimal=6)
    # check kinetic energy density
    result = tool.weizsacker_kinetic_energy_density
    assert_array_almost_equal(result, data["ke_weizsacker"], decimal=6)
    result = tool.thomas_fermi_kinetic_energy_density
    assert_allclose(result, data["ke_thomas_fermi"], rtol=1e-08, atol=1e-08)
    result = tool.kinetic_energy_density
    assert_array_almost_equal(result, data["ke_positive_definite"], decimal=6)
    # check HOMO & LUMO density
    result = tool.compute_orbital_expression(index_homo)
    assert_array_almost_equal(result[:, 0], data["density_homo"], decimal=6)
    result = tool.compute_orbital_expression(index_lumo)
    assert_array_almost_equal(result[:, 0], data["density_lumo"], decimal=6)
    result = tool.compute_orbital_expression([index_homo, index_lumo])
    assert_array_almost_equal(result[:, 0], data["density_homo"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["density_lumo"], decimal=6)
    result = tool.compute_orbital_expression((index_homo, index_lumo))
    assert_array_almost_equal(result[:, 0], data["density_homo"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["density_lumo"], decimal=6)
    result = tool.compute_orbital_expression(np.array([index_homo, index_lumo]))
    assert_array_almost_equal(result[:, 0], data["density_homo"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["density_lumo"], decimal=6)


def test_orbital_based_from_file_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    data = np.load(context.get_fn("test/data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz"))
    # test from_file initialization
    tool = OrbitalLocalTool.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"), data["points"])
    # check against Fortran code
    check_orbital_based(tool, data, 8, 9)


def test_orbital_based_from_molecule_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    data = np.load(context.get_fn("test/data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz"))
    # test from_molecule initialization with exp_beta
    molecule = Molecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    tool = OrbitalLocalTool(molecule, data["points"])
    # check against Fortran code
    check_orbital_based(tool, data, 8, 9)
    # test from_molecule initialization without exp_beta
    del molecule._exp_beta
    molecule = Molecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    tool = OrbitalLocalTool(molecule, data["points"])
    # check against Fortran code
    check_orbital_based(tool, data, 8, 9)


def test_orbital_based_from_file_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    data = np.load(context.get_fn("test/data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz"))
    # test from_file initialization
    tool = OrbitalLocalTool.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"), data["points"])
    # check against Fortran code
    check_orbital_based(tool, data, 8, 9)


def test_orbital_based_from_molecule_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    data = np.load(context.get_fn("test/data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz"))
    # test from_file initialization
    molecule = Molecule.from_file(context.get_fn("test/ch4_rhf_ccpvdz.fchk"))
    tool = OrbitalLocalTool(molecule, data["points"])
    # check against Fortran code
    check_orbital_based(tool, data, 8, 9)
