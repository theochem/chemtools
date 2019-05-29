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
# pragma pylint: disable=invalid-name,bad-whitespace,bad-continuation
"""Test chemtools.toolbox.dftbased."""


import numpy as np

from numpy.testing import assert_raises, assert_array_almost_equal

from chemtools.wrappers.molecule import Molecule
from chemtools.toolbox.dftbased import DFTBasedTool
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_orbital_based_raises():
    # check file name
    assert_raises(ValueError, DFTBasedTool.from_file, "gibberish", np.array([0.0, 1.0]))
    with path("chemtools.data", "h2o_dimer_pbe_sto3g.wfn") as fname:
        assert_raises(ValueError, DFTBasedTool.from_file, fname, np.array([0.0, 1.0]))
        assert_raises(ValueError, DFTBasedTool.from_file, fname, np.array([0.0, 1.0]))
    with path("chemtools.data", "h2o_dimer_pbe_sto3g.wfn") as fname:
        assert_raises(ValueError, DFTBasedTool.from_file, fname, np.array([0.0, 1.0]))
        assert_raises(ValueError, DFTBasedTool.from_file, fname, np.array([0.0, 1.0, 0.0]))
        assert_raises(ValueError, DFTBasedTool.from_file, fname, np.array([[0, 1], [1, 0.]]))
        # check spin argument
        tool = DFTBasedTool.from_file(fname, np.array([[0., 0., 0.]]))
    assert_raises(KeyError, tool._compute_orbital_expression, np.array([9]), spin="error")
    assert_raises(KeyError, tool._compute_orbital_expression, np.array([9]), spin="alph")
    assert_raises(KeyError, tool._compute_orbital_expression, np.array([9]), spin="bet")


def check_orbital_based_properties(tool, data):
    """Test OrbitalLocalTool against stored data arrays."""
    # check spin chemical potential at T=25000K
    result = tool.compute_spin_chemical_potential(25000.0)
    assert_array_almost_equal(result, data["spin_mu_temp_25000"], decimal=6)
    # check density at T=25000K
    result = tool.compute_temperature_dependent_density(25000.0)
    assert_array_almost_equal(result[0], 0.5 * data["density_temp_25000"], decimal=6)
    assert_array_almost_equal(result[1], 0.5 * data["density_temp_25000"], decimal=6)
    # check local density of state at T=25000K
    result = tool.compute_temperature_dependent_state(25000)
    assert_array_almost_equal(result[0], 0.5 * data["density_states_temp_25000"], decimal=6)
    assert_array_almost_equal(result[1], 0.5 * data["density_states_temp_25000"], decimal=6)
    # check local ionization potential
    result = tool.average_local_ionization_energy
    assert_array_almost_equal(result[0], 0.5 * data["local_ip"], decimal=4)
    assert_array_almost_equal(result[1], 0.5 * data["local_ip"], decimal=4)


def test_orbital_based_from_file_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = DFTBasedTool.from_file(fname, data["points"])
    check_orbital_based_properties(tool, data)


def test_orbital_based_from_molecule_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_molecule initialization with exp_beta & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_based_properties(tool, data)
    # test from_molecule initialization without exp_beta & check against Fortran code
    del molecule._exp_beta
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_based_properties(tool, data)


def test_orbital_based_from_file_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = DFTBasedTool.from_file(fname, data["points"])
    check_orbital_based_properties(tool, data)


def test_orbital_based_from_molecule_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_based_properties(tool, data)


def check_orbital_expression(tool, data):
    """Check OrbitalLocalTool.compute_orbital_expression against stored data array."""
    result = tool._compute_orbital_expression(5)
    assert_array_almost_equal(result[:, 0], data["orbital_05"], decimal=6)
    result = tool._compute_orbital_expression(np.array([5]), spin="beta")
    assert_array_almost_equal(result[:, 0], data["orbital_05"], decimal=6)
    result = tool._compute_orbital_expression(6)
    assert_array_almost_equal(result[:, 0], data["orbital_06"], decimal=6)
    result = tool._compute_orbital_expression(8)
    assert_array_almost_equal(result[:, 0], data["orbital_08"], decimal=6)
    result = tool._compute_orbital_expression(9)
    assert_array_almost_equal(result[:, 0], data["orbital_09"], decimal=6)
    result = tool._compute_orbital_expression([8, 9])
    assert_array_almost_equal(result[:, 0], data["orbital_08"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_09"], decimal=6)
    result = tool._compute_orbital_expression([5, 8, 9], spin="beta")
    assert_array_almost_equal(result[:, 0], data["orbital_05"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_08"], decimal=6)
    assert_array_almost_equal(result[:, 2], data["orbital_09"], decimal=6)
    result = tool._compute_orbital_expression((8, 9, 5))
    assert_array_almost_equal(result[:, 0], data["orbital_08"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_09"], decimal=6)
    assert_array_almost_equal(result[:, 2], data["orbital_05"], decimal=6)
    result = tool._compute_orbital_expression(np.array([8, 9]))
    assert_array_almost_equal(result[:, 0], data["orbital_08"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_09"], decimal=6)
    result = tool._compute_orbital_expression(np.array([6, 5]), spin="b")
    assert_array_almost_equal(result[:, 0], data["orbital_06"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_05"], decimal=6)
    result = tool._compute_orbital_expression(np.array([5, 6, 8, 9]), spin="a")
    assert_array_almost_equal(result[:, 0], data["orbital_05"], decimal=6)
    assert_array_almost_equal(result[:, 1], data["orbital_06"], decimal=6)
    assert_array_almost_equal(result[:, 2], data["orbital_08"], decimal=6)
    assert_array_almost_equal(result[:, 3], data["orbital_09"], decimal=6)


def test_orbital_based_from_file_orbital_expression_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = DFTBasedTool.from_file(fname, data["points"])
    check_orbital_expression(tool, data)


def test_orbital_based_from_molecule_orbital_expression_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_molecule initialization with exp_beta & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_expression(tool, data)
    # test from_molecule initialization without exp_beta & check against Fortran code
    del molecule._exp_beta
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_expression(tool, data)


def test_orbital_based_from_file_orbital_expression_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = DFTBasedTool.from_file(fname, data["points"])
    check_orbital_expression(tool, data)


def test_orbital_based_from_molecule_orbital_expression_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_orbitalbased_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = DFTBasedTool(molecule, data["points"])
    check_orbital_expression(tool, data)


def test_orbital_based_h2o_b3lyp_sto3g():
    # points array
    points = np.array([[-3.,-3.,-3.], [-3.,-3., 0.], [-3.,-3., 3.], [-3., 0.,-3.], [-3., 0., 0.],
                       [-3., 0., 3.], [-3., 3.,-3.], [-3., 3., 0.], [-3., 3., 3.], [ 0.,-3.,-3.],
                       [ 0.,-3., 0.], [ 0.,-3., 3.], [ 0., 0.,-3.], [ 0., 0., 0.], [ 0., 0., 3.],
                       [ 0., 3.,-3.], [ 0., 3., 0.], [ 0., 3., 3.], [ 3.,-3.,-3.], [ 3.,-3., 0.],
                       [ 3.,-3., 3.], [ 3., 0.,-3.], [ 3., 0., 0.], [ 3., 0., 3.], [ 3., 3.,-3.],
                       [ 3., 3., 0.], [ 3., 3., 3.]])
    # initialize OrbitalLocalTool from_file
    with path("chemtools.data", "water_b3lyp_sto3g.fchk") as fname:
        tool = DFTBasedTool.from_file(fname, points)
    # check mep against Fortran code
    expected = np.array([-0.01239766, -0.02982537, -0.02201149,   -0.01787292, -0.05682143,
                         -0.02503563, -0.00405942, -0.00818772,   -0.00502268,  0.00321181,
                         -0.03320573, -0.02788605,  0.02741914, 1290.21135500, -0.03319778,
                          0.01428660,  0.10127092,  0.01518299,    0.01530548,  0.00197975,
                         -0.00894206,  0.04330806,  0.03441681,   -0.00203017,  0.02272626,
                          0.03730846,  0.01463959])
    assert_array_almost_equal(tool.electrostatic_potential, expected, decimal=6)
    # check spin chemical potential against manual calculation
    result = tool.compute_spin_chemical_potential(25000.0)
    assert_array_almost_equal(result, [0.10821228040, 0.10821228040], decimal=6)
