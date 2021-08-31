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


import sys
import numpy as np

from numpy.testing import assert_array_almost_equal, assert_allclose, assert_raises

from chemtools.toolbox.kinetic import KED

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule
else:
    from chemtools.wrappers3.molecule import Molecule

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def check_kinetic_energy_density(tool, data):
    """Test KineticEnergyDensity against stored data arrays."""
    # check density & gradient
    assert_array_almost_equal(tool.density, data["dens"], decimal=6)
    # check kinetic energy density
    assert_array_almost_equal(tool.ked_positive_definite, data["ked_pd"], decimal=6)
    assert_array_almost_equal(tool.ked_weizsacker, data["ked_w"], decimal=6)
    # TODO: Why assert_array_almost_equal doesn't work for ke_thomas_fermi
    assert_allclose(tool.ked_thomas_fermi, data["ked_tf"], rtol=1e-08, atol=1e-08)


def test_kinetic_from_file_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = KED.from_file(fname, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_molecule_ch4_uhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_molecule initialization with exp_beta & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_file_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        tool = KED.from_file(fname, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_from_molecule_ch4_rhf_ccpvdz():
    # load data computed with Fortran code
    with path("chemtools.data", "data_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    # test from_file initialization & check against Fortran code
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['points'])
    check_kinetic_energy_density(tool, data)


def test_kinetic_raises_lap_ked_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    tool = KED(data['nuc_dens'], data['nuc_grad'], None, None)
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_raises(ValueError, lambda: tool.laplacian)
    assert_raises(ValueError, lambda: tool.ked_positive_definite)
    assert_raises(ValueError, lambda: tool.ked_gradient_expansion)
    assert_raises(ValueError, lambda: tool.ked_gradient_expansion_empirical)
    assert_raises(ValueError, tool.ked_gradient_expansion_general, 1.0, 1.0)
    assert_raises(ValueError, lambda: tool.ked_hamiltonian)
    assert_raises(ValueError, tool.ked_general, 1.0)


def test_kinetic_raises_lap_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    tool = KED(data['nuc_dens'], data['nuc_grad'], None, data['nuc_ked_pd'])
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_raises(ValueError, lambda: tool.laplacian)
    assert_raises(ValueError, lambda: tool.ked_gradient_expansion)
    assert_raises(ValueError, lambda: tool.ked_gradient_expansion_empirical)
    assert_raises(ValueError, tool.ked_gradient_expansion_general, 1.0, 1.0)
    assert_raises(ValueError, lambda: tool.ked_hamiltonian)
    assert_raises(ValueError, tool.ked_general, 1.0)


def test_kinetic_raises_ked_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    tool = KED(data['nuc_dens'], data['nuc_grad'], data['nuc_lap'], None)
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_raises(ValueError, lambda: tool.ked_positive_definite)
    assert_raises(ValueError, lambda: tool.ked_hamiltonian)
    assert_raises(ValueError, tool.ked_general, 1.0)


def test_kinetic_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    tool = KED(data['nuc_dens'], data['nuc_grad'], data['nuc_lap'], data['nuc_ked_pd'])
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.laplacian, data['nuc_lap'], rtol=1.e-6, atol=0.)
    # check ked at the position of nuclei
    assert_allclose(tool.ked_positive_definite, data['nuc_ked_pd'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_hamiltonian, data['nuc_ked_ham'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_general(alpha=0.), data['nuc_ked_ham'], rtol=1.e-5, atol=0.)


def test_kinetic_from_file_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        tool = KED.from_file(fname, data['coords'])
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.laplacian, data['nuc_lap'], rtol=1.e-6, atol=0.)
    # check ked at the position of nuclei
    assert_allclose(tool.ked_positive_definite, data['nuc_ked_pd'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_hamiltonian, data['nuc_ked_ham'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_general(alpha=0.), data['nuc_ked_ham'], rtol=1.e-5, atol=0.)


def test_kinetic_from_molecule_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    tool = KED.from_molecule(molecule, data['coords'])
    # check attributes
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-6, atol=0.)
    assert_allclose(tool.laplacian, data['nuc_lap'], rtol=1.e-6, atol=0.)
    # check ked at the position of nuclei
    assert_allclose(tool.ked_positive_definite, data['nuc_ked_pd'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_hamiltonian, data['nuc_ked_ham'], rtol=1.e-5, atol=0.)
    assert_allclose(tool.ked_general(alpha=0.), data['nuc_ked_ham'], rtol=1.e-5, atol=0.)
