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
# pragma pylint: disable=invalid-name
"""Test chemtools.analysis.conceptual.GlobalConceptualDFT."""


import sys
import numpy as np

from numpy.testing import assert_raises, assert_almost_equal

from chemtools.toolbox.conceptual import GlobalConceptualDFT

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule
else:
    from chemtools.wrappers3.molecule import Molecule

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_global_conceptual_raises():
    # dictionary of energy values
    values = {0.0: 0.0, 1.0: -0.5, 2.0: -0.45}
    # check invalid model
    assert_raises(ValueError, GlobalConceptualDFT, values, "quad")
    assert_raises(ValueError, GlobalConceptualDFT, values, "Gibberish")
    assert_raises(NotImplementedError, GlobalConceptualDFT, values, "general")
    # check invalid coordinates
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", np.array([0., 0., 0.]))
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", np.array([[0., 0.]]))
    # check invalid atomic numbers
    coord = np.array([[0., 0., 0.]])
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", coord, [])
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", coord, [1, 1])
    # check invalid number of electrons
    values = {"0": 0.0, 1: -0.5, 2: -0.45}
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", coord)
    assert_raises(ValueError, GlobalConceptualDFT, values, "quadratic", coord)
    values = {0.5: 0.0, 1.3: -0.5, 2: -0.45}
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", coord)
    assert_raises(ValueError, GlobalConceptualDFT, values, "quadratic", coord)
    values = {0: 0.0, 1: -0.5, 3: -0.45}
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear")
    assert_raises(ValueError, GlobalConceptualDFT, values, "quadratic")
    values = {0: 0.0, 1: -0.5, 2: -0.45, 3: -0.4}
    assert_raises(ValueError, GlobalConceptualDFT, values, "linear", coord)
    assert_raises(ValueError, GlobalConceptualDFT, values, "quadratic", coord)
    # check non-existing attribute
    model = GlobalConceptualDFT({0.0: 0.0, 1.0: -0.5, 2.0: -0.45}, "quadratic")
    assert_raises(AttributeError, getattr, model, "mu_plus")
    assert_raises(AttributeError, getattr, model, "gibberish")
    # check molecule file inconsistency
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as file1:
        with path('chemtools.data', 'o2_uhf.fchk') as file2:
            fnames = [file1, file2]
            assert_raises(ValueError, GlobalConceptualDFT.from_file, fnames, "linear")
            assert_raises(ValueError, GlobalConceptualDFT.from_file, fnames, "quadratic")
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        assert_raises(ValueError, GlobalConceptualDFT.from_file, [fname, fname], "linear")
        assert_raises(ValueError, GlobalConceptualDFT.from_file, [fname, fname], "quadratic")
        assert_raises(ValueError, GlobalConceptualDFT.from_file, [fname, fname, fname], "linear")
        assert_raises(ValueError, GlobalConceptualDFT.from_file, [fname, fname, fname],
                      "quadratic")


def check_global_reactivity_linear(model, ip, ea, energy, n):
    """Check expected linear global reactivity descriptors."""
    # check print statement
    assert isinstance(model.__repr__(), str)
    # check energy values
    assert_almost_equal(model.energy(n), energy, decimal=6)
    assert_almost_equal(model.energy(n - 1), energy + ip, decimal=6)
    assert_almost_equal(model.energy(n + 1), energy - ea, decimal=6)
    # check ionization potential and electron affinity
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    # check N0
    assert_almost_equal(model.n0, n, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    assert model.mu is None
    assert model.chemical_potential is None
    assert model.eta is None
    assert model.chemical_hardness is None
    assert model.hyper_hardness(2) is None
    assert model.hyper_hardness(3) is None
    assert model.hyper_hardness(4) is None
    # check mu+, mu-, mu0
    assert_almost_equal(model.mu_plus, -ea, decimal=6)
    assert_almost_equal(model.mu_minus, -ip, decimal=6)
    assert_almost_equal(model.mu_zero, -0.5 * (ip + ea), decimal=6)
    # check derivatives of energy w.r.t. number of electrons
    assert model.energy_derivative(n, 1) is None
    assert model.energy_derivative(n, 2) is None
    assert model.energy_derivative(n, 3) is None
    assert_almost_equal(model.energy_derivative(n + 1, 1), -ea, decimal=6)
    assert_almost_equal(model.energy_derivative(n - 1, 1), -ip, decimal=6)
    assert_almost_equal(model.energy_derivative(n - 0.3, 1), -ip, decimal=6)
    assert_almost_equal(model.energy_derivative(n + 0.5, 1), -ea, decimal=6)
    assert_almost_equal(model.energy_derivative(n + 1, 4), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(n - 1, 3), 0.0, decimal=6)
    assert_almost_equal(model.energy_derivative(n + 0.4, 2), 0.0, decimal=6)
    # check electrophilicity, nucleofugality & electrofugality
    if model.n_max is not None:
        assert_almost_equal(model.electrophilicity, 0., decimal=6)
        assert_almost_equal(model.nucleofugality, -ea, decimal=6)
        assert_almost_equal(model.electrofugality, ip, decimal=6)


def test_global_linear_from_file_fmo_ch4_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # check from_file
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        model = GlobalConceptualDFT.from_file(fname, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)
    # check from_file given as a list
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        model = GlobalConceptualDFT.from_file([fname], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)


def test_global_linear_from_molecule_fmo_ch4_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)


def test_global_linear_from_file_fmo_ch4_wfn():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "linear")
        model2 = GlobalConceptualDFT.from_file([fname], "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 10)
    # check from_file given as a list
    check_global_reactivity_linear(model2, ip, ea, energy, 10)


def test_global_linear_from_molecule_fmo_ch4_wfn():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)


def test_global_linear_from_file_fmo_h2o_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-3.09871604E-01), -2.48704636E-02, -7.645980351270224E+01
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "linear")
        model2 = GlobalConceptualDFT.from_file([fname], "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 10)
    # check from_file given as a list
    check_global_reactivity_linear(model2, ip, ea, energy, 10)


def test_global_linear_from_molecule_fmo_h2o_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-3.09871604E-01), -2.48704636E-02, -7.645980351270224E+01
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 10)


def test_global_linear_from_file_fmo_h2o_cation_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-8.47044131E-01), -(-6.19391831E-01), -7.599493522312368E+01
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "linear")
        model2 = GlobalConceptualDFT.from_file([fname], "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 9)
    # check from_file given as a list
    check_global_reactivity_linear(model2, ip, ea, energy, 9)


def test_global_linear_from_molecule_fmo_h2o_cation_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-8.47044131E-01), -(-6.19391831E-01), -7.599493522312368E+01
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 9)
    # check from_molecule given as a file
    model = GlobalConceptualDFT.from_molecule([molecule], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 9)


def test_global_linear_from_file_fmo_h2o_anion_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -1.93118022E-01, -2.69116912E-01, -7.635212549312298E+01
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "linear")
        model2 = GlobalConceptualDFT.from_file([fname], "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 11)
    # check from_file given as a list
    check_global_reactivity_linear(model2, ip, ea, energy, 11)


def test_global_linear_from_molecule_fmo_h2o_anion_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -1.93118022E-01, -2.69116912E-01, -7.635212549312298E+01
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 11)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "linear")
    check_global_reactivity_linear(model, ip, ea, energy, 11)


def test_global_linear_fd_h2o_fchk():
    ep, e0, en = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea, energy = ep - e0, e0 - en, e0
    # check linear global conceptual DFT model
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file3:
                fname = [file1, file2, file3]
                model1 = GlobalConceptualDFT.from_file(fname, "linear")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_linear(model2, ip, ea, energy, 10)
    # rearrange input files
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
                fname = [file1, file2, file3]
                model1 = GlobalConceptualDFT.from_file(fname, "linear")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_linear(model2, ip, ea, energy, 10)
    # rearrange input files
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file3:
                fname = [file1, file2, file3]
                model1 = GlobalConceptualDFT.from_file(fname, "linear")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "linear")
    # check from_file
    check_global_reactivity_linear(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_linear(model2, ip, ea, energy, 10)


def check_global_reactivity_quadratic(model, ip, ea, energy, n):
    """Check expected quadratic global reactivity descriptors."""
    # check print statement
    assert isinstance(model.__repr__(), str)
    # check energy
    assert_almost_equal(model.energy(n), energy, decimal=6)
    assert_almost_equal(model.energy(n - 1), energy + ip, decimal=6)
    assert_almost_equal(model.energy(n + 1), energy - ea, decimal=6)
    # check ionization-potential & electron-affinity
    assert_almost_equal(model.ip, ip, decimal=6)
    assert_almost_equal(model.ionization_potential, ip, decimal=6)
    assert_almost_equal(model.ea, ea, decimal=6)
    assert_almost_equal(model.electron_affinity, ea, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    mu, eta = -0.5 * (ip + ea), ip - ea
    assert_almost_equal(model.mu, mu, decimal=6)
    assert_almost_equal(model.chemical_potential, mu, decimal=6)
    assert_almost_equal(model.eta, eta, decimal=6)
    assert_almost_equal(model.chemical_hardness, eta, decimal=6)
    assert_almost_equal(model.hyper_hardness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_hardness(4), 0.0, decimal=6)
    # check softness & hyper-softness
    assert_almost_equal(model.softness, 1.0 / eta, decimal=5)
    assert_almost_equal(model.hyper_softness(2), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(3), 0.0, decimal=6)
    assert_almost_equal(model.hyper_softness(4), 0.0, decimal=6)
    # check N0 & N_max
    assert_almost_equal(model.n0, n, decimal=6)
    assert_almost_equal(model.n_max, n - mu / eta, decimal=6)
    # check electrophilicity
    sgn = np.sign(n - mu / eta - n)
    assert_almost_equal(model.electrophilicity, sgn * 0.5 * mu**2 / eta, decimal=6)
    assert_almost_equal(model.electrophilicity, sgn * (ip + ea)**2 / (8 * (ip - ea)), decimal=6)
    # check nucleofugality
    sgn = np.sign(n + 1 - n + mu / eta)
    assert_almost_equal(model.nucleofugality, sgn * (ip - 3 * ea)**2 / (8 * (ip - ea)), decimal=6)
    assert_almost_equal(model.nucleofugality, sgn * (mu + eta)**2 / (2 * eta), decimal=6)
    assert_almost_equal(model.nucleofugality, sgn * (-ea + 0.5 * mu**2 / eta), decimal=6)
    # check electrofugality
    sgn = np.sign(n - mu / eta - n + 1)
    assert_almost_equal(model.electrofugality, sgn * (3 * ip - ea)**2 / (8 * (ip - ea)), decimal=6)
    assert_almost_equal(model.electrofugality, sgn * (mu - eta)**2 / (2 * eta), decimal=6)
    assert_almost_equal(model.electrofugality, sgn * (ip + 0.5 * mu**2 / eta), decimal=6)


def test_global_quadratic_from_file_fmo_ch4_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
        model2 = GlobalConceptualDFT.from_file([fname], "quadratic")
    # check from_file
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_file given as a list
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)


def test_global_quadratic_from_molecule_fmo_ch4_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)


def test_global_quadratic_from_file_fmo_ch4_wfn():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
        model2 = GlobalConceptualDFT.from_file([fname], "quadratic")
    # check from_file
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_file given as a list
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)


def test_global_quadratic_from_molecule_fmo_ch4_wfn():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as fname:
        molecule = Molecule.from_file(fname)
    # check from_molecule
    model = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)
    # check from_molecule given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)


def test_global_quadratic_from_file_fmo_h2o_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-3.09871604E-01), -2.48704636E-02, -7.645980351270224E+01
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
        model2 = GlobalConceptualDFT.from_file([fname], "quadratic")
    # check from_file
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_file given as a list
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)


def test_global_quadratic_from_molecule_fmo_h2o_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-3.09871604E-01), -2.48704636E-02, -7.645980351270224E+01
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        molecule = Molecule.from_file(fname)
    # check from_file
    model = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)
    # check from_file given as a list
    model = GlobalConceptualDFT.from_molecule([molecule], "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 10)


def test_global_quadratic_from_file_fmo_h2o_cation_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-8.47044131E-01), -(-6.19391831E-01), -7.599493522312368E+01
    with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
        model2 = GlobalConceptualDFT.from_file([fname], "quadratic")
    # check quadratic global conceptual DFT model from a fname given as string
    check_global_reactivity_quadratic(model1, ip, ea, energy, 9)
    # check quadratic global conceptual DFT model from a fname given as a list of string
    check_global_reactivity_quadratic(model2, ip, ea, energy, 9)


def test_global_quadratic_from_molecule_fmo_h2o_cation_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -(-8.47044131E-01), -(-6.19391831E-01), -7.599493522312368E+01
    with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    # check quadratic global conceptual DFT model from a fname given as string
    model = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 9)
    # check quadratic global conceptual DFT model from a fname given as a list of string
    model = GlobalConceptualDFT.from_molecule([molecule], "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 9)


def test_global_quadratic_from_file_fmo_h2o_anion_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -1.93118022E-01, -2.69116912E-01, -7.635212549312298E+01
    with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as fname:
        model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
        model2 = GlobalConceptualDFT.from_file([fname], "quadratic")
    # check quadratic global conceptual DFT model from a fname given as string
    check_global_reactivity_quadratic(model1, ip, ea, energy, 11)
    # check quadratic global conceptual DFT model from a fname given as a list of string
    check_global_reactivity_quadratic(model2, ip, ea, energy, 11)


def test_global_quadratic_from_molecule_fmo_h2o_anion_fchk():
    # FMO: ip = -E(HOMO) & ea = -E(LUMO)
    ip, ea, energy = -1.93118022E-01, -2.69116912E-01, -7.635212549312298E+01
    with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    # check quadratic global conceptual DFT model from a fname given as string
    model = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 11)
    # check quadratic global conceptual DFT model from a fname given as a list of string
    model = GlobalConceptualDFT.from_molecule([molecule], "quadratic")
    check_global_reactivity_quadratic(model, ip, ea, energy, 11)


def test_global_quadratic_fd_h2o_fchk():
    ep, e0, en = -7.599493522312368E+01, -7.645980351270224E+01, -7.635212549312298E+01
    ip, ea, energy = ep - e0, e0 - en, e0
    # check quadratic global conceptual DFT model
    with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as file1:
        with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as file2:
            with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as file3:
                fname = [file1, file2, file3]
                model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    # cehck from_file
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)
    # rearrange input files
    with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as file1:
        with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as file2:
            with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as file3:
                fname = [file1, file2, file3]
                # check from_file
                model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)
    # rearrange input files
    with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as file1:
        with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as file2:
            with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as file3:
                fname = [file1, file2, file3]
                model1 = GlobalConceptualDFT.from_file(fname, "quadratic")
                molecule = [Molecule.from_file(item) for item in fname]
                model2 = GlobalConceptualDFT.from_molecule(molecule, "quadratic")
    # check from_file
    check_global_reactivity_quadratic(model1, ip, ea, energy, 10)
    # check from_molecule
    check_global_reactivity_quadratic(model2, ip, ea, energy, 10)


# def test_global_rational_ch4_fchk():
#     file_path = context.get_fn("test/ch4_uhf_ccpvdz.fchk")
#     # ip = -E(HOMO) & ea = -E(LUMO)
#     ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
#     # build global conceptual DFT desp
#     desp = GlobalConceptualDFT.from_file(file_path, desp="rational")
#     assert_almost_equal(desp.energy(10.), energy, decimal=6)
#     assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
#     assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
#     # check ionization potential and electron affinity
#     assert_almost_equal(desp.ip, ip, decimal=6)
#     assert_almost_equal(desp.ionization_potential, ip, decimal=6)
#     assert_almost_equal(desp.ea, ea, decimal=6)
#     assert_almost_equal(desp.electron_affinity, ea, decimal=6)
#     # check parameters; expected parameters were solved for with np.linalg.solve
#     a0, a1, b1 = -39.8993538175, 4.19041197078, -0.104987142594
#     assert_almost_equal(desp.params[0], a0, decimal=6)
#     assert_almost_equal(desp.params[1], a1, decimal=6)
#     assert_almost_equal(desp.params[2], b1, decimal=6)
#     # Check chemical-potential & chemical-hardness; expected derivative values were computed
#     # symbolically
#     mu, eta = 0.600211746260527, -2.52707898353937
#     assert_almost_equal(desp.mu, mu, decimal=6)
#     assert_almost_equal(desp.chemical_potential, mu, decimal=6)
#     assert_almost_equal(desp.eta, eta, decimal=6)
#     assert_almost_equal(desp.chemical_hardness, eta, decimal=6)
#     # Check derivative of E(N); expected derivative values were computed symbolically
#     assert_almost_equal(desp.energy_derivative(10, 1), mu, decimal=6)
#     assert_almost_equal(desp.energy_derivative(10, 2), eta, decimal=6)
#     assert_almost_equal(desp.energy_derivative(9.60, 1), 24.0621211372066, decimal=6)
#     assert_almost_equal(desp.energy_derivative(9.10, 2), 3.52917348462076, decimal=6)
#     assert_almost_equal(desp.energy_derivative(8.50, 1), 0.128916512086747, decimal=6)
#     assert_almost_equal(desp.energy_derivative(10.3, 1), 0.225478627953592, decimal=6)
#     assert_almost_equal(desp.energy_derivative(10.8, 2), -0.130680448326993, decimal=6)
#     assert_almost_equal(desp.energy_derivative(11.5, 1), 0.0347209036078056, decimal=6)
#     assert_almost_equal(desp.energy_derivative(10.2, 3), 3.91390151199886, decimal=6)
#     # Check hyper-hardness; expected derivative values were computed symbolically
#     assert_almost_equal(desp.hyper_hardness(2), 15.9596881321556, decimal=6)
#     assert_almost_equal(desp.hyper_hardness(3), -134.390547049114, decimal=6)
#     assert_almost_equal(desp.hyper_hardness(4), 1414.56548105706, decimal=6)
#     assert_almost_equal(desp.hyper_hardness(5), -17867.2879377639, decimal=6)
#     # Check softness & hyper-softness; expected derivative values were computed symbolically
#     assert_almost_equal(desp.softness, 1.0 / eta, decimal=6)
#     # assert_almost_equal(desp.hyper_softness(2), value, decimal=6)
#     # assert_almost_equal(desp.hyper_softness(3), value, decimal=6)
#     # Check N_max and related descriptors
