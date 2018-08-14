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
"""Test chemtools.analysis.conceptual."""

import numpy as np

from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import IOData, BeckeMolGrid
from chemtools import context
from chemtools.toolbox.conceptual import GlobalConceptualDFT, LocalConceptualDFT
from chemtools.toolbox.conceptual import CondensedConceptualDFT


def test_global_conceptual_raises():
    # dictionary of energy values
    values = {0.0: 0.0, 1.0: -0.5, 2.0: -0.45}
    # check invalid model
    assert_raises(ValueError, GlobalConceptualDFT, values, 'quad')
    assert_raises(ValueError, GlobalConceptualDFT, values, 'gibberish')
    assert_raises(NotImplementedError, GlobalConceptualDFT, values, 'general')
    # check invalid coordinates
    assert_raises(ValueError, GlobalConceptualDFT, values, 'linear', np.array([0., 0., 0.]))
    assert_raises(ValueError, GlobalConceptualDFT, values, 'linear', np.array([[0., 0.]]))
    # check invalid atomic numbers
    coord = np.array([[0., 0., 0.]])
    assert_raises(ValueError, GlobalConceptualDFT, values, 'linear', coord, [])
    assert_raises(ValueError, GlobalConceptualDFT, values, 'linear', coord, [1, 1])
    # check invalid number of electrons
    values = {'0': 0.0, 1: -0.5, 2: -0.45}
    assert_raises(ValueError, GlobalConceptualDFT, values, 'quadratic', coord)
    values = {0: 0.0, 1: -0.5, 3: -0.45}
    assert_raises(ValueError, GlobalConceptualDFT, values, 'quadratic')
    values = {0: 0.0, 1: -0.5, 2: -0.45, 3: -0.4}
    assert_raises(ValueError, GlobalConceptualDFT, values, 'quadratic', coord)
    # check non-existing attribute
    model = GlobalConceptualDFT({0.0: 0.0, 1.0: -0.5, 2.0: -0.45}, 'quadratic')
    assert_raises(AttributeError, getattr, model, 'mu_plus')
    assert_raises(AttributeError, getattr, model, 'gibberish')
    # check molecule file inconsistency
    fnames = [context.get_fn('test/ch4_uhf_ccpvdz.fchk'), context.get_fn('test/o2_uhf.fchk')]
    assert_raises(ValueError, GlobalConceptualDFT.from_file, fnames, 'linear')
    fname = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    assert_raises(ValueError, GlobalConceptualDFT.from_file, [fname, fname], 'linear')


def check_global_linear_fmo_ch4_uhf_ccpvdz(filename):
    """Check expected linear global indicators for ch4_uhf_ccpvdz within FMO approach."""
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(context.get_fn(filename), model='linear')
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # check energy values
    assert_almost_equal(desp.energy(10.), energy, decimal=6)
    assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
    assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
    # check ionization potential and electron affinity
    assert_almost_equal(desp.ip, ip, decimal=6)
    assert_almost_equal(desp.ionization_potential, ip, decimal=6)
    assert_almost_equal(desp.ea, ea, decimal=6)
    assert_almost_equal(desp.electron_affinity, ea, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    assert_equal(desp.mu, None)
    assert_equal(desp.chemical_potential, None)
    assert_equal(desp.eta, None)
    assert_equal(desp.chemical_hardness, None)
    assert_equal(desp.hyper_hardness(2), None)
    assert_equal(desp.hyper_hardness(3), None)
    assert_equal(desp.hyper_hardness(4), None)
    # check mu+, mu-, mu0
    assert_almost_equal(desp.mu_plus, -ea, decimal=6)
    assert_almost_equal(desp.mu_minus, -ip, decimal=6)
    assert_almost_equal(desp.mu_zero, -0.5 * (ip + ea), decimal=6)
    # check derivatives of energy w.r.t. number of electrons
    assert_equal(desp.energy_derivative(10, 3), None)
    assert_almost_equal(desp.energy_derivative(9.5, 4), 0.0, decimal=6)
    assert_almost_equal(desp.energy_derivative(9.0, 3), 0.0, decimal=6)
    assert_almost_equal(desp.energy_derivative(10.4, 2), 0.0, decimal=6)
    assert_almost_equal(desp.energy_derivative(11, 1), -ea, decimal=6)
    assert_almost_equal(desp.energy_derivative(9.0, 1), -ip, decimal=6)
    assert_almost_equal(desp.energy_derivative(9.7, 1), -ip, decimal=6)
    assert_almost_equal(desp.energy_derivative(10.5, 1), -ea, decimal=6)


def test_global_linear_fmo_ch4_uhf_ccpvdz_fchk():
    check_global_linear_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.fchk')


def test_global_linear_fmo_ch4_uhf_ccpvdz_wfn():
    check_global_linear_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.wfn')


def check_global_quadratic_fmo_ch4_uhf_ccpvdz(filename):
    """Check expected quadratic global indicators for ch4_uhf_ccpvdz within FMO approach."""
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(context.get_fn(filename), model='quadratic')
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # check energy
    assert_almost_equal(desp.energy(10.), energy, decimal=6)
    assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
    assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
    # check ionization-potential & electron-affinity
    assert_almost_equal(desp.ip, ip, decimal=6)
    assert_almost_equal(desp.ionization_potential, ip, decimal=6)
    assert_almost_equal(desp.ea, ea, decimal=6)
    assert_almost_equal(desp.electron_affinity, ea, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    mu, eta = -0.5 * (ip + ea), ip - ea
    assert_almost_equal(desp.mu, mu, decimal=6)
    assert_almost_equal(desp.chemical_potential, mu, decimal=6)
    assert_almost_equal(desp.eta, eta, decimal=6)
    assert_almost_equal(desp.chemical_hardness, eta, decimal=6)
    assert_almost_equal(desp.hyper_hardness(2), 0.0, decimal=6)
    assert_almost_equal(desp.hyper_hardness(3), 0.0, decimal=6)
    assert_almost_equal(desp.hyper_hardness(4), 0.0, decimal=6)
    # check softness & hyper-softness
    assert_almost_equal(desp.softness, 1.0 / eta, decimal=5)
    # assert_almost_equal(desp.hyper_softness(2), 0.0, decimal=6)
    # assert_almost_equal(desp.hyper_softness(3), 0.0, decimal=6)
    # assert_almost_equal(desp.hyper_softness(4), 0.0, decimal=6)
    # check N_max and related descriptors
    assert_almost_equal(desp.n0, 10, decimal=6)
    assert_almost_equal(desp.n_max, 10 - mu / eta, decimal=6)
    # check Electrophilicity
    value = 0.5 * mu * mu / eta
    assert_almost_equal(desp.electrophilicity, value, decimal=6)
    value = (ip + ea)**2 / (8 * (ip - ea))
    assert_almost_equal(desp.electrophilicity, value, decimal=6)
    # check Nucleofugality
    value = (ip - 3 * ea)**2 / (8 * (ip - ea))
    assert_almost_equal(desp.nucleofugality, value, decimal=6)
    value = (mu + eta)**2 / (2 * eta)
    assert_almost_equal(desp.nucleofugality, value, decimal=6)
    value = - ea + 0.5 * mu * mu / eta
    assert_almost_equal(desp.nucleofugality, value, decimal=6)
    # check Electrofugality
    value = (3 * ip - ea)**2 / (8 * (ip - ea))
    assert_almost_equal(desp.electrofugality, value, decimal=6)
    value = (mu - eta)**2 / (2 * eta)
    assert_almost_equal(desp.electrofugality, value, decimal=6)
    value = ip + 0.5 * mu * mu / eta
    assert_almost_equal(desp.electrofugality, value, decimal=6)


def test_global_quadratic_fmo_ch4_uhf_ccpvdz_fchk():
    check_global_quadratic_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.fchk')


def test_global_quadratic_fmo_ch4_uhf_ccpvdz_wfn():
    check_global_quadratic_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.wfn')


# def test_global_rational_ch4_fchk():
#     file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
#     # ip = -E(homo) & ea = E(lumo)
#     ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
#     # build global conceptual DFT tool
#     desp = GlobalConceptualDFT.from_file(file_path, model='rational')
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
