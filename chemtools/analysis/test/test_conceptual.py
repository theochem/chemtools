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
#pylint: skip-file


import os
from chemtools import *


def test_global_linear_ch4_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(file_path, model='linear')
    # check energy values
    np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=8)
    # check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=8)
    # check chemical-potential, chemical-hardness & hyper-hardness
    np.testing.assert_equal(desp.mu, None)
    np.testing.assert_equal(desp.chemical_potential, None)
    np.testing.assert_equal(desp.eta, None)
    np.testing.assert_equal(desp.chemical_hardness, None)
    np.testing.assert_equal(desp.hyper_hardness(2), None)
    np.testing.assert_equal(desp.hyper_hardness(3), None)
    np.testing.assert_equal(desp.hyper_hardness(4), None)
    # check mu+, mu-, mu0
    np.testing.assert_almost_equal(desp.mu_plus, -ea, decimal=8)
    np.testing.assert_almost_equal(desp.mu_minus, -ip, decimal=8)
    np.testing.assert_almost_equal(desp.mu_zero, -0.5 * (ip + ea), decimal=8)
    # check derivatives of energy w.r.t. number of electrons
    np.testing.assert_equal(desp.energy_derivative(10, 3), None)
    np.testing.assert_almost_equal(desp.energy_derivative(9.5, 4), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(9.0, 3), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10.4, 2), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(11, 1), -ea, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(9.0, 1), -ip, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(9.7, 1), -ip, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10.5, 1), -ea, decimal=8)


def test_local_linear_ch4_fchk():
    # Check softness & hyper-softness
    # Check N_max and related descriptors
    # temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec='exp:5e-4:2e1:175:434',
                        random_rotate=False, mode='keep')
    # build local conceptual DFT tool
    desp = LocalConceptualDFT.from_file(file_path, model='linear', points=grid.points)
    # check shape of density
    np.testing.assert_equal(desp.density_zero.shape, grid.shape)
    np.testing.assert_equal(desp.density_plus.shape, grid.shape)
    np.testing.assert_equal(desp.density_minus.shape, grid.shape)
    # check integral of density
    np.testing.assert_almost_equal(grid.integrate(desp.density_zero), 10., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density_plus), 11., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density_minus), 9., decimal=4)
    # check shape of Fukui functions & dual descriptor
    np.testing.assert_equal(desp.ff_zero.shape, grid.shape)
    np.testing.assert_equal(desp.ff_plus.shape, grid.shape)
    np.testing.assert_equal(desp.ff_minus.shape, grid.shape)
    # check Fukui functions & dual descriptor
    np.testing.assert_almost_equal(grid.integrate(desp.ff_plus), 1., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.ff_minus), 1., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.ff_zero), 1., decimal=4)
    # check dual descriptor


def test_gloabl_quadratic_ch4_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(file_path, model='quadratic')
    # check energy
    np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=8)
    # check ionization-potential & electron-affinity
    np.testing.assert_almost_equal(desp.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=8)
    # check chemical-potential, chemical-hardness & hyper-hardness
    mu, eta = -0.5 * (ip + ea), ip - ea
    np.testing.assert_almost_equal(desp.mu, mu, decimal=8)
    np.testing.assert_almost_equal(desp.chemical_potential, mu, decimal=8)
    np.testing.assert_almost_equal(desp.eta, eta, decimal=8)
    np.testing.assert_almost_equal(desp.chemical_hardness, eta, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(2), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(3), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(4), 0.0, decimal=8)
    # check softness & hyper-softness
    np.testing.assert_almost_equal(desp.softness, 1.0/eta, decimal=8)
    # np.testing.assert_almost_equal(desp.hyper_softness(2), 0.0, decimal=8)
    # np.testing.assert_almost_equal(desp.hyper_softness(3), 0.0, decimal=8)
    # np.testing.assert_almost_equal(desp.hyper_softness(4), 0.0, decimal=8)
    # check N_max and related descriptors
    np.testing.assert_almost_equal(desp.n0, 10, decimal=8)
    np.testing.assert_almost_equal(desp.n_max, 10 - mu/eta, decimal=8)
    # check Electrophilicity
    value = 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.electrophilicity, value, decimal=8)
    value = (ip + ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.electrophilicity, value, decimal=8)
    # check Nucleofugality
    value = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=8)
    value = (mu + eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=8)
    value = - ea + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=8)
    # check Electrofugality
    value = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=8)
    value = (mu - eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=8)
    value = ip + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=8)


def test_local_quadratic_ch4_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    mu, eta = -0.5 * (ip + ea), ip - ea
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec='exp:5e-4:2e1:175:434',
                        random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = LocalConceptualDFT.from_file(file_path, model='quadratic', points=grid.points)
    # check shape of density
    np.testing.assert_equal(desp.density_zero.shape, grid.shape)
    np.testing.assert_equal(desp.density_plus.shape, grid.shape)
    np.testing.assert_equal(desp.density_minus.shape, grid.shape)
    # check integral of density
    np.testing.assert_almost_equal(grid.integrate(desp.density_zero), 10., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density_plus), 11., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density_minus), 9., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density(10.62)), 10.62, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density(9.78)), 9.78, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.density(10.0)), 10.0, decimal=4)
    # Check Fukui function, dual descriptor & softness
    np.testing.assert_almost_equal(grid.integrate(desp.fukui_function()), 1., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.fukui_function(10.5)), 1., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.fukui_function(9.50)), 1., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.dual_descriptor()), 0., decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.dual_descriptor(10.79)), 0., decimal=4)
    # Check local softness
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta, 10.3)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta, 9.1)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.hyper_softness(eta)), 0., decimal=3)
    np.testing.assert_almost_equal(grid.integrate(desp.hyper_softness(eta, 9.91)), 0., decimal=3)


# def test_condense_quadratic_ch4_fchk():
#     # Temporary trick to find the data files
#     path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
#     file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
#     # make molecular grid
#     mol = IOData.from_file(file_path)
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec='exp:5e-4:2e1:175:434',
#                         random_rotate=False, mode='keep')
#     # build global conceptual DFT tool
#     desp = LocalConceptualDFT.from_file([file_path], model='quadratic', points=grid.points)
#     # Check condensed dual descriptors (Becke part only)
#     # TODO: How were the expected values calculated?
#     c, h1, h2, h3, h4 = -0.26854311,  0.05276027,  0.09886118, -0.03029482,  0.14726817
#     condens = desp.condensedtool.dual_descriptor()
#     np.testing.assert_almost_equal(condens[0], c, decimal=4)
#     np.testing.assert_almost_equal(condens[1], h1, decimal=4)
#     np.testing.assert_almost_equal(condens[2], h2, decimal=4)
#     np.testing.assert_almost_equal(condens[3], h3, decimal=4)
#     np.testing.assert_almost_equal(condens[4], h4, decimal=4)


def test_global_rational_ch4_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(file_path, model='rational')
    np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=8)
    # check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=8)
    # check parameters; expected parameters were solved for with np.linalg.solve
    a0, a1, b1 = -39.8993538175, 4.19041197078, -0.104987142594
    np.testing.assert_almost_equal(desp._a0, a0, decimal=8)
    np.testing.assert_almost_equal(desp._a1, a1, decimal=8)
    np.testing.assert_almost_equal(desp._b1, b1, decimal=8)
    # Check chemical-potential & chemical-hardness; expected derivative values were computed symbolically
    mu, eta = 0.600211746260527, -2.52707898353937
    np.testing.assert_almost_equal(desp.mu, mu, decimal=8)
    np.testing.assert_almost_equal(desp.chemical_potential, mu, decimal=8)
    np.testing.assert_almost_equal(desp.eta, eta, decimal=8)
    np.testing.assert_almost_equal(desp.chemical_hardness, eta, decimal=8)
    # Check derivative of E(N); expected derivative values were computed symbolically
    np.testing.assert_almost_equal(desp.energy_derivative(10, 1), mu, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10, 2), eta, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(9.60, 1), 24.0621211372066, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(9.10, 2), 3.52917348462076, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(8.50, 1), 0.128916512086747, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10.3, 1), 0.225478627953592, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10.8, 2), -0.130680448326993, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(11.5, 1), 0.0347209036078056, decimal=8)
    np.testing.assert_almost_equal(desp.energy_derivative(10.2, 3), 3.91390151199886, decimal=8)
    # Check hyper-hardness; expected derivative values were computed symbolically
    np.testing.assert_almost_equal(desp.hyper_hardness(2), 15.9596881321556, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(3), -134.390547049114, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(4), 1414.56548105706, decimal=8)
    np.testing.assert_almost_equal(desp.hyper_hardness(5), -17867.2879377639, decimal=8)
    # Check softness & hyper-softness; expected derivative values were computed symbolically
    np.testing.assert_almost_equal(desp.softness, 1.0 / eta, decimal=8)
    # np.testing.assert_almost_equal(desp.hyper_softness(2), value, decimal=8)
    # np.testing.assert_almost_equal(desp.hyper_softness(3), value, decimal=8)
    # Check N_max and related descriptors
