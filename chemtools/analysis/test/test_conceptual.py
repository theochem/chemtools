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
from horton import IOData, BeckeMolGrid
from chemtools import context
from chemtools.analysis.conceptual import (GlobalConceptualDFT, LocalConceptualDFT,
                                           CondensedConceptualDFT)


def test_global_linear_ch4_fchk():
    # use context to get path
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(file_path, model='linear')
    # check energy values
    np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=6)
    np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
    np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
    # check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.ip, ip, decimal=6)
    np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(desp.ea, ea, decimal=6)
    np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    np.testing.assert_equal(desp.mu, None)
    np.testing.assert_equal(desp.chemical_potential, None)
    np.testing.assert_equal(desp.eta, None)
    np.testing.assert_equal(desp.chemical_hardness, None)
    np.testing.assert_equal(desp.hyper_hardness(2), None)
    np.testing.assert_equal(desp.hyper_hardness(3), None)
    np.testing.assert_equal(desp.hyper_hardness(4), None)
    # check mu+, mu-, mu0
    np.testing.assert_almost_equal(desp.mu_plus, -ea, decimal=6)
    np.testing.assert_almost_equal(desp.mu_minus, -ip, decimal=6)
    np.testing.assert_almost_equal(desp.mu_zero, -0.5 * (ip + ea), decimal=6)
    # check derivatives of energy w.r.t. number of electrons
    np.testing.assert_equal(desp.energy_derivative(10, 3), None)
    np.testing.assert_almost_equal(desp.energy_derivative(9.5, 4), 0.0, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(9.0, 3), 0.0, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(10.4, 2), 0.0, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(11, 1), -ea, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(9.0, 1), -ip, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(9.7, 1), -ip, decimal=6)
    np.testing.assert_almost_equal(desp.energy_derivative(10.5, 1), -ea, decimal=6)


def test_local_linear_ch4_fchk():
    # Check softness & hyper-softness
    # Check N_max and related descriptors
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='exp:5e-4:2e1:175:434', random_rotate=False, mode='keep')
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


def test_global_quadratic_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
    # build global conceptual DFT tool
    desp = GlobalConceptualDFT.from_file(file_path, model='quadratic')
    # check energy
    np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=6)
    np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
    np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
    # check ionization-potential & electron-affinity
    np.testing.assert_almost_equal(desp.ip, ip, decimal=6)
    np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=6)
    np.testing.assert_almost_equal(desp.ea, ea, decimal=6)
    np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=6)
    # check chemical-potential, chemical-hardness & hyper-hardness
    mu, eta = -0.5 * (ip + ea), ip - ea
    np.testing.assert_almost_equal(desp.mu, mu, decimal=6)
    np.testing.assert_almost_equal(desp.chemical_potential, mu, decimal=6)
    np.testing.assert_almost_equal(desp.eta, eta, decimal=6)
    np.testing.assert_almost_equal(desp.chemical_hardness, eta, decimal=6)
    np.testing.assert_almost_equal(desp.hyper_hardness(2), 0.0, decimal=6)
    np.testing.assert_almost_equal(desp.hyper_hardness(3), 0.0, decimal=6)
    np.testing.assert_almost_equal(desp.hyper_hardness(4), 0.0, decimal=6)
    # check softness & hyper-softness
    np.testing.assert_almost_equal(desp.softness, 1.0 / eta, decimal=6)
    # np.testing.assert_almost_equal(desp.hyper_softness(2), 0.0, decimal=6)
    # np.testing.assert_almost_equal(desp.hyper_softness(3), 0.0, decimal=6)
    # np.testing.assert_almost_equal(desp.hyper_softness(4), 0.0, decimal=6)
    # check N_max and related descriptors
    np.testing.assert_almost_equal(desp.n0, 10, decimal=6)
    np.testing.assert_almost_equal(desp.n_max, 10 - mu / eta, decimal=6)
    # check Electrophilicity
    value = 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.electrophilicity, value, decimal=6)
    value = (ip + ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.electrophilicity, value, decimal=6)
    # check Nucleofugality
    value = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=6)
    value = (mu + eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=6)
    value = - ea + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.nucleofugality, value, decimal=6)
    # check Electrofugality
    value = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=6)
    value = (mu - eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=6)
    value = ip + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.electrofugality, value, decimal=6)


def test_local_quadratic_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # ip = -E(homo) & ea = E(lumo)
    ip, ea = -(-5.43101269E-01), -1.93295185E-01
    eta = ip - ea
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='exp:5e-4:2e1:175:434', random_rotate=False, mode='keep')
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
    # Check local softness
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta, 10.3)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.softness(1./eta, 9.1)), 1./eta, decimal=4)
    np.testing.assert_almost_equal(grid.integrate(desp.hyper_softness(eta)), 0., decimal=3)


def test_condense_mbis_quadratic_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec='insane',
                        random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file([file_path], 'quadratic', grid=grid, scheme='mbis')
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    np.testing.assert_almost_equal(desp.density_zero, expected, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 11., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 10., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 9.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)
    # check condensed dual descriptor
    np.testing.assert_almost_equal(np.sum(desp.dual_descriptor()), 0.0, decimal=2)


def test_condense_mbis_linear_fmr_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'linear', 'FMR', 'mbis', grid=grid)
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    np.testing.assert_almost_equal(desp.density_zero, expected, decimal=4)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 11., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 10., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 9.0, decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_mbis_linear_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'linear', 'FMR', 'mbis', grid=grid)
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    np.testing.assert_almost_equal(desp.density_zero, expected, decimal=4)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 11., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 10., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 9.0, decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_h_linear_fd_rmf_ch2o_fchk():
    file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # make molecular grid
    mol = IOData.from_file(file_path[0])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, model='linear', grid=grid,
                                            scheme='h', approach='RMF')
    expectedm = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
    expected0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
    expectedp = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
    # check charges
    np.testing.assert_almost_equal(desp.density_plus, expectedp, decimal=2)
    np.testing.assert_almost_equal(desp.density_zero, expected0, decimal=2)
    np.testing.assert_almost_equal(desp.density_minus, expectedm, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 17., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 16., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 15., decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    np.testing.assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    np.testing.assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_h_linear_fd_fmr_ch20_fchk():
    file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # make molecular grid
    mol = IOData.from_file(file_path[0])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, model='linear', grid=grid,
                                            scheme='h', approach='FMR')
    expectedm = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
    expected0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
    expectedp = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
    # check charges
    np.testing.assert_almost_equal(desp.density_plus, expectedp, decimal=2)
    np.testing.assert_almost_equal(desp.density_zero, expected0, decimal=2)
    np.testing.assert_almost_equal(desp.density_minus, expectedm, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 17., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 16., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 15., decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    np.testing.assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    np.testing.assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_mbis_linear_fd_rmf_ch2o_fchk():
    file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # make molecular grid
    mol = IOData.from_file(file_path[0])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, model='linear', grid=grid, scheme='mbis',
                                            approach='RMF')
    expectedm = np.array([7.8580338, 5.70425809, 0.71885554, 0.71884404])
    expected0 = np.array([8.41149, 5.66445074, 0.96204946, 0.96202722])
    expectedp = np.array([8.13881352, 6.81770852, 1.28123219, 0.76225513])
    # check charges
    np.testing.assert_almost_equal(desp.density_plus, expectedp, decimal=2)
    np.testing.assert_almost_equal(desp.density_zero, expected0, decimal=2)
    np.testing.assert_almost_equal(desp.density_minus, expectedm, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 17., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 16., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 15., decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    np.testing.assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    np.testing.assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_mbis_linear_fd_fmr_ch4_fchk():
    file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # make molecular grid
    mol = IOData.from_file(file_path[0])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, model='linear', grid=grid, scheme='mbis',
                                            approach='FMR')
    expected0 = np.array([8.41149, 5.66445074, 0.96204946, 0.96202722])
    np.testing.assert_almost_equal(desp.density_zero, expected0, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 17., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 16., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 15., decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)


def test_condense_h_quadratic_fd_ch4_fchk():
    file_path = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
                 context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # make molecular grid
    mol = IOData.from_file(file_path[0])
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'quadratic', 'RMF', 'h', grid=grid)
    # computed with horton separately
    expectedm = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
    expected0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
    expectedp = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
    np.testing.assert_almost_equal(desp.density_plus, expectedp, decimal=2)
    np.testing.assert_almost_equal(desp.density_zero, expected0, decimal=2)
    np.testing.assert_almost_equal(desp.density_minus, expectedm, decimal=2)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 17., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 16., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 15., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)
    # check condensed dual descriptor
    np.testing.assert_almost_equal(np.sum(desp.dual_descriptor()), 0.0, decimal=2)


def test_condense_h_linear_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'linear', 'FMR', 'h', grid=grid)
    # computed with horton separately
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    np.testing.assert_almost_equal(desp.density_zero, expected, decimal=4)
    # check condensed density
    np.testing.assert_almost_equal(np.sum(desp.density_plus), 11., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_zero), 10., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.density_minus), 9.0, decimal=2)
    # check condensed Fukui function
    np.testing.assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    np.testing.assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    condense = lambda x: np.sum(desp.density(x))
    np.testing.assert_almost_equal(condense(15.5), 15.5, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 16.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    condense = lambda x: np.sum(desp.fukui_function(x))
    np.testing.assert_almost_equal(condense(15.5), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.0), 1.0, decimal=2)
    np.testing.assert_almost_equal(condense(16.5), 1.0, decimal=2)

# def test_global_rational_ch4_fchk():
#     file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
#     # ip = -E(homo) & ea = E(lumo)
#     ip, ea, energy = -(-5.43101269E-01), -1.93295185E-01, -4.019868797400735E+01
#     # build global conceptual DFT tool
#     desp = GlobalConceptualDFT.from_file(file_path, model='rational')
#     np.testing.assert_almost_equal(desp.energy(10.), energy, decimal=6)
#     np.testing.assert_almost_equal(desp.energy(9.), energy + ip, decimal=6)
#     np.testing.assert_almost_equal(desp.energy(11.), energy - ea, decimal=6)
#     # check ionization potential and electron affinity
#     np.testing.assert_almost_equal(desp.ip, ip, decimal=6)
#     np.testing.assert_almost_equal(desp.ionization_potential, ip, decimal=6)
#     np.testing.assert_almost_equal(desp.ea, ea, decimal=6)
#     np.testing.assert_almost_equal(desp.electron_affinity, ea, decimal=6)
#     # check parameters; expected parameters were solved for with np.linalg.solve
#     a0, a1, b1 = -39.8993538175, 4.19041197078, -0.104987142594
#     np.testing.assert_almost_equal(desp.params[0], a0, decimal=6)
#     np.testing.assert_almost_equal(desp.params[1], a1, decimal=6)
#     np.testing.assert_almost_equal(desp.params[2], b1, decimal=6)
#     # Check chemical-potential & chemical-hardness; expected derivative values were computed
#     # symbolically
#     mu, eta = 0.600211746260527, -2.52707898353937
#     np.testing.assert_almost_equal(desp.mu, mu, decimal=6)
#     np.testing.assert_almost_equal(desp.chemical_potential, mu, decimal=6)
#     np.testing.assert_almost_equal(desp.eta, eta, decimal=6)
#     np.testing.assert_almost_equal(desp.chemical_hardness, eta, decimal=6)
#     # Check derivative of E(N); expected derivative values were computed symbolically
#     np.testing.assert_almost_equal(desp.energy_derivative(10, 1), mu, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(10, 2), eta, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(9.60, 1), 24.0621211372066, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(9.10, 2), 3.52917348462076, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(8.50, 1), 0.128916512086747, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(10.3, 1), 0.225478627953592, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(10.8, 2), -0.130680448326993, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(11.5, 1), 0.0347209036078056, decimal=6)
#     np.testing.assert_almost_equal(desp.energy_derivative(10.2, 3), 3.91390151199886, decimal=6)
#     # Check hyper-hardness; expected derivative values were computed symbolically
#     np.testing.assert_almost_equal(desp.hyper_hardness(2), 15.9596881321556, decimal=6)
#     np.testing.assert_almost_equal(desp.hyper_hardness(3), -134.390547049114, decimal=6)
#     np.testing.assert_almost_equal(desp.hyper_hardness(4), 1414.56548105706, decimal=6)
#     np.testing.assert_almost_equal(desp.hyper_hardness(5), -17867.2879377639, decimal=6)
#     # Check softness & hyper-softness; expected derivative values were computed symbolically
#     np.testing.assert_almost_equal(desp.softness, 1.0 / eta, decimal=6)
#     # np.testing.assert_almost_equal(desp.hyper_softness(2), value, decimal=6)
#     # np.testing.assert_almost_equal(desp.hyper_softness(3), value, decimal=6)
#     # Check N_max and related descriptors
