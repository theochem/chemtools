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
"""Test chemtools.analysis.conceptual.CondensedConceptualDFT."""


import numpy as np

from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import IOData, BeckeMolGrid
from chemtools import context
from chemtools.toolbox.conceptual import CondensedConceptualDFT


def test_condensed_conceptual_raises():
    # file for FMO
    fname = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # files for FD
    fnames = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
              context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
              context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]
    # check invalid scheme
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, 'quadratic', scheme='gib')
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, 'linear', scheme='err')
    # check invalid grid type
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, 'quadratic', grid='str')
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, 'linear', grid='gibb')
    # check invalid grid coordinates
    grid = BeckeMolGrid(np.array([[0., 0., 0.]]), np.array([1]))
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, 'quadratic', grid=grid)
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, 'quadratic', grid=grid)
    # check invalid grid atomic numbers
    coord = np.array([[-3.77945227E-05,  3.77945227E-05, -1.88972613E-05],
                      [ 1.04290206E+00,  1.50497789E+00,  9.34507367E-01],
                      [ 1.28607202E+00, -1.53098052E+00, -4.77307027E-01],
                      [-1.46467003E+00, -7.02997019E-01,  1.25954026E+00],
                      [-8.64474117E-01,  7.29131931E-01, -1.71670281E+00]])
    grid = BeckeMolGrid(coord, np.array([7, 1, 1, 1, 1]))
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, 'quadratic', grid=grid)
    assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, 'quadratic', grid=grid)


def test_condense_mbis_quadratic_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec='insane',
                        random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file([file_path], 'quadratic', grid=grid, scheme='mbis')
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    assert_almost_equal(desp.density(10.), expected, decimal=4)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(11.)), 11., decimal=2)
    assert_almost_equal(np.sum(desp.density(10.)), 10., decimal=2)
    assert_almost_equal(np.sum(desp.density(9.0)), 9.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)
    # check condensed dual descriptor
    assert_almost_equal(np.sum(desp.dual_descriptor), 0.0, decimal=2)


def test_condense_mbis_linear_fmr_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'linear', 'FMR', 'mbis', grid=grid)
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    assert_almost_equal(desp.density(10.), expected, decimal=4)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(11.)), 11., decimal=2)
    assert_almost_equal(np.sum(desp.density(10.)), 10., decimal=2)
    assert_almost_equal(np.sum(desp.density(9.0)), 9.0, decimal=2)
    # check condensed Fukui function
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


def test_condense_mbis_linear_ch4_fchk():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # make molecular grid
    mol = IOData.from_file(file_path)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='insane', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = CondensedConceptualDFT.from_file(file_path, 'linear', 'FMR', 'mbis', grid=grid)
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # computed with horton separately
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    assert_almost_equal(desp.density(10.), expected, decimal=4)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(11.)), 11., decimal=2)
    assert_almost_equal(np.sum(desp.density(10.)), 10., decimal=2)
    assert_almost_equal(np.sum(desp.density(9.0)), 9.0, decimal=2)
    # check condensed Fukui function
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1.0, decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


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
    assert_almost_equal(desp.density(17.), expectedp, decimal=2)
    assert_almost_equal(desp.density(16.), expected0, decimal=2)
    assert_almost_equal(desp.density(15.), expectedm, decimal=2)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(17.)), 17., decimal=2)
    assert_almost_equal(np.sum(desp.density(16.)), 16., decimal=2)
    assert_almost_equal(np.sum(desp.density(15.)), 15., decimal=2)
    # check condensed Fukui function
    assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


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
    assert_almost_equal(desp.density(17.), expectedp, decimal=2)
    assert_almost_equal(desp.density(16.), expected0, decimal=2)
    assert_almost_equal(desp.density(15.), expectedm, decimal=2)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(17.)), 17., decimal=2)
    assert_almost_equal(np.sum(desp.density(16.)), 16., decimal=2)
    assert_almost_equal(np.sum(desp.density(15.)), 15., decimal=2)
    # check condensed Fukui function
    assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


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
    assert_almost_equal(desp.density(17.), expectedp, decimal=2)
    assert_almost_equal(desp.density(16.), expected0, decimal=2)
    assert_almost_equal(desp.density(15.), expectedm, decimal=2)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(17.)), 17., decimal=2)
    assert_almost_equal(np.sum(desp.density(16.)), 16., decimal=2)
    assert_almost_equal(np.sum(desp.density(15.)), 15., decimal=2)
    # check condensed Fukui function
    assert_almost_equal(desp.ff_plus, expectedp - expected0, decimal=2)
    assert_almost_equal(desp.ff_zero, 0.5 * (expectedp - expectedm), decimal=2)
    assert_almost_equal(desp.ff_minus, expected0 - expectedm, decimal=2)
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


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
    assert_almost_equal(desp.density(16.), expected0, decimal=2)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(17.)), 17., decimal=2)
    assert_almost_equal(np.sum(desp.density(16.)), 16., decimal=2)
    assert_almost_equal(np.sum(desp.density(15.)), 15., decimal=2)
    # check condensed Fukui function
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)


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
    assert_almost_equal(desp.density(17.), expectedp, decimal=2)
    assert_almost_equal(desp.density(16.), expected0, decimal=2)
    assert_almost_equal(desp.density(15.), expectedm, decimal=2)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(17.)), 17., decimal=2)
    assert_almost_equal(np.sum(desp.density(16.)), 16., decimal=2)
    assert_almost_equal(np.sum(desp.density(15.)), 15., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)
    # check condensed dual descriptor
    assert_almost_equal(np.sum(desp.dual_descriptor), 0.0, decimal=2)


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
    assert_almost_equal(desp.density(10.), expected, decimal=4)
    # check condensed density
    assert_almost_equal(np.sum(desp.density(11.)), 11., decimal=2)
    assert_almost_equal(np.sum(desp.density(10.)), 10., decimal=2)
    assert_almost_equal(np.sum(desp.density(9.0)), 9.0, decimal=2)
    # check condensed Fukui function
    assert_almost_equal(np.sum(desp.ff_plus), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_zero), 1., decimal=2)
    assert_almost_equal(np.sum(desp.ff_minus), 1., decimal=2)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.density(15.5)), 15.5, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.0)), 16.0, decimal=2)
    assert_almost_equal(np.sum(desp.density(16.5)), 16.5, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(desp.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(15.5, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(desp.density_derivative(16.5, 1)), 1.0, decimal=2)

