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
"""Test chemtools.analysis.conceptual.LocalConceptualDFT"""


import numpy as np

from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import IOData, BeckeMolGrid
from chemtools import context
from chemtools.toolbox.conceptual import LocalConceptualDFT


def test_local_conceptual_raises():
    # check invalid densities
    values = {1.0: np.array([0.0, 0.5]), 2.0: np.array([1.0]), 3.0: np.array([[2.0]])}
    assert_raises(ValueError, LocalConceptualDFT, values)
    # check in valid model
    values = {1.0: np.array([0.0, 0.5]), 2.0: np.array([1.0, 1.2]), 3.0: np.array([2.0, 2.2])}
    assert_raises(ValueError, LocalConceptualDFT, values, 'rational')
    # check in valid points
    fname = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    assert_raises(ValueError, LocalConceptualDFT.from_file, fname, 'linear', np.array([0., 0., 0.]))
    assert_raises(ValueError, LocalConceptualDFT.from_file, fname, 'linear', np.array([[0., 0.]]))
    # check molecule file inconsistency
    points = np.array([[0., 0., 0.]])
    fnames = [context.get_fn('test/ch4_uhf_ccpvdz.fchk'), context.get_fn('test/o2_uhf.fchk')]
    assert_raises(ValueError, LocalConceptualDFT.from_file, fnames, 'linear', points)
    fname = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    assert_raises(ValueError, LocalConceptualDFT.from_file, [fname, fname], 'linear', points)
    # check invalid grid
    fnames = [context.get_fn('examples/ch2o_q+0_ub3lyp_augccpvtz.fchk'),
              context.get_fn('examples/ch2o_q+1_ub3lyp_augccpvtz.fchk'),
              context.get_fn('examples/ch2o_q-1_ub3lyp_augccpvtz.fchk')]


def check_local_linear_fmo_ch4_uhf_ccpvdz(filename):
    """Check expected linear local indicators for ch4_uhf_ccpvdz within FMO approach."""
    # Check softness & hyper-softness
    # make molecular grid
    mol = IOData.from_file(context.get_fn(filename))
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='exp:5e-4:2e1:175:434', random_rotate=False, mode='keep')
    # build local conceptual DFT tool
    desp = LocalConceptualDFT.from_file(context.get_fn(filename), 'linear', points=grid.points)
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # check integral of density
    assert_almost_equal(grid.integrate(desp.density(10.)), 10., decimal=4)
    assert_almost_equal(grid.integrate(desp.density(11.)), 11., decimal=4)
    assert_almost_equal(grid.integrate(desp.density(9.)), 9., decimal=4)
    # check shape of Fukui functions & dual descriptor
    assert_equal(desp.ff_zero.shape, grid.shape)
    assert_equal(desp.ff_plus.shape, grid.shape)
    assert_equal(desp.ff_minus.shape, grid.shape)
    # check Fukui functions & dual descriptor
    assert_almost_equal(grid.integrate(desp.ff_plus), 1., decimal=4)
    assert_almost_equal(grid.integrate(desp.ff_minus), 1., decimal=4)
    assert_almost_equal(grid.integrate(desp.ff_zero), 1., decimal=4)
    # check dual descriptor


def test_local_linear_fmo_ch4_uhf_ccpvdz_fchk():
    check_local_linear_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.fchk')


def test_local_linear_fmo_ch4_uhf_ccpvdz_wfn():
    check_local_linear_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.wfn')


def check_local_quadratic_fmo_ch4_uhf_ccpvdz(filename):
    """Check expected quadratic local indicators for ch4_uhf_ccpvdz within FMO approach."""
    # ip = -E(homo) & ea = E(lumo)
    ip, ea = -(-5.43101269E-01), -1.93295185E-01
    eta = ip - ea
    # make molecular grid
    mol = IOData.from_file(context.get_fn(filename))
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                        agspec='exp:5e-4:2e1:175:434', random_rotate=False, mode='keep')
    # build global conceptual DFT tool
    desp = LocalConceptualDFT.from_file(context.get_fn(filename), 'quadratic', points=grid.points)
    # check print statement
    assert_equal(type(desp.__repr__()), str)
    # check shape of density
    assert_equal(desp.density_zero.shape, grid.shape)
    assert_equal(desp.density_plus.shape, grid.shape)
    assert_equal(desp.density_minus.shape, grid.shape)
    # check integral of density
    assert_almost_equal(grid.integrate(desp.density(10.)), 10., decimal=4)
    assert_almost_equal(grid.integrate(desp.density(11.)), 11., decimal=4)
    assert_almost_equal(grid.integrate(desp.density(9.)), 9., decimal=4)
    assert_almost_equal(grid.integrate(desp.density(10.62)), 10.62, decimal=4)
    assert_almost_equal(grid.integrate(desp.density(9.78)), 9.78, decimal=4)
    assert_almost_equal(grid.integrate(desp.density(10.0)), 10.0, decimal=4)
    # Check Fukui function, dual descriptor & softness
    assert_almost_equal(grid.integrate(desp.fukui_function(10.)), 1., decimal=4)
    assert_almost_equal(grid.integrate(desp.fukui_function(10.5)), 1., decimal=4)
    assert_almost_equal(grid.integrate(desp.fukui_function(9.50)), 1., decimal=4)
    assert_almost_equal(grid.integrate(desp.dual_descriptor()), 0., decimal=4)
    # Check local softness
    assert_almost_equal(grid.integrate(desp.softness(10.0, 1./eta)), 1./eta, decimal=4)
    assert_almost_equal(grid.integrate(desp.softness(10.3, 1./eta)), 1./eta, decimal=4)
    assert_almost_equal(grid.integrate(desp.softness(9.10, 1./eta)), 1./eta, decimal=4)
    assert_almost_equal(grid.integrate(desp.hyper_softness(eta)), 0., decimal=3)


def test_local_quadratic_fmo_ch4_uhf_ccpvdz_fchk():
    check_local_quadratic_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.fchk')


def test_local_quadratic_fmo_ch4_uhf_ccpvdz_wfn():
    check_local_quadratic_fmo_ch4_uhf_ccpvdz('test/ch4_uhf_ccpvdz.wfn')
