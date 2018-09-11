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
# pragma pylint: disable=invalid-name,bad-whitespace
"""Test chemtools.analysis.conceptual.LocalConceptualDFT."""


import numpy as np

from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import BeckeMolGrid
from chemtools import context
from chemtools.toolbox.molecule import make_molecule
from chemtools.toolbox.conceptual import LocalConceptualDFT


def get_data_ch4():
    coord = np.array([[-3.77945227e-05,  3.77945227e-05, -1.88972613e-05],
                      [ 1.04290206,  1.50497789,  0.934507367],
                      [ 1.28607202, -1.53098052, -0.477307027],
                      [-1.46467003, -0.702997019,  1.25954026],
                      [-0.864474117, 0.729131931, -1.71670281]])
    nums = np.array([6, 1, 1, 1, 1])
    return coord, nums


def test_local_conceptual_raises():
    # check invalid densities
    values = {1.0: np.array([0.0, 0.5]), 2.0: np.array([1.0]), 3.0: np.array([[2.0]])}
    assert_raises(ValueError, LocalConceptualDFT, values)
    # check in valid model
    values = {1.0: np.array([0.0, 0.5]), 2.0: np.array([1.0, 1.2]), 3.0: np.array([2.0, 2.2])}
    assert_raises(ValueError, LocalConceptualDFT, values, "rational")
    assert_raises(ValueError, LocalConceptualDFT, values, "Rational")
    # check in valid points
    fn = context.get_fn("test/ch4_uhf_ccpvdz.fchk")
    assert_raises(ValueError, LocalConceptualDFT.from_file, fn, "linear", np.array([0., 0., 0.]))
    assert_raises(ValueError, LocalConceptualDFT.from_file, fn, "linear", np.array([[0., 0.]]))
    assert_raises(ValueError, LocalConceptualDFT.from_file, fn, "quadratic", np.array([[0., 0.]]))
    # check molecule file inconsistency
    points = np.array([[0., 0., 0.]])
    fns = [context.get_fn("test/ch4_uhf_ccpvdz.fchk"), context.get_fn("test/o2_uhf.fchk")]
    assert_raises(ValueError, LocalConceptualDFT.from_file, fns, "linear", points)
    assert_raises(ValueError, LocalConceptualDFT.from_file, fns, "quadratic", points)
    fn = context.get_fn("test/ch4_uhf_ccpvdz.fchk")
    assert_raises(ValueError, LocalConceptualDFT.from_file, [fn, fn], "linear", points)
    assert_raises(ValueError, LocalConceptualDFT.from_file, [fn, fn], "quadratic", points)
    # check invalid files
    fns = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
           context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk"),
           context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
    assert_raises(ValueError, LocalConceptualDFT.from_file, fns, "linear", points)
    assert_raises(ValueError, LocalConceptualDFT.from_file, fns, "quadratic", points)


def check_local_reactivity(model, energy_model, grid, n0, eta):
    """Check expected linear local reactivity descriptors."""
    # build local conceptual DFT tool

    # check print statement
    assert_equal(type(model.__repr__()), str)
    # check integral of density
    assert_almost_equal(grid.integrate(model.density(n0 - 1)), n0 - 1, decimal=4)
    assert_almost_equal(grid.integrate(model.density(n0)), n0, decimal=4)
    assert_almost_equal(grid.integrate(model.density(n0 + 1)), n0 + 1, decimal=4)
    assert_almost_equal(grid.integrate(model.density(0.75 * n0)), 0.75 * n0, decimal=4)
    assert_almost_equal(grid.integrate(model.density(1.25 * n0)), 1.25 * n0, decimal=4)
    # Check Fukui function, dual descriptor & softness
    assert_almost_equal(grid.integrate(model.fukui_function), 1., decimal=4)
    assert_almost_equal(grid.integrate(model.density_derivative(n0, 1)), 1., decimal=4)
    assert_almost_equal(grid.integrate(model.density_derivative(n0 - 1, 1)), 1., decimal=4)
    assert_almost_equal(grid.integrate(model.density_derivative(n0 + 1, 1)), 1., decimal=4)
    assert_almost_equal(grid.integrate(model.density_derivative(1.5 * n0, 1)), 1., decimal=4)
    assert_almost_equal(grid.integrate(model.density_derivative(0.8 * n0, 1)), 1., decimal=4)
    if energy_model == "linear":
        # check shape of Fukui functions & dual descriptor
        assert_equal(model.ff_zero.shape, grid.shape)
        assert_equal(model.ff_plus.shape, grid.shape)
        assert_equal(model.ff_minus.shape, grid.shape)
        # check Fukui functions & dual descriptor
        assert_almost_equal(grid.integrate(model.ff_plus), 1., decimal=4)
        assert_almost_equal(grid.integrate(model.ff_minus), 1., decimal=4)
        assert_almost_equal(grid.integrate(model.ff_zero), 1., decimal=4)
    if energy_model == "quadratic":
        assert_almost_equal(grid.integrate(model.dual_descriptor), 0., decimal=4)
        # assert_almost_equal(grid.integrate(model.softness(1./eta)), 1./eta, decimal=4)
        # assert_almost_equal(grid.integrate(model.softness(1./eta, 10.3)), 1./eta, decimal=4)
        # assert_almost_equal(grid.integrate(model.softness(1./eta, 9.1)), 1./eta, decimal=4)
        # assert_almost_equal(grid.integrate(model.hyper_softness(eta)), 0., decimal=3)


def test_local_linear_from_file_fmo_ch4_uhf_ccpvdz_fchk():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    file_path = context.get_fn("test/ch4_uhf_ccpvdz.fchk")
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_file passing a grid
    model = LocalConceptualDFT.from_file(file_path, "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)
    # check from_file given as a list passing a grid
    model = LocalConceptualDFT.from_file([file_path], "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)


def test_local_linear_from_molecule_fmo_ch4_uhf_ccpvdz_fchk():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    molecule = make_molecule(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_molecule passing a grid
    model = LocalConceptualDFT.from_molecule(molecule, "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)
    # check from_molecule given as a list passing a grid
    model = LocalConceptualDFT.from_molecule([molecule], "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)


def test_local_linear_from_file_fmo_ch4_uhf_ccpvdz_wfn():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    file_path = context.get_fn("test/ch4_uhf_ccpvdz.wfn")
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_file passing a grid
    model = LocalConceptualDFT.from_file(file_path, "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)
    # check from_file given as a list passing a grid
    model = LocalConceptualDFT.from_file([file_path], "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)


def test_local_linear_from_molecule_fmo_ch4_uhf_ccpvdz_wfn():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    molecule = make_molecule(context.get_fn("test/ch4_uhf_ccpvdz.wfn"))
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_molecule passing a grid
    model = LocalConceptualDFT.from_molecule(molecule, "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)
    # check from_molecule given as a list passing a grid
    model = LocalConceptualDFT.from_molecule([molecule], "linear", grid.points)
    check_local_reactivity(model, "linear", grid, 10, None)


def test_local_quadratic_from_file_fmo_ch4_uhf_ccpvdz_fchk():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    # ip = -E(homo) & ea = E(lumo) & eta = ip - ea
    eta = -(-5.43101269E-01) - (-1.93295185E-01)
    file_path = context.get_fn("test/ch4_uhf_ccpvdz.fchk")
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_file passing grid
    model = LocalConceptualDFT.from_file(file_path, "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)
    # check from_file given as a list passing grid
    model = LocalConceptualDFT.from_file([file_path], "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)


def test_local_quadratic_from_molecule_fmo_ch4_uhf_ccpvdz_fchk():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    # ip = -E(homo) & ea = E(lumo) & eta = ip - ea
    eta = -(-5.43101269E-01) - (-1.93295185E-01)
    molecule = make_molecule(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_molecule passing grid
    model = LocalConceptualDFT.from_molecule(molecule, "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)
    # check from_molecule given as a list passing grid
    model = LocalConceptualDFT.from_molecule([molecule], "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)


def test_local_quadratic_from_file_fmo_ch4_uhf_ccpvdz_wfn():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    # ip = -E(homo) & ea = E(lumo) & eta = ip - ea
    eta = -(-5.43101269E-01) - (-1.93295185E-01)
    file_path = context.get_fn("test/ch4_uhf_ccpvdz.wfn")
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_file passing grid
    model = LocalConceptualDFT.from_file(file_path, "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)
    # check from_file given as a list passing grid
    model = LocalConceptualDFT.from_file([file_path], "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)


def test_local_quadratic_from_molecule_fmo_ch4_uhf_ccpvdz_wfn():
    # atomic coordinates and numbers of CH4
    coord, nums = get_data_ch4()
    # ip = -E(homo) & ea = E(lumo) & eta = ip - ea
    eta = -(-5.43101269E-01) - (-1.93295185E-01)
    molecule = make_molecule(context.get_fn("test/ch4_uhf_ccpvdz.wfn"))
    grid = BeckeMolGrid(coord, nums, nums, agspec="insane", random_rotate=False, mode="keep")
    # check from_molecule passing grid
    model = LocalConceptualDFT.from_molecule(molecule, "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)
    # check from_molecule given as a list passing grid
    model = LocalConceptualDFT.from_molecule([molecule], "quadratic", grid.points)
    check_local_reactivity(model, "quadratic", grid, 10, eta)
