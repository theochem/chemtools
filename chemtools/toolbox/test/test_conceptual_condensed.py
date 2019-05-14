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
# pragma pylint: disable=invalid-name,bad-whitespace
"""Test chemtools.analysis.conceptual.CondensedConceptualDFT."""


import numpy as np

from numpy.testing import assert_raises, assert_equal, assert_almost_equal

from horton import BeckeMolGrid
from chemtools.toolbox.conceptual import CondensedConceptualDFT
from chemtools.wrappers.molecule import Molecule
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


# def test_condensed_conceptual_raises():
#     # file for FMO
#     with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as file:
#         fname = file
#     # files for FD
#     fnames = []
#     with path('chemtools.data.examples', "ch2o_q+0_ub3lyp_augccpvtz.fchk") as file1:
#         fnames.append(file1)
#     with path('chemtools.data.examples', "ch2o_q+1_ub3lyp_augccpvtz.fchk") as file2:
#         fnames.append(file2)
#     with path('chemtools.data.examples', "ch2o_q-1_ub3lyp_augccpvtz.fchk") as file3:
#         fnames.append(file3)
#     # check invalid scheme
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, "quadratic", scheme="gib")
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, "linear", scheme="err")
#     # check invalid grid type
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, "quadratic", grid="str")
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, "linear", grid="gibb")
#     # check invalid grid coordinates
#     grid = BeckeMolGrid(np.array([[0., 0., 0.]]), np.array([1]))
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, "quadratic", grid=grid)
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, "quadratic", grid=grid)
#     # check invalid grid atomic numbers
#     coord = np.array([[-3.77945227E-05,  3.77945227E-05, -1.88972613E-05],
#                       [ 1.04290206E+00,  1.50497789E+00,  9.34507367E-01],
#                       [ 1.28607202E+00, -1.53098052E+00, -4.77307027E-01],
#                       [-1.46467003E+00, -7.02997019E-01,  1.25954026E+00],
#                       [-8.64474117E-01,  7.29131931E-01, -1.71670281E+00]])
#     grid = BeckeMolGrid(coord, np.array([7, 1, 1, 1, 1]))
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, "quadratic", grid=grid)
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, "quadratic", grid=grid)
#     # check number of molecules
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fname, "linear", "RMF", "esp")
#     # check approach
#     assert_raises(ValueError, CondensedConceptualDFT.from_file, fnames, "linear", "FMR", "esp")


def check_condensed_reactivity(model, energy_model, pop_0, pop_p, pop_m, n0):
    """Check expected condensed reactivity descriptors."""
    # check print statement
    assert_equal(type(model.__repr__()), str)
    # check expected charges
    assert_almost_equal(model.population(n0), pop_0, decimal=3)
    if pop_m is not None:
        assert_almost_equal(model.population(n0 - 1), pop_m, decimal=3)
    if pop_p is not None:
        assert_almost_equal(model.population(n0 + 1), pop_p, decimal=3)
    # check condensed density
    assert_almost_equal(np.sum(model.population(n0 + 1)), n0 + 1, decimal=2)
    assert_almost_equal(np.sum(model.population(n0)), n0, decimal=3)
    assert_almost_equal(np.sum(model.population(n0 - 1)), n0 - 1, decimal=3)
    # check condensed density with arbitrary number of electrons
    assert_almost_equal(np.sum(model.population(0.90 * n0)), 0.90 * n0, decimal=2)
    assert_almost_equal(np.sum(model.population(1.15 * n0)), 1.15 * n0, decimal=2)
    # check condensed fukui function with arbitrary number of electrons
    assert_almost_equal(np.sum(model.fukui_function), 1.0, decimal=2)
    assert_almost_equal(np.sum(model.population_derivative(n0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(model.population_derivative(0.85 * n0, 1)), 1.0, decimal=2)
    assert_almost_equal(np.sum(model.population_derivative(1.20 * n0, 1)), 1.0, decimal=2)
    if energy_model == "linear":
        if pop_m is not None and pop_p is not None:
            assert_almost_equal(model.ff_plus, pop_p - pop_0, decimal=3)
            assert_almost_equal(model.ff_zero, 0.5 * (pop_p - pop_m), decimal=3)
            assert_almost_equal(model.ff_minus, pop_0 - pop_m, decimal=3)
        assert_almost_equal(np.sum(model.ff_plus), 1., decimal=2)
        assert_almost_equal(np.sum(model.ff_zero), 1., decimal=2)
        assert_almost_equal(np.sum(model.ff_minus), 1.0, decimal=2)
    if energy_model == "quadratic":
        # check condensed dual descriptor
        assert_almost_equal(np.sum(model.dual_descriptor), 0.0, decimal=2)


def test_condense_linear_from_file_fmr_h_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "h")
        model2 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "h")
        mol = Molecule.from_file(filename)
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                            agspec="insane", random_rotate=False, mode="keep")
        model3 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "h", grid=grid)
        model4 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "h", grid=grid)
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    # check using filename given as a string
    check_condensed_reactivity(model1, "linear", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "linear", expected, None, None, 10)
    # check using filename as a string & passing grid
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec="insane",
                        random_rotate=False, mode="keep")
    check_condensed_reactivity(model3, "linear", expected, None, None, 10)
    # check using filename as a list & passing grid
    check_condensed_reactivity(model4, "linear", expected, None, None, 10)


def test_condense_linear_from_molecule_fmr_h_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    # check from_molecule
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "h")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule & passing grid
    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                        agspec="insane", random_rotate=False, mode="keep")
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list & passing grid
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "h", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)


def test_condense_linear_from_file_fmr_h_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "h")
        model2 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "h")
        mol = Molecule.from_file(filename)
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec="insane",
                            random_rotate=False, mode="keep")
        model3 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "h", grid=grid)
        model4 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "h", grid=grid)
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    # check using filename given as a string
    check_condensed_reactivity(model1, "linear", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "linear", expected, None, None, 10)
    # check using filename as a string & passing grid
    mol = Molecule.from_file(filename)
    grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec="insane",
                        random_rotate=False, mode="keep")
    check_condensed_reactivity(model3, "linear", expected, None, None, 10)
    # check using filename as a list & passing grid
    check_condensed_reactivity(model4, "linear", expected, None, None, 10)


def test_condense_linear_from_molecule_fmr_h_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.11301651, 0.97175462, 0.97175263, 0.9717521, 0.97174353])
    # check from_molecule
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "h")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule & passing grid
    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                        agspec="insane", random_rotate=False, mode="keep")
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list & passing grid
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "h", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)


def test_condense_linear_from_file_fmr_mbis_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "mbis")
        model2 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "mbis")
        mol = Molecule.from_file(filename)
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec="insane",
                            random_rotate=False, mode="keep")
        model3 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "mbis", grid=grid)
        model4 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "mbis", grid=grid)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check using filename given as a string
    check_condensed_reactivity(model1, "linear", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "linear", expected, None, None, 10)
    # check using filename as a string & passing grid
    check_condensed_reactivity(model3, "linear", expected, None, None, 10)
    # check using filename as a list & passing grid
    check_condensed_reactivity(model4, "linear", expected, None, None, 10)


def test_condense_linear_from_molecule_fmr_mbis_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check from_molecule
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "mbis")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule & passing grid
    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                        agspec="insane", random_rotate=False, mode="keep")
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list & passing grid
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)


def test_condense_linear_from_file_fmr_mbis_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "linear", "FMR", "mbis")
        model2 = CondensedConceptualDFT.from_file([filename], "linear", "FMR", "mbis")
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check using filename given as a string
    check_condensed_reactivity(model1, "linear", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "linear", expected, None, None, 10)


def test_condense_linear_from_molecule_fmr_mbis_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check from_molecule
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "mbis")
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule & passing grid
    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                        agspec="insane", random_rotate=False, mode="keep")
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)
    # check from_molecule given as a list & passing grid
    model = CondensedConceptualDFT.from_molecule([molecule], "linear", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "linear", expected, None, None, 10)


def test_condense_quadratic_from_file_fmr_mbis_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "quadratic", "FMR", "mbis")
        model2 = CondensedConceptualDFT.from_file([filename], "quadratic", "FMR", "mbis")
        mol = Molecule.from_file(filename)
        grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, agspec="insane",
                            random_rotate=False, mode="keep")
        model3 = CondensedConceptualDFT.from_file(filename, "quadratic", "FMR", "mbis", grid=grid)
        model4 = CondensedConceptualDFT.from_file([filename], "quadratic", "FMR", "mbis", grid=grid)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check using filename given as a string
    check_condensed_reactivity(model1, "quadratic", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "quadratic", expected, None, None, 10)
    # check using filename as a string & passing grid
    check_condensed_reactivity(model3, "quadratic", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model4, "quadratic", expected, None, None, 10)


def test_condense_quadratic_from_molecule_fmr_mbis_ch4_fchk():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check from_molecule
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "FMR", "mbis")
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "quadratic", "FMR", "mbis")
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)
    # check from_molecule & passing grid
    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                        agspec="insane", random_rotate=False, mode="keep")
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)
    # check from_molecule given as a list & passing grid
    model = CondensedConceptualDFT.from_molecule([molecule], "quadratic", "FMR", "mbis", grid=grid)
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)


def test_condense_quadratic_from_file_fmr_mbis_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        model1 = CondensedConceptualDFT.from_file(filename, "quadratic", "FMR", "mbis")
        model2 = CondensedConceptualDFT.from_file([filename], "quadratic", "FMR", "mbis")
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check using filename given as a string
    check_condensed_reactivity(model1, "quadratic", expected, None, None, 10)
    # check using filename given as a list
    check_condensed_reactivity(model2, "quadratic", expected, None, None, 10)


def test_condense_quadratic_from_molecule_fmr_mbis_ch4_wfn():
    # expected populations of CH4 computed with HORTON
    with path('chemtools.data', 'ch4_uhf_ccpvdz.wfn') as filename:
        molecule = Molecule.from_file(filename)
    expected = np.array([6.46038055, 0.88489494, 0.88492901, 0.88493897, 0.88492396])
    # check from_molecule given as a string
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "FMR", "mbis")
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)
    # check from_molecule given as a list
    model = CondensedConceptualDFT.from_molecule([molecule], "quadratic", "FMR", "mbis")
    check_condensed_reactivity(model, "quadratic", expected, None, None, 10)


def test_condense_from_file_fd_rmf_esp_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([ 4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-5.81613550E-01,-2.09193820E-01,-2.09192630E-01])
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
                filename = [file1, file2, file3]
                # check from_file linear
                model1 = CondensedConceptualDFT.from_file(filename, "linear", "RMF", "esp")
                model2 = CondensedConceptualDFT.from_file(filename, "quadratic", "RMF", "esp")
    check_condensed_reactivity(model1, "linear", expected_0, expected_p, expected_m, 10)
    # check from_file quadratic
    check_condensed_reactivity(model2, "quadratic", expected_0, expected_p, expected_m, 10)


def test_condense_from_molecule_fd_rmf_esp_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([ 4.14233893E-02, 4.79288419E-01, 4.79288192E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-7.00779373E-01, 3.50389629E-01, 3.50389744E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-5.81613550E-01,-2.09193820E-01,-2.09192630E-01])
    molecule = []
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        molecule.append(Molecule.from_file(file1))
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
        molecule.append(Molecule.from_file(file2))
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
        molecule.append(Molecule.from_file(file3))
    # check from_molecule linear
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "RMF", "esp")
    check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 10)
    # check from_molecule quadratic
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "RMF", "esp")
    check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 10)


def test_condense_from_file_fd_rmf_npa_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([-3.64452391E-02, 5.18222784E-01, 5.18222455E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-9.00876494E-01, 4.50438267E-01, 4.50438227E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-1.11332869E+00, 5.66635486E-02, 5.66651430E-02])
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
                filename = [file1, file2, file3]
                # check from_file linear
                model1 = CondensedConceptualDFT.from_file(filename, "linear", "RMF", "npa")
                model2 = CondensedConceptualDFT.from_file(filename, "quadratic", "RMF", "npa")
    check_condensed_reactivity(model1, "linear", expected_0, expected_p, expected_m, 10)
    # check from_file quadratic
    check_condensed_reactivity(model2, "quadratic", expected_0, expected_p, expected_m, 10)


def test_condense_from_molecule_fd_rmf_npa_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([-3.64452391E-02, 5.18222784E-01, 5.18222455E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-9.00876494E-01, 4.50438267E-01, 4.50438227E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-1.11332869E+00, 5.66635486E-02, 5.66651430E-02])
    molecule = []
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        molecule.append(Molecule.from_file(file1))
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
        molecule.append(Molecule.from_file(file2))
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
        molecule.append(Molecule.from_file(file3))
    # check from_molecule linear
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "RMF", "npa")
    check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 10)
    # check from_molecule quadratic
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "RMF", "npa")
    check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 10)


def test_condense_from_file_fd_rmf_mulliken_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([ 3.49417097E-01, 3.25291762E-01, 3.25291141E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-4.32227787E-01, 2.16114060E-01, 2.16113727E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-2.64833827E-01,-3.67583325E-01,-3.67582849E-01])
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
            with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
                filename = [file1, file2, file3]
                model1 = CondensedConceptualDFT.from_file(filename, "linear", "RMF", "mulliken")
                model2 = CondensedConceptualDFT.from_file(filename, "quadratic", "RMF", "mulliken")
    # check from_file linear
    check_condensed_reactivity(model1, "linear", expected_0, expected_p, expected_m, 10)
    # check from_file quadratic
    check_condensed_reactivity(model2, "quadratic", expected_0, expected_p, expected_m, 10)


def test_condense_from_molecule_fd_rmf_mulliken_h2o_fchk():
    # expected populations of H2O from fchk file
    expected_m = np.array([8, 1, 1]) - np.array([ 3.49417097E-01, 3.25291762E-01, 3.25291141E-01])
    expected_0 = np.array([8, 1, 1]) - np.array([-4.32227787E-01, 2.16114060E-01, 2.16113727E-01])
    expected_p = np.array([8, 1, 1]) - np.array([-2.64833827E-01,-3.67583325E-01,-3.67582849E-01])
    molecule = []
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as file1:
        molecule.append(Molecule.from_file(file1))
    with path('chemtools.data', 'h2o_q-1_ub3lyp_ccpvtz.fchk') as file2:
        molecule.append(Molecule.from_file(file2))
    with path('chemtools.data', 'h2o_q+1_ub3lyp_ccpvtz.fchk') as file3:
        molecule.append(Molecule.from_file(file3))
    # check from_molecule linear
    model = CondensedConceptualDFT.from_molecule(molecule, "linear", "RMF", "mulliken")
    check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 10)
    # check from_molecule quadratic
    model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "RMF", "mulliken")
    check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 10)


# def test_condense_linear_from_file_fd_rmf_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     file_path = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
#     # check from_file
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "RMF", "h")
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)
#     # check from_file passing a grid
#     mol = Molecule.from_file(file_path[0])
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "RMF", "h", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)


# def test_condense_linear_from_molecule_fd_rmf_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     molecule = [Molecule.from_file(context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk"))]
#     # check from_molecule
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "RMF", "h")
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)
#     # check from_molecule passing as a grid
#     grid = BeckeMolGrid(molecule[0].coordinates, molecule[0].numbers, molecule[0].pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "RMF", "h", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)


# def test_condense_linear_from_file_fd_fmr_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     file_path = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
#     # check from_file
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "FMR", "h")
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)
#     # check from_file passing as a grid
#     mol = Molecule.from_file(file_path[0])
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "FMR", "h", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)


# def test_condense_linear_from_molecule_fd_fmr_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     molecule = [Molecule.from_file(context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk"))]
#     # check from_file
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h")
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)
#     # check from_file passing as a grid
#     grid = BeckeMolGrid(molecule[0].coordinates, molecule[0].numbers, molecule[0].pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "h", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)


# def test_condense_linear_from_file_fd_rmf_mbis_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.8580338, 5.70425809, 0.71885554, 0.71884404])
#     expected_0 = np.array([8.41149, 5.66445074, 0.96204946, 0.96202722])
#     expected_p = np.array([8.13881352, 6.81770852, 1.28123219, 0.76225513])
#     file_path = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
#     # check from_file
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "RMF", "mbis")
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)
#     # check from_file passing grid
#     mol = Molecule.from_file(file_path[0])
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "RMF", "mbis", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, expected_p, expected_m, 16)


# def test_condense_linear_from_file_fd_fmr_mbis_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_0 = np.array([8.41149, 5.66445074, 0.96204946, 0.96202722])
#     file_path = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
#     # check from_file
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "FMR", "mbis")
#     check_condensed_reactivity(model, "linear", expected_0, None, None, 16)
#     # check from_file passing grid
#     mol = Molecule.from_file(file_path[0])
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_file(file_path, "linear", "FMR", "mbis", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, None, None, 16)


# def test_condense_linear_from_molecule_fd_fmr_mbis_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_0 = np.array([8.41149, 5.66445074, 0.96204946, 0.96202722])
#     molecule = [Molecule.from_file(context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk"))]
#     # check from_molecule
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis")
#     check_condensed_reactivity(model, "linear", expected_0, None, None, 16)
#     # check from_molecule passing grid
#     grid = BeckeMolGrid(molecule[0].coordinates, molecule[0].numbers, molecule[0].pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_molecule(molecule, "linear", "FMR", "mbis", grid=grid)
#     check_condensed_reactivity(model, "linear", expected_0, None, None, 16)


# def test_condense_quadratic_from_file_fd_rmf_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     file_path = [context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk"),
#                  context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk")]
#     # check from_file
#     model = CondensedConceptualDFT.from_file(file_path, "quadratic", "RMF", "h")
#     check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 16)
#     # check from_file passing grid
#     mol = Molecule.from_file(file_path[0])
#     grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_file(file_path, "quadratic", "RMF", "h", grid=grid)
#     check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 16)


# def test_condense_quadratic_from_molecule_fd_rmf_h_ch2o_fchk():
#     # expected populations of CH2O computed with HORTON
#     expected_m = np.array([7.98237872, 5.47698573, 0.77030456, 0.77031781])
#     expected_0 = np.array([8.46718639, 5.67524299, 0.92860658, 0.92866436])
#     expected_p = np.array([8.76534627, 6.18498153, 1.02517556, 1.02513059])
#     molecule = [Molecule.from_file(context.get_fn("examples/ch2o_q+0_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q+1_ub3lyp_augccpvtz.fchk")),
#                 Molecule.from_file(context.get_fn("examples/ch2o_q-1_ub3lyp_augccpvtz.fchk"))]
#     # check from_file
#     model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "RMF", "h")
#     check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 16)
#     # check from_file passing grid
#     grid = BeckeMolGrid(molecule[0].coordinates, molecule[0].numbers, molecule[0].pseudo_numbers,
#                         agspec="insane", random_rotate=False, mode="keep")
#     model = CondensedConceptualDFT.from_molecule(molecule, "quadratic", "RMF", "h", grid=grid)
#     check_condensed_reactivity(model, "quadratic", expected_0, expected_p, expected_m, 16)
