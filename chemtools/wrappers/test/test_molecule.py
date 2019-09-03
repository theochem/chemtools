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
"""Test chemtools.wrappers.molecule."""


import numpy as np
from numpy.testing import assert_raises, assert_equal, assert_almost_equal
from chemtools.wrappers import Molecule

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def check_molecule_raises(mol):
    """Check expected raised error messages by HortonWaveFunction class."""
    # example point array
    points = np.array([[0., 0., 0.], [1., 1., 1.]])
    # check invalid orbital spin argument
    assert_raises(ValueError, mol.compute_density_matrix, "alphabeta")
    assert_raises(ValueError, mol.compute_density, points, spin="alpha")
    assert_raises(ValueError, mol.compute_gradient, points, spin="beta")
    assert_raises(ValueError, mol.compute_hessian, points, spin="betaalpha")
    assert_raises(ValueError, mol.compute_ked, points, spin="balpha")
    # check invalid points argument
    assert_raises(ValueError, mol.compute_molecular_orbital, [0.1, 0.5, 0.7], "a")
    assert_raises(ValueError, mol.compute_molecular_orbital, np.array([0.1, 0.5, 0.7]), "b")
    assert_raises(ValueError, mol.compute_molecular_orbital, np.array([[0.4, 0.2]]), "a")
    assert_raises(ValueError, mol.compute_molecular_orbital, np.array([[5, 10, 15]]), "b")
    assert_raises(ValueError, mol.compute_density, [0., 0., 0.])
    assert_raises(ValueError, mol.compute_density, np.array([0., 0., 0.]))
    assert_raises(ValueError, mol.compute_density, np.array([[0., 0.]]))
    assert_raises(ValueError, mol.compute_density, np.array([[0, 0, 0]]))
    assert_raises(ValueError, mol.compute_gradient, [.1, .2, .3])
    assert_raises(ValueError, mol.compute_gradient, np.array([.1, .2, .3]))
    assert_raises(ValueError, mol.compute_gradient, np.array([[.1, 2., .3, .4]]))
    assert_raises(ValueError, mol.compute_gradient, np.array([[1, 2, 3]]))
    assert_raises(ValueError, mol.compute_hessian, [.5, .5, .5])
    assert_raises(ValueError, mol.compute_hessian, np.array([.5, .5, .5]))
    assert_raises(ValueError, mol.compute_hessian, np.array([[.5, 5.]]))
    assert_raises(ValueError, mol.compute_hessian, np.array([[5, 5, 5]]))
    assert_raises(ValueError, mol.compute_esp, [1., .5])
    assert_raises(ValueError, mol.compute_esp, np.array([1., .5, .25]))
    assert_raises(ValueError, mol.compute_esp, np.array([[1., .25]]))
    assert_raises(ValueError, mol.compute_esp, np.array([[1, 25, 10]]))
    assert_raises(ValueError, mol.compute_ked, [.5, 0., .2])
    assert_raises(ValueError, mol.compute_ked, np.array([.5, 0., .2]))
    assert_raises(ValueError, mol.compute_ked, np.array([[.5, 0., .2, .1, .3]]))
    assert_raises(ValueError, mol.compute_ked, np.array([[5, 0, 2]]))
    # check invalid charges argument
    assert_raises(ValueError, mol.compute_esp, points, charges=np.array([6., 1., 1.]))
    assert_raises(ValueError, mol.compute_esp, points, charges=[6., 1., 1., 1., 1.])


def test_molecule_raises_fchk_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_raises(molecule)


def test_molecule_raises_fchk_rhf_ch4():
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_raises(molecule)


def test_molecule_raises_wfn_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.wfn") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_raises(molecule)


def check_molecule_basics(mol):
    """Check expected basic attributes of HortonWaveFunction class."""
    # check basic numerics
    assert_equal(mol.natom, 5)
    assert_equal(mol.nelectrons, (5, 5))
    # assert_equal(mol.ao.nbasis, 34)
    assert_equal(mol.numbers, [6, 1, 1, 1, 1])
    assert_equal(mol.pseudo_numbers, [6, 1, 1, 1, 1])
    assert_equal(mol.mo.homo_index, (5, 5))
    assert_equal(mol.mo.lumo_index, (6, 6))
    assert_equal(mol.mo.occupation[0], np.array([1] * 5 + [0] * 29))
    assert_equal(mol.mo.occupation[1], np.array([1] * 5 + [0] * 29))
    assert_almost_equal(mol.energy, -4.019868797400735E+01, decimal=8)
    # check coordinates
    coord = np.array([[-3.77945227E-05,  3.77945227E-05, -1.88972613E-05],
                      [ 1.04290206E+00,  1.50497789E+00,  9.34507367E-01],
                      [ 1.28607202E+00, -1.53098052E+00, -4.77307027E-01],
                      [-1.46467003E+00, -7.02997019E-01,  1.25954026E+00],
                      [-8.64474117E-01,  7.29131931E-01, -1.71670281E+00]])
    assert_almost_equal(mol.coordinates, coord, decimal=6)


def test_molecule_basics_fchk_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    # check basics
    check_molecule_basics(molecule)
    # check charges
    esp = np.array([-0.502277518, 0.125567970, 0.125569655, 0.125566743, 0.125573150])
    assert_almost_equal(molecule.esp_charges, esp, decimal=6)
    npa = np.array([-0.791299253, 0.197824989, 0.197825250, 0.197824326, 0.197824689])
    assert_almost_equal(molecule.npa_charges, npa, decimal=6)
    mul = np.array([-0.139702704, 0.0349253868, 0.0349266071, 0.0349235395, 0.0349271707])
    assert_almost_equal(molecule.mulliken_charges, mul, decimal=6)


def test_molecule_basics_wfn_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.wfn") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_basics(molecule)


def test_molecule_orbitals_fchk_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        mol = Molecule.from_file(fname)
    # check orbital energy
    orb_energy = np.array([-1.12152085E+01, -9.42914385E-01, -5.43117091E-01, -5.43114279E-01,
                           -5.43101269E-01,  1.93295185E-01,  2.74358942E-01,  2.74359310E-01,
                            2.74359740E-01,  5.89328697E-01,  5.89342443E-01,  5.89343893E-01,
                            8.90214386E-01,  8.90219069E-01,  8.90222480E-01,  9.36275219E-01,
                            1.13182813E+00,  1.13184117E+00,  1.25675685E+00,  1.68763897E+00,
                            1.68764372E+00,  1.68765502E+00,  1.89570058E+00,  1.89570452E+00,
                            1.89571385E+00,  2.21323213E+00,  2.21324619E+00,  2.21328532E+00,
                            2.54691042E+00,  2.54694190E+00,  2.75532231E+00,  2.79906776E+00,
                            2.79907762E+00,  2.79908651E+00])
    assert_almost_equal(mol.mo.energy[0], orb_energy, decimal=6)
    assert_almost_equal(mol.mo.energy[1], orb_energy, decimal=6)
    assert_almost_equal(mol.mo.homo_energy[0], orb_energy[4], decimal=6)
    assert_almost_equal(mol.mo.homo_energy[1], orb_energy[4], decimal=6)
    assert_almost_equal(mol.mo.lumo_energy[0], orb_energy[5], decimal=6)
    assert_almost_equal(mol.mo.lumo_energy[1], orb_energy[5], decimal=6)
    # check orbital coefficients
    orb_coeff = np.array([9.97287609E-01, 1.86004593E-02, -8.24772487E-03])
    assert_almost_equal(mol.mo.coefficient[1][:3, 0], orb_coeff, decimal=6)
    assert_almost_equal(mol.mo.coefficient[0][:3, 0], orb_coeff, decimal=6)
    assert_almost_equal(mol.mo.coefficient[1][0, 1], -0.188285003, decimal=6)
    assert_almost_equal(mol.mo.coefficient[0][0, 1], -0.188285003, decimal=6)
    assert_almost_equal(mol.mo.coefficient[1][-1, -1], 1.02960200, decimal=6)
    assert_almost_equal(mol.mo.coefficient[0][-1, -1], 1.02960200, decimal=6)
    # check overlap matrix
    overlap = mol.compute_orbital_overlap()
    assert_equal(overlap.shape, (34, 34))
    assert_almost_equal(np.diag(overlap), np.ones(34), decimal=6)
    assert_almost_equal(overlap, overlap.T, decimal=6)


def test_molecule_density_matrix_fchk_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        mol = Molecule.from_file(fname)
    # check density matrix against Gaussian (printed in log fname)
    expected_diag = np.array([1.03003, 0.13222, 0.05565, 0.17944, 0.17944, 0.17944, 0.03891,
                              0.03891, 0.03891, 0.00028, 0.00064, 0.00045, 0.00011, 0.00072,
                              0.18575, 0.02710, 0.00044, 0.00072, 0.00038, 0.18575, 0.02710,
                              0.00057, 0.00074, 0.00023, 0.18575, 0.02710, 0.00069, 0.00029,
                              0.00056, 0.18575, 0.02710, 0.00035, 0.00030, 0.00089])
    # column 5, rows 15-21
    expected_05 = np.array([0.12066, 0.05027, -0.00560, -0.00252, -0.00502, -0.12275, -0.05114])
    # column 15, rows 16-34
    expected_15 = np.array([ 0.06838, -0.00691, -0.00996, -0.00619, -0.01612, -0.01573, 0.00244,
                             0.00393,  0.00238, -0.01612, -0.01573,  0.00277,  0.00383, 0.00217,
                            -0.01612, -0.01573,  0.00270,  0.00366,  0.00253])
    # column 19, rows 20-34
    expected_19 = np.array([-0.00130,  0.00004,  0.00004, -0.00033,  0.00002,  0.00302, 0.00184,
                             0.00001, -0.00010, -0.00003, -0.00438, -0.00125, -0.00032, 0.00003,
                            -0.00035])
    # column 29, rows 30-34
    expected_29 = np.array([-0.00442, -0.00106, -0.00003, 0.00029, -0.00047])
    # check alpha density matrix
    dm_array_a = mol.compute_density_matrix(spin="a")
    assert_almost_equal(np.diag(dm_array_a), expected_diag, decimal=5)
    assert_almost_equal(dm_array_a, dm_array_a.T, decimal=5)
    assert_almost_equal(dm_array_a[0, 1:3], np.array([-0.04982, -0.05262]), decimal=5)
    assert_almost_equal(dm_array_a[4, 14:21], expected_05, decimal=5)
    assert_almost_equal(dm_array_a[14, 15:], expected_15, decimal=5)
    assert_almost_equal(dm_array_a[18, 19:], expected_19, decimal=5)
    assert_almost_equal(dm_array_a[28, 29:], expected_29, decimal=5)
    # check beta density matrix
    dm_array_b = mol.compute_density_matrix(spin="b")
    assert_almost_equal(np.diag(dm_array_b), expected_diag, decimal=5)
    assert_almost_equal(dm_array_b, dm_array_b.T, decimal=5)
    assert_almost_equal(dm_array_b[0, 1:3], np.array([-0.04982, -0.05262]), decimal=5)
    assert_almost_equal(dm_array_b[4, 14:21], expected_05, decimal=5)
    assert_almost_equal(dm_array_b[14, 15:], expected_15, decimal=5)
    assert_almost_equal(dm_array_b[18, 19:], expected_19, decimal=5)
    assert_almost_equal(dm_array_b[28, 29:], expected_29, decimal=5)
    # check total density matrix
    dm_array_ab = mol.compute_density_matrix(spin="ab")
    assert_almost_equal(np.diag(dm_array_ab), 2 * expected_diag, decimal=5)
    assert_almost_equal(dm_array_ab, dm_array_ab.T, decimal=5)
    assert_almost_equal(dm_array_ab[0, 1:3], 2*np.array([-0.04982, -0.05262]), decimal=5)
    assert_almost_equal(dm_array_ab[4, 14:21], 2*expected_05, decimal=5)
    assert_almost_equal(dm_array_ab[14, 15:], 2*expected_15, decimal=5)
    assert_almost_equal(dm_array_ab[18, 19:], 2*expected_19, decimal=5)
    assert_almost_equal(dm_array_ab[28, 29:], 2*expected_29, decimal=5)


def test_molecule_esp_fchk_uhf_ch4():
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        mol = Molecule.from_file(fname)
    # check esp against Gaussian (printed in log file)
    # check esp at the position of each nuclei (1.e-14 is added to avoid division by zero)
    # excluding the nucleus itself.
    point = mol.coordinates[0].reshape(1, 3) + 1.e-14
    charge = np.array([0., 1., 1., 1., 1.])
    assert_almost_equal(mol.compute_esp(point, charges=charge), [-14.745629], decimal=5)
    point = mol.coordinates[1].reshape(1, 3) + 1.e-14
    charge = np.array([6., 0., 1., 1., 1.])
    assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
    point = mol.coordinates[2].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 0., 1., 1.])
    assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
    point = mol.coordinates[3].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 1., 0., 1.])
    assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116067], decimal=5)
    point = mol.coordinates[4].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 1., 1., 0.])
    assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
    # check esp at non-nuclei points
    points = np.array([[ 0.5,  0.5,  0.5],
                       [-0.5, -0.5, -0.5],
                       [-0.5,  0.5,  0.5],
                       [-0.5, -0.5,  0.5],
                       [-0.5,  0.5, -0.5],
                       [ 0.5, -0.5, -0.5],
                       [ 0.5, -0.5,  0.5],
                       [ 0.5,  0.5, -0.5]]) / 0.529177
    expected_esp = np.array([0.895650, 0.237257, 0.234243, 0.708301,
                             0.499083, 0.479275, 0.241434, 0.235102])
    assert_almost_equal(mol.compute_esp(points), expected_esp, decimal=5)


def check_molecule_against_gaussian_ch4(mol):
    """Check local properties of HortonWaveFunction class against Gaussian09_C.01"s cubegen."""
    # get expected data computed by Gaussian09_C.01 cubegen
    with path("chemtools.data", "data_cubegen_g09_C01_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
        points, dens, grad = data["points"], data["dens"], data["grad"]
        lap, hess_xx, esp = data["lap"], data["hess_xx"], data["esp"]
    # check density, gradient, esp & hessian
    assert_almost_equal(mol.compute_density(points, "ab"), dens, decimal=5)
    assert_almost_equal(mol.compute_density(points, "a"), 0.5 * dens, decimal=5)
    assert_almost_equal(mol.compute_density(points, "b"), 0.5 * dens, decimal=5)
    assert_almost_equal(mol.compute_gradient(points, "ab"), grad, decimal=5)
    assert_almost_equal(mol.compute_esp(points, "ab"), esp, decimal=5)
    assert_almost_equal(mol.compute_laplacian(points, "ab", None), lap, decimal=5)
    assert_almost_equal(mol.compute_hessian(points, "ab", None)[:, 0, 0], hess_xx, decimal=5)
    # density computed by summing squared mo expressions
    assert_almost_equal(mol.compute_density(points, "ab", range(1, 6)), dens, decimal=5)
    assert_almost_equal(mol.compute_density(points, "a", range(1, 6)), 0.5 * dens, decimal=5)


def test_molecule_grid_g09_fchk_uhf_ch4():
    # make an instance of molecule from fchk file
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_against_gaussian_ch4(molecule)


def test_molecule_grid_g09_fchk_rhf_ch4():
    # make an instance of molecule from fchk file
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_against_gaussian_ch4(molecule)


def test_molecule_grid_g09_wfn_uhf_ch4():
    # make an instance of molecule from wfn file
    with path("chemtools.data", "ch4_uhf_ccpvdz.wfn") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_against_gaussian_ch4(molecule)


def check_molecule_against_fortran_ch4(mol):
    # get expected data computed by Fortran code
    with path("chemtools.data", "data_fortran_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
        points, exp8, exp9 = data["points"], data["orb_08"], data["orb_09"]
        dens, grad, ke = data["dens"], data["grad"], data["ked_pd"]
    # check density & gradient
    assert_almost_equal(mol.compute_density(points, "ab", None), dens, decimal=6)
    assert_almost_equal(mol.compute_gradient(points, "ab", None), grad, decimal=6)
    assert_almost_equal(mol.compute_density(points, "a", None), 0.5 * dens, decimal=6)
    assert_almost_equal(mol.compute_density(points, "b", None), 0.5 * dens, decimal=6)
    # check density computed by summing squared mo expressions
    assert_almost_equal(mol.compute_density(points, "ab", range(1, 6)), dens, decimal=6)
    assert_almost_equal(mol.compute_density(points, "a", range(1, 6)), 0.5 * dens, decimal=6)
    assert_almost_equal(mol.compute_density(points, "b", range(1, 6)), 0.5 * dens, decimal=6)
    # check mo expression
    assert_almost_equal(mol.compute_molecular_orbital(points, "a", 8)[:, 0], exp8, decimal=6)
    assert_almost_equal(mol.compute_molecular_orbital(points, "b", 8)[:, 0], exp8, decimal=6)
    assert_almost_equal(mol.compute_molecular_orbital(points, "a", 9)[:, 0], exp9, decimal=6)
    assert_almost_equal(mol.compute_molecular_orbital(points, "b", 9)[:, 0], exp9, decimal=6)
    # check positive definite ke
    assert_almost_equal(mol.compute_ked(points, "ab"), ke, decimal=6)


def test_molecule_grid_fortran_fchk_uhf_ch4():
    # make an instance of molecule
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_against_fortran_ch4(molecule)


def test_molecule_grid_fortran_fchk_rhf_ch4():
    # make an instance of molecule
    with path("chemtools.data", "ch4_rhf_ccpvdz.fchk") as fname:
        molecule = Molecule.from_file(fname)
    check_molecule_against_fortran_ch4(molecule)


def test_molecule_basic_fchk_uhf_o2():
    with path("chemtools.data", "o2_uhf_virtual.fchk") as fname:
        mol = Molecule.from_file(fname)
    # print mol.nelectrons
    # check basic numerics
    assert_equal(mol.natom, 2)
    assert_equal(mol.nelectrons, (9, 7))
    assert_equal(mol.ao.nbasis, 44)
    assert_equal(mol.numbers, [8, 8])
    assert_equal(mol.pseudo_numbers, [8, 8])
    assert_equal(mol.mo.homo_index, (9, 7))
    assert_equal(mol.mo.lumo_index, (10, 8))
    assert_equal(mol.mo.occupation[0], np.array([1] * 9 + [0] * 35))
    assert_equal(mol.mo.occupation[1], np.array([1] * 7 + [0] * 37))
    assert_almost_equal(mol.energy, -1.496641407696776E+02, decimal=8)
    # check coordinates
    coord = np.array([[0.0, 0.0, 1.09452498], [0.0, 0.0, -1.09452498]])
    assert_almost_equal(mol.coordinates, coord, decimal=6)
    # energy of alpha orbitals
    orb_energy_a = np.array([-2.07520003E+01, -2.07511624E+01, -1.77073830E+00, -1.19176563E+00,
                             -8.67505123E-01, -8.67505123E-01, -7.86590608E-01, -5.41367609E-01,
                             -5.41367609E-01,  1.75551613E-01,  1.79578059E-01,  2.12136240E-01,
                              2.12136240E-01,  2.82746427E-01,  2.82746427E-01,  2.85693824E-01,
                              4.25803814E-01,  5.48137441E-01,  1.13330121E+00,  1.13563801E+00,
                              1.13563801E+00,  1.21130153E+00,  1.21130153E+00,  1.21869395E+00,
                              1.42273629E+00,  1.74301463E+00,  2.54671907E+00,  2.54671907E+00,
                              2.83314193E+00,  2.83314193E+00,  3.16996201E+00,  3.16996201E+00,
                              3.35513494E+00,  3.91418793E+00,  3.91418793E+00,  4.32509152E+00,
                              5.22427895E+00,  5.22427895E+00,  5.43561012E+00,  5.51122695E+00,
                              5.51122695E+00,  6.78081131E+00,  5.12932061E+01,  5.15031931E+01])
    # energy of beta orbitals
    orb_energy_b = np.array([-2.06970266E+01, -2.06955835E+01, -1.64286825E+00, -9.81871414E-01,
                             -7.18265821E-01, -6.01903968E-01, -6.01903968E-01,  1.04744787E-01,
                              1.04744787E-01,  1.82240025E-01,  1.84775146E-01,  2.25175431E-01,
                              2.25175431E-01,  2.81988319E-01,  3.22590360E-01,  3.22590360E-01,
                              4.31522630E-01,  6.25012892E-01,  1.14414927E+00,  1.21805822E+00,
                              1.21805822E+00,  1.24392742E+00,  1.30587863E+00,  1.30587863E+00,
                              1.45650216E+00,  1.79253113E+00,  2.62926234E+00,  2.62926234E+00,
                              2.95890359E+00,  2.95890359E+00,  3.32630000E+00,  3.32630000E+00,
                              3.42846106E+00,  3.99627997E+00,  3.99627997E+00,  4.36808390E+00,
                              5.33007026E+00,  5.33007026E+00,  5.45827702E+00,  5.61601037E+00,
                              5.61601037E+00,  6.81582045E+00,  5.13257489E+01,  5.15352582E+01])
    # check orbital energy
    assert_almost_equal(mol.mo.energy[0], orb_energy_a, decimal=6)
    assert_almost_equal(mol.mo.energy[1], orb_energy_b, decimal=6)
    assert_almost_equal(mol.mo.homo_energy[0], orb_energy_a[8], decimal=6)
    assert_almost_equal(mol.mo.homo_energy[1], orb_energy_b[6], decimal=6)
    assert_almost_equal(mol.mo.lumo_energy[0], orb_energy_a[9], decimal=6)
    assert_almost_equal(mol.mo.lumo_energy[1], orb_energy_b[7], decimal=6)
    assert_almost_equal(mol.mulliken_charges, 0.0, decimal=6)
    # check orbital coefficients
    assert_almost_equal(mol.mo.coefficient[0][:3, 0],
                        np.array([0.389497609, 0.333421243, 0.]), decimal=6)


def test_molecule_density_matrix_index_fchk_uhf_ch4():
    # check get_density_matrix for different values of the index/
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        mol = Molecule.from_file(fname)
    dm_full = mol.mo.compute_dm("a")._array
    # errors
    assert_raises(ValueError, mol.compute_dm, "a", [[1]])
    assert_raises(ValueError, mol.compute_dm, "a", [0])
    # one index
    for i in range(1, mol.ao.nbasis + 1):
        assert np.allclose(dm_full[i - 1, i - 1], mol.mo.compute_dm("a", i)._array)
    # multiple indices
    for i in range(1, mol.ao.nbasis + 1):
        for j in range(1, mol.ao.nbasis + 1):
            # NOTE: indices can be repeated
            indices = np.array([i -1, j - 1])
            assert np.allclose(dm_full[indices[:, None], indices[None, :]],
                               mol.mo.compute_dm("a", [i, j])._array)


def test_molecule_horton_h2o():
    with path("chemtools.data", "data_horton_fchk_h2o_ub3lyp_ccpvtz.npz") as fname:
        data = np.load(str(fname))
    with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as fname:
        mol = Molecule.from_file(fname)
    # check properties computed at nucleus against HORTON
    points = data["coords"]
    assert np.allclose(mol.compute_density(points), data["nuc_dens"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_gradient(points), data["nuc_grad"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_hessian(points), data["nuc_hess"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_ked(points), data["nuc_ked_pd"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_esp(points), data["nuc_esp"], rtol=0., atol=1.e-6)
    # check properties computed on a grid against HORTON
    assert np.allclose(mol.compute_density(data["points"]), data["dens"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_gradient(data["points"]), data["grad"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_hessian(data["points"]), data["hess"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_ked(data["points"]), data["ked_pd"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_esp(data["points"]), data["esp"], rtol=0., atol=1.e-6)


def test_molecule_horton_ch4():
    with path("chemtools.data", "data_horton_fchk_ch4_uhf_ccpvdz.npz") as fname:
        data = np.load(str(fname))
    with path("chemtools.data", "ch4_uhf_ccpvdz.fchk") as fname:
        mol = Molecule.from_file(fname)
    # check properties computed at nucleus against HORTON
    points = data["coords"]
    assert np.allclose(mol.compute_density(points), data["nuc_dens"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_gradient(points), data["nuc_grad"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_hessian(points), data["nuc_hess"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_ked(points), data["nuc_ked_pd"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_esp(points), data["nuc_esp"], rtol=0., atol=1.e-6)
    # check properties computed on a grid against HORTON
    assert np.allclose(mol.compute_density(data["points"]), data["dens"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_gradient(data["points"]), data["grad"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_hessian(data["points"]), data["hess"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_ked(data["points"]), data["ked_pd"], rtol=0., atol=1.e-6)
    assert np.allclose(mol.compute_esp(data["points"]), data["esp"], rtol=0., atol=1.e-6)
