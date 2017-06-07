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
"""Test chemtools.utils.molecule."""


import numpy as np
from numpy.testing import assert_raises
from test_data import get_cube_data_ch4_uhf_ccpvdz
from chemtools import context, HortonMolecule


def test_horton_molecule_check_raises():
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
    points = np.array([[0., 0., 0.]])
    # check invalid orbital spin argument
    assert_raises(ValueError, mol.compute_density_matrix_array, 'alphabeta')
    assert_raises(ValueError, mol.compute_density, points, spin='alph')
    assert_raises(ValueError, mol.compute_gradient, points, spin='ba')
    assert_raises(ValueError, mol.compute_hessian, points, spin='betaalpha')
    assert_raises(ValueError, mol.compute_kinetic_energy_density, points, spin='balpha')
    # check invalid points argument
    assert_raises(ValueError, mol.compute_density, [0., 0., 0.])
    assert_raises(ValueError, mol.compute_density, np.array([0., 0., 0.]))
    assert_raises(ValueError, mol.compute_density, np.array([[0., 0.]]))
    assert_raises(ValueError, mol.compute_gradient, [.1, .2, .3])
    assert_raises(ValueError, mol.compute_gradient, np.array([.1, .2, .3]))
    assert_raises(ValueError, mol.compute_gradient, np.array([[.1, 2., .3, .4]]))
    assert_raises(ValueError, mol.compute_hessian, [.5, .5, .5])
    assert_raises(ValueError, mol.compute_hessian, np.array([.5, .5, .5]))
    assert_raises(ValueError, mol.compute_hessian, np.array([[.5, 5.]]))
    assert_raises(ValueError, mol.compute_esp, [1., .5])
    assert_raises(ValueError, mol.compute_esp, np.array([1., .5, .25]))
    assert_raises(ValueError, mol.compute_esp, np.array([[1., .25]]))
    assert_raises(ValueError, mol.compute_kinetic_energy_density, [.5, 0., .2])
    assert_raises(ValueError, mol.compute_kinetic_energy_density, np.array([.5, 0., .2]))
    assert_raises(ValueError, mol.compute_kinetic_energy_density, np.array([[.5, 0., .2, .1, .3]]))
    # check invalid charges argument
    assert_raises(ValueError, mol.compute_esp, points, charges=np.array([6., 1., 1.]))
    assert_raises(ValueError, mol.compute_esp, points, charges=[6., 1., 1., 1., 1.])
    # check not implemented erros
    assert_raises(NotImplementedError, mol.compute_gradient, points, 'ab', [0, 1], None)
    assert_raises(NotImplementedError, mol.compute_gradient, points, 'a', [5, 7], None)
    assert_raises(NotImplementedError, mol.compute_gradient, points, 'b', [20], None)
    assert_raises(NotImplementedError, mol.compute_hessian, points, 'ab', [32, 33], None)
    assert_raises(NotImplementedError, mol.compute_esp, points, 'a', [33], None)
    assert_raises(NotImplementedError, mol.compute_kinetic_energy_density, points, 'a', [30])


def test_horton_molecule_basics_fchk_ch4_uhf():
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
    # check basic numerics
    np.testing.assert_equal(mol.natom, 5)
    np.testing.assert_equal(mol.nelectrons, (5, 5))
    np.testing.assert_equal(mol.nbasis, 34)
    np.testing.assert_equal(mol.numbers, [6, 1, 1, 1, 1])
    np.testing.assert_equal(mol.pseudo_numbers, [6, 1, 1, 1, 1])
    np.testing.assert_equal(mol.homo_index, (5, 5))
    np.testing.assert_equal(mol.lumo_index, (6, 6))
    np.testing.assert_equal(mol.orbital_occupation[0], np.array([1] * 5 + [0] * 29))
    np.testing.assert_equal(mol.orbital_occupation[1], np.array([1] * 5 + [0] * 29))
    np.testing.assert_almost_equal(mol.energy, -4.019868797400735E+01, decimal=8)
    # check coordinates
    coord = np.array([[-3.77945227E-05,  3.77945227E-05, -1.88972613E-05],
                      [ 1.04290206E+00,  1.50497789E+00,  9.34507367E-01],
                      [ 1.28607202E+00, -1.53098052E+00, -4.77307027E-01],
                      [-1.46467003E+00, -7.02997019E-01,  1.25954026E+00],
                      [-8.64474117E-01,  7.29131931E-01, -1.71670281E+00]])
    np.testing.assert_almost_equal(mol.coordinates, coord, decimal=6)
    # check charges
    esp = np.array([-0.502277518, 0.125567970, 0.125569655, 0.125566743, 0.125573150])
    np.testing.assert_almost_equal(mol.esp_charges, esp, decimal=6)
    npa = np.array([-0.791299253, 0.197824989, 0.197825250, 0.197824326, 0.197824689])
    np.testing.assert_almost_equal(mol.npa_charges, npa, decimal=6)
    mul = np.array([-0.139702704, 0.0349253868, 0.0349266071, 0.0349235395, 0.0349271707])
    np.testing.assert_almost_equal(mol.mulliken_charges, mul, decimal=6)


def test_horton_molecule_orbitals_fchk_ch4_uhf():
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
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
    np.testing.assert_almost_equal(mol.orbital_energy[0], orb_energy, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_energy[1], orb_energy, decimal=6)
    np.testing.assert_almost_equal(mol.homo_energy[0], orb_energy[4], decimal=6)
    np.testing.assert_almost_equal(mol.homo_energy[1], orb_energy[4], decimal=6)
    np.testing.assert_almost_equal(mol.lumo_energy[0], orb_energy[5], decimal=6)
    np.testing.assert_almost_equal(mol.lumo_energy[1], orb_energy[5], decimal=6)
    # check orbital coefficients
    orb_coeff = np.array([9.97287609E-01, 1.86004593E-02, -8.24772487E-03])
    np.testing.assert_almost_equal(mol.orbital_coefficient[1][:3, 0], orb_coeff, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_coefficient[0][:3, 0], orb_coeff, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_coefficient[1][0, 1], -0.188285003, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_coefficient[0][0, 1], -0.188285003, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_coefficient[1][-1, -1], 1.02960200, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_coefficient[0][-1, -1], 1.02960200, decimal=6)
    # check overlap matrix
    overlap = mol.compute_overlap()
    np.testing.assert_equal(overlap.shape, (34, 34))
    np.testing.assert_almost_equal(np.diag(overlap), np.ones(34), decimal=6)
    np.testing.assert_almost_equal(overlap, overlap.T, decimal=6)


def test_horton_molecule_density_matrix_fchk_ch4_uhf():
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
    # check density matrix against Gaussian (printed in log file)
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
    dm_array_a = mol.compute_density_matrix_array(spin='a')
    np.testing.assert_almost_equal(np.diag(dm_array_a), expected_diag, decimal=5)
    np.testing.assert_almost_equal(dm_array_a, dm_array_a.T, decimal=5)
    np.testing.assert_almost_equal(dm_array_a[0, 1:3], np.array([-0.04982, -0.05262]), decimal=5)
    np.testing.assert_almost_equal(dm_array_a[4, 14:21], expected_05, decimal=5)
    np.testing.assert_almost_equal(dm_array_a[14, 15:], expected_15, decimal=5)
    np.testing.assert_almost_equal(dm_array_a[18, 19:], expected_19, decimal=5)
    np.testing.assert_almost_equal(dm_array_a[28, 29:], expected_29, decimal=5)
    # check beta density matrix
    dm_array_b = mol.compute_density_matrix_array(spin='b')
    np.testing.assert_almost_equal(np.diag(dm_array_b), expected_diag, decimal=5)
    np.testing.assert_almost_equal(dm_array_b, dm_array_b.T, decimal=5)
    np.testing.assert_almost_equal(dm_array_b[0, 1:3], np.array([-0.04982, -0.05262]), decimal=5)
    np.testing.assert_almost_equal(dm_array_b[4, 14:21], expected_05, decimal=5)
    np.testing.assert_almost_equal(dm_array_b[14, 15:], expected_15, decimal=5)
    np.testing.assert_almost_equal(dm_array_b[18, 19:], expected_19, decimal=5)
    np.testing.assert_almost_equal(dm_array_b[28, 29:], expected_29, decimal=5)
    # check total density matrix
    dm_array_ab = mol.compute_density_matrix_array(spin='ab')
    np.testing.assert_almost_equal(np.diag(dm_array_ab), 2 * expected_diag, decimal=5)
    np.testing.assert_almost_equal(dm_array_ab, dm_array_ab.T, decimal=5)
    np.testing.assert_almost_equal(dm_array_ab[0, 1:3], 2*np.array([-0.04982, -0.05262]), decimal=5)
    np.testing.assert_almost_equal(dm_array_ab[4, 14:21], 2*expected_05, decimal=5)
    np.testing.assert_almost_equal(dm_array_ab[14, 15:], 2*expected_15, decimal=5)
    np.testing.assert_almost_equal(dm_array_ab[18, 19:], 2*expected_19, decimal=5)
    np.testing.assert_almost_equal(dm_array_ab[28, 29:], 2*expected_29, decimal=5)


def test_horton_molecule_grid_esp_fchk_ch4_uhf():
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
    # check esp against Gaussian (printed in log file)
    # check esp at the position of each nuclei (1.e-14 is added to avoid division by zero)
    # excluding the nucleus itself.
    point = mol.coordinates[0].reshape(1, 3) + 1.e-14
    charge = np.array([0., 1., 1., 1., 1.])
    np.testing.assert_almost_equal(mol.compute_esp(point, charges=charge), [-14.745629], decimal=5)
    point = mol.coordinates[1].reshape(1, 3) + 1.e-14
    charge = np.array([6., 0., 1., 1., 1.])
    np.testing.assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
    point = mol.coordinates[2].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 0., 1., 1.])
    np.testing.assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
    point = mol.coordinates[3].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 1., 0., 1.])
    np.testing.assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116067], decimal=5)
    point = mol.coordinates[4].reshape(1, 3) + 1.e-14
    charge = np.array([6., 1., 1., 1., 0.])
    np.testing.assert_almost_equal(mol.compute_esp(point, charges=charge), [-1.116065], decimal=5)
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
    np.testing.assert_almost_equal(mol.compute_esp(points), expected_esp, decimal=5)


def test_horton_molecule_grid_dens_gradient_fchk_ch4_uhf():
    # make an instance of molecule
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.fchk'))
    # check properties against Gaussian
    points, density, gradient, esp = get_cube_data_ch4_uhf_ccpvdz()
    np.testing.assert_almost_equal(mol.compute_density(points, 'ab'), density, decimal=5)
    np.testing.assert_almost_equal(mol.compute_gradient(points, 'ab'), gradient, decimal=5)
    np.testing.assert_almost_equal(mol.compute_esp(points, 'ab'), esp, decimal=5)


def test_horton_molecule_grid_dens_gradient_wfn_ch4_uhf():
    # make an instance of molecule
    mol = HortonMolecule.from_file(context.get_fn('test/ch4_uhf_ccpvdz.wfn'))
    # check properties against Gaussian
    points, density, gradient, esp = get_cube_data_ch4_uhf_ccpvdz()
    np.testing.assert_almost_equal(mol.compute_density(points, 'ab'), density, decimal=5)
    np.testing.assert_almost_equal(mol.compute_gradient(points, 'ab'), gradient, decimal=5)
    np.testing.assert_almost_equal(mol.compute_esp(points, 'ab'), esp, decimal=5)


def test_horton_molecule_basic_fchk_o2_uhf():
    mol = HortonMolecule.from_file(context.get_fn('test/o2_uhf_virtual.fchk'))
    print mol.nelectrons
    # check basic numerics
    np.testing.assert_equal(mol.natom, 2)
    np.testing.assert_equal(mol.nelectrons, (9, 7))
    np.testing.assert_equal(mol.nbasis, 44)
    np.testing.assert_equal(mol.numbers, [8, 8])
    np.testing.assert_equal(mol.pseudo_numbers, [8, 8])
    np.testing.assert_equal(mol.homo_index, (9, 7))
    np.testing.assert_equal(mol.lumo_index, (10, 8))
    np.testing.assert_equal(mol.orbital_occupation[0], np.array([1] * 9 + [0] * 35))
    np.testing.assert_equal(mol.orbital_occupation[1], np.array([1] * 7 + [0] * 37))
    np.testing.assert_almost_equal(mol.energy, -1.496641407696776E+02, decimal=8)
    # check coordinates
    coord = np.array([[0.0, 0.0, 1.09452498], [0.0, 0.0, -1.09452498]])
    np.testing.assert_almost_equal(mol.coordinates, coord, decimal=6)
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
    np.testing.assert_almost_equal(mol.orbital_energy[0], orb_energy_a, decimal=6)
    np.testing.assert_almost_equal(mol.orbital_energy[1], orb_energy_b, decimal=6)
    np.testing.assert_almost_equal(mol.homo_energy[0], orb_energy_a[8], decimal=6)
    np.testing.assert_almost_equal(mol.homo_energy[1], orb_energy_b[6], decimal=6)
    np.testing.assert_almost_equal(mol.lumo_energy[0], orb_energy_a[9], decimal=6)
    np.testing.assert_almost_equal(mol.lumo_energy[1], orb_energy_b[7], decimal=6)
    np.testing.assert_almost_equal(mol.mulliken_charges, 0.0, decimal=6)
    # check orbital coefficients
    np.testing.assert_almost_equal(mol.orbital_coefficient[0][:3, 0],
                                   np.array([0.389497609, 0.333421243, 0.]), decimal=6)
