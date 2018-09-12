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
"""Test chemtools.analysis.orbitalbased."""


import numpy as np

from numpy.testing import assert_raises, assert_array_almost_equal, assert_allclose

from horton import IOData
from chemtools.wrappers.molecule import Molecule
from chemtools import context, CubeGen, OrbitalAnalysis


def test_toolbox_orbitalbased_raises():
    # check file name
    assert_raises(ValueError, OrbitalAnalysis.from_file, "gibberish", np.array([0.0, 1.0]))
    filename = context.get_fn('test/h2o_dimer_pbe_sto3g.wf')
    assert_raises(ValueError, OrbitalAnalysis.from_file, filename, np.array([0.0, 1.0]))


def test_orbital_analysis_from_file_ch4_uhf_ccpvdz():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    mol = IOData.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalAnalysis.from_file(file_path, cube.points)

    # density results obtained from Fortran code:
    result = [0.00003304, 0.00053319, 0.00019292, 0.00111552, 0.00679461,
              0.00153604, 0.00015922, 0.00030448, 0.00003973, 0.00045413,
              0.00754940, 0.00043585, 0.01189345, 120.661406, 0.00488532,
              0.00085596, 0.00715178, 0.00084528, 0.00015549, 0.00192313,
              0.00004713, 0.00034775, 0.00541748, 0.00042815, 0.00003358,
              0.00103735, 0.00021200]
    # check density array
    test = orbtool.density
    assert_array_almost_equal(test, result, decimal=6)

    # gradient results obtained from Fortran code:
    result = [[ 0.00004568,  0.00005560,  0.00004170], [ 0.00071421,  0.00090481,  0.00031958],
              [ 0.00022178,  0.00030427, -0.00024400], [ 0.00214134,  0.00057641,  0.00147697],
              [ 0.01414091, -0.00274820,  0.00462187], [ 0.00242062, -0.00084759, -0.00266456],
              [ 0.00022914, -0.00024055,  0.00015719], [ 0.00043634, -0.00045833, -0.00009951],
              [ 0.00005359, -0.00006511, -0.00004867], [ 0.00020864,  0.00060902,  0.00073946],
              [ 0.00682003,  0.01554390, -0.00220110], [-0.00020073,  0.00066217, -0.00062608],
              [-0.00943599,  0.00784003,  0.02386404], [-8.40090166,  8.40090054, -4.20041995],
              [-0.00254422, -0.00035925, -0.01023804], [-0.00045752, -0.00160984,  0.00110110],
              [ 0.00513284, -0.01477175,  0.00413479], [ 0.00054025, -0.00119568, -0.00147869],
              [-0.00018878,  0.00016871,  0.00025415], [-0.00340647,  0.00301688, -0.00076147],
              [-0.00005689,  0.00005325, -0.00008253], [-0.00048891, -0.00014071,  0.00053314],
              [-0.01141306, -0.00213317, -0.00009663], [-0.00061368,  0.00022324, -0.00065764],
              [-0.00005222, -0.00004309,  0.00005006], [-0.00184030, -0.00150716,  0.00067124],
              [-0.00029870, -0.00024190, -0.00031205]]
    # check gradient
    test = orbtool.gradient
    assert_array_almost_equal(test, result, decimal=6)

    # Weizsacker KE results obtained from Fortran code:
    result = [0.00002617, 0.00033546, 0.00013043, 0.00079549, 0.00421069,
              0.00111306, 0.00010605, 0.00016847, 0.00002982, 0.00026458,
              0.00485089, 0.00024972, 0.00756715, 0.16450352, 0.00285088,
              0.00058608, 0.00457311, 0.00057792, 0.00010346, 0.00138352,
              0.00003417, 0.00019521, 0.00311071, 0.00025077, 0.00002639,
              0.00073611, 0.00014453]
    # check Weizsacker kinetic energy density
    test = orbtool.weizsacker_kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # Thomas-Fermi KE results obtained from Fortran code:
    result = [0.00000010, 0.00001007, 0.00000185,    0.00003445, 0.00069986,
              0.00005871, 0.00000134, 0.00000396,    0.00000013, 0.00000770,
              0.00083417, 0.00000719, 0.00177930, 8459.58828066, 0.00040385,
              0.00002216, 0.00076224, 0.00002170,    0.00000129, 0.00008539,
              0.00000018, 0.00000494, 0.00047980,    0.00000698, 0.00000010,
              0.00003052, 0.00000216]
    # check Thomas-Fermi kinetic energy density
    test = orbtool.thomas_fermi_kinetic_energy_density
    assert_allclose(test, result, rtol=1e-08, atol=1e-08)

    # Positive Definite KE results obtained from Fortran code:
    result = [0.00002941, 0.00036577, 0.00013191, 0.00082047, 0.00530631,
              0.00113415, 0.00010757, 0.00020285, 0.00003336, 0.00030256,
              0.00574122, 0.00028868, 0.00824769, 5.43470706, 0.00392851,
              0.00061677, 0.00549151, 0.00060841, 0.00010513, 0.00140232,
              0.00003700, 0.00023018, 0.00437588, 0.00028569, 0.00002990,
              0.00076019, 0.00014593]
    # check positive definite kinetic energy density
    test = orbtool.kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # HOMO orbital expansion results obtained from Fortran code:
    homo_result = [ 0.00375367,  0.03178078,  0.05348131,  0.02726712,  0.01400675,
                    0.07445751,  0.02320186,  0.00425437,  0.00577119, -0.02781335,
                   -0.04082439,  0.03324904,  0.03689038,  0.00002702,  0.02901776,
                    0.03288809, -0.04124146, -0.05165096, -0.02892164, -0.05043417,
                   -0.00997580, -0.00024128,  0.00307327, -0.02688246,  0.00031094,
                   -0.04755402, -0.04872100]
    # LUMO orbital expansion results obtained from Fortran code:
    lumo_result = [ 0.01019334, -0.03985858, -0.04740753,  0.08131806, -0.02842964,
                   -0.06900307,  0.06019631,  0.01863115, -0.01397425,  0.00146082,
                   -0.00519832, -0.01734235,  0.04867160,  0.00004154,  0.01719274,
                    0.06485919, -0.01288644, -0.01262115, -0.00838637, -0.00849781,
                    0.00009801, -0.00910972, -0.01612369,  0.00810636,  0.00485805,
                   -0.02029558, -0.00936400]
    # check homo expansion
    test = orbtool.orbitals_exp(8)
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    # check lumo expansion
    test = orbtool.orbitals_exp(9)
    assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    # check homo & lumo expansion with list of orbital indices
    test = orbtool.orbitals_exp([8, 9])
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with tuple of orbital indices
    test = orbtool.orbitals_exp((8, 9))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with array of orbital indices
    test = orbtool.orbitals_exp(np.array([8, 9]))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)


def test_orbital_analysis_from_file_ch4_rhf_ccpvdz():
    file_path = context.get_fn('test/ch4_rhf_ccpvdz.fchk')
    mol = IOData.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalAnalysis.from_file(file_path, cube.points)

    # density results obtained from Fortran code:
    result = [0.00003304, 0.00053319, 0.00019292, 0.00111552, 0.00679461,
              0.00153604, 0.00015922, 0.00030448, 0.00003973, 0.00045413,
              0.00754940, 0.00043585, 0.01189345, 120.661406, 0.00488532,
              0.00085596, 0.00715178, 0.00084528, 0.00015549, 0.00192313,
              0.00004713, 0.00034775, 0.00541748, 0.00042815, 0.00003358,
              0.00103735, 0.00021200]
    # check density array
    test = orbtool.density
    assert_array_almost_equal(test, result, decimal=6)

    # gradient results obtained from Fortran code:
    result = [[ 0.00004568,  0.00005560,  0.00004170], [ 0.00071421,  0.00090481,  0.00031958],
              [ 0.00022178,  0.00030427, -0.00024400], [ 0.00214134,  0.00057641,  0.00147697],
              [ 0.01414091, -0.00274820,  0.00462187], [ 0.00242062, -0.00084759, -0.00266456],
              [ 0.00022914, -0.00024055,  0.00015719], [ 0.00043634, -0.00045833, -0.00009951],
              [ 0.00005359, -0.00006511, -0.00004867], [ 0.00020864,  0.00060902,  0.00073946],
              [ 0.00682003,  0.01554390, -0.00220110], [-0.00020073,  0.00066217, -0.00062608],
              [-0.00943599,  0.00784003,  0.02386404], [-8.40090166,  8.40090054, -4.20041995],
              [-0.00254422, -0.00035925, -0.01023804], [-0.00045752, -0.00160984,  0.00110110],
              [ 0.00513284, -0.01477175,  0.00413479], [ 0.00054025, -0.00119568, -0.00147869],
              [-0.00018878,  0.00016871,  0.00025415], [-0.00340647,  0.00301688, -0.00076147],
              [-0.00005689,  0.00005325, -0.00008253], [-0.00048891, -0.00014071,  0.00053314],
              [-0.01141306, -0.00213317, -0.00009663], [-0.00061368,  0.00022324, -0.00065764],
              [-0.00005222, -0.00004309,  0.00005006], [-0.00184030, -0.00150716,  0.00067124],
              [-0.00029870, -0.00024190, -0.00031205]]
    # check gradient
    test = orbtool.gradient
    assert_array_almost_equal(test, result, decimal=6)

    # Weizsacker KE results obtained from Fortran code:
    result = [0.00002617, 0.00033546, 0.00013043, 0.00079549, 0.00421069,
              0.00111306, 0.00010605, 0.00016847, 0.00002982, 0.00026458,
              0.00485089, 0.00024972, 0.00756715, 0.16450352, 0.00285088,
              0.00058608, 0.00457311, 0.00057792, 0.00010346, 0.00138352,
              0.00003417, 0.00019521, 0.00311071, 0.00025077, 0.00002639,
              0.00073611, 0.00014453]
    # check Weizsacker kinetic energy density
    test = orbtool.weizsacker_kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # Thomas-Fermi KE results obtained from Fortran code:
    result = [0.00000010, 0.00001007, 0.00000185,    0.00003445, 0.00069986,
              0.00005871, 0.00000134, 0.00000396,    0.00000013, 0.00000770,
              0.00083417, 0.00000719, 0.00177930, 8459.58828066, 0.00040385,
              0.00002216, 0.00076224, 0.00002170,    0.00000129, 0.00008539,
              0.00000018, 0.00000494, 0.00047980,    0.00000698, 0.00000010,
              0.00003052, 0.00000216]
    # check Thomas-Fermi kinetic energy density
    test = orbtool.thomas_fermi_kinetic_energy_density
    assert_allclose(test, result, rtol=1e-08, atol=1e-08)

    # Positive Definite KE results obtained from Fortran code:
    result = [0.00002941, 0.00036577, 0.00013191, 0.00082047, 0.00530631,
              0.00113415, 0.00010757, 0.00020285, 0.00003336, 0.00030256,
              0.00574122, 0.00028868, 0.00824769, 5.43470706, 0.00392851,
              0.00061677, 0.00549151, 0.00060841, 0.00010513, 0.00140232,
              0.00003700, 0.00023018, 0.00437588, 0.00028569, 0.00002990,
              0.00076019, 0.00014593]
    # check positive definite kinetic energy density
    test = orbtool.kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # HOMO orbital expansion results obtained from Fortran code:
    homo_result = [ 0.00375367,  0.03178078,  0.05348131,  0.02726712,  0.01400675,
                    0.07445751,  0.02320186,  0.00425437,  0.00577119, -0.02781335,
                   -0.04082439,  0.03324904,  0.03689038,  0.00002702,  0.02901776,
                    0.03288809, -0.04124146, -0.05165096, -0.02892164, -0.05043417,
                   -0.00997580, -0.00024128,  0.00307327, -0.02688246,  0.00031094,
                   -0.04755402, -0.04872100]
    # LUMO orbital expansion results obtained from Fortran code:
    lumo_result = [ 0.01019334, -0.03985858, -0.04740753,  0.08131806, -0.02842964,
                   -0.06900307,  0.06019631,  0.01863115, -0.01397425,  0.00146082,
                   -0.00519832, -0.01734235,  0.04867160,  0.00004154,  0.01719274,
                    0.06485919, -0.01288644, -0.01262115, -0.00838637, -0.00849781,
                    0.00009801, -0.00910972, -0.01612369,  0.00810636,  0.00485805,
                   -0.02029558, -0.00936400]
    # check homo expansion
    test = orbtool.orbitals_exp(8)
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    # check lumo expansion
    test = orbtool.orbitals_exp(9)
    assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    # check homo & lumo expansion with list of orbital indices
    test = orbtool.orbitals_exp([8, 9])
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with tuple of orbital indices
    test = orbtool.orbitals_exp((8, 9))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with array of orbital indices
    test = orbtool.orbitals_exp(np.array([8, 9]))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)


def test_orbital_analysis_from_molecule_ch4_uhf_ccpvdz():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    mol = IOData.from_file(file_path)
    rmol = IOData.from_file(file_path)
    del rmol.exp_beta

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalAnalysis.from_molecule(Molecule(mol), cube.points)
    rorbtool = OrbitalAnalysis.from_molecule(Molecule(rmol), cube.points)

    # density results obtained from Fortran code:
    result = [0.00003304, 0.00053319, 0.00019292, 0.00111552, 0.00679461,
              0.00153604, 0.00015922, 0.00030448, 0.00003973, 0.00045413,
              0.00754940, 0.00043585, 0.01189345, 120.661406, 0.00488532,
              0.00085596, 0.00715178, 0.00084528, 0.00015549, 0.00192313,
              0.00004713, 0.00034775, 0.00541748, 0.00042815, 0.00003358,
              0.00103735, 0.00021200]
    # check density array
    test = orbtool.density
    assert_array_almost_equal(test, result, decimal=6)
    test = rorbtool.density
    assert_array_almost_equal(test, result, decimal=6)

    # gradient results obtained from Fortran code:
    result = [[ 0.00004568,  0.00005560,  0.00004170], [ 0.00071421,  0.00090481,  0.00031958],
              [ 0.00022178,  0.00030427, -0.00024400], [ 0.00214134,  0.00057641,  0.00147697],
              [ 0.01414091, -0.00274820,  0.00462187], [ 0.00242062, -0.00084759, -0.00266456],
              [ 0.00022914, -0.00024055,  0.00015719], [ 0.00043634, -0.00045833, -0.00009951],
              [ 0.00005359, -0.00006511, -0.00004867], [ 0.00020864,  0.00060902,  0.00073946],
              [ 0.00682003,  0.01554390, -0.00220110], [-0.00020073,  0.00066217, -0.00062608],
              [-0.00943599,  0.00784003,  0.02386404], [-8.40090166,  8.40090054, -4.20041995],
              [-0.00254422, -0.00035925, -0.01023804], [-0.00045752, -0.00160984,  0.00110110],
              [ 0.00513284, -0.01477175,  0.00413479], [ 0.00054025, -0.00119568, -0.00147869],
              [-0.00018878,  0.00016871,  0.00025415], [-0.00340647,  0.00301688, -0.00076147],
              [-0.00005689,  0.00005325, -0.00008253], [-0.00048891, -0.00014071,  0.00053314],
              [-0.01141306, -0.00213317, -0.00009663], [-0.00061368,  0.00022324, -0.00065764],
              [-0.00005222, -0.00004309,  0.00005006], [-0.00184030, -0.00150716,  0.00067124],
              [-0.00029870, -0.00024190, -0.00031205]]
    # check gradient
    test = orbtool.gradient
    assert_array_almost_equal(test, result, decimal=6)
    test = rorbtool.gradient
    assert_array_almost_equal(test, result, decimal=6)

    # Weizsacker KE results obtained from Fortran code:
    result = [0.00002617, 0.00033546, 0.00013043, 0.00079549, 0.00421069,
              0.00111306, 0.00010605, 0.00016847, 0.00002982, 0.00026458,
              0.00485089, 0.00024972, 0.00756715, 0.16450352, 0.00285088,
              0.00058608, 0.00457311, 0.00057792, 0.00010346, 0.00138352,
              0.00003417, 0.00019521, 0.00311071, 0.00025077, 0.00002639,
              0.00073611, 0.00014453]
    # check Weizsacker kinetic energy density
    test = orbtool.weizsacker_kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)
    test = rorbtool.weizsacker_kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # Thomas-Fermi KE results obtained from Fortran code:
    result = [0.00000010, 0.00001007, 0.00000185,    0.00003445, 0.00069986,
              0.00005871, 0.00000134, 0.00000396,    0.00000013, 0.00000770,
              0.00083417, 0.00000719, 0.00177930, 8459.58828066, 0.00040385,
              0.00002216, 0.00076224, 0.00002170,    0.00000129, 0.00008539,
              0.00000018, 0.00000494, 0.00047980,    0.00000698, 0.00000010,
              0.00003052, 0.00000216]
    # check Thomas-Fermi kinetic energy density
    test = orbtool.thomas_fermi_kinetic_energy_density
    assert_allclose(test, result, rtol=1e-08, atol=1e-08)
    test = rorbtool.thomas_fermi_kinetic_energy_density
    assert_allclose(test, result, rtol=1e-08, atol=1e-08)

    # Positive Definite KE results obtained from Fortran code:
    result = [0.00002941, 0.00036577, 0.00013191, 0.00082047, 0.00530631,
              0.00113415, 0.00010757, 0.00020285, 0.00003336, 0.00030256,
              0.00574122, 0.00028868, 0.00824769, 5.43470706, 0.00392851,
              0.00061677, 0.00549151, 0.00060841, 0.00010513, 0.00140232,
              0.00003700, 0.00023018, 0.00437588, 0.00028569, 0.00002990,
              0.00076019, 0.00014593]
    # check positive definite kinetic energy density
    test = orbtool.kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)
    test = rorbtool.kinetic_energy_density
    assert_array_almost_equal(test, result, decimal=6)

    # HOMO orbital expansion results obtained from Fortran code:
    homo_result = [ 0.00375367,  0.03178078,  0.05348131,  0.02726712,  0.01400675,
                    0.07445751,  0.02320186,  0.00425437,  0.00577119, -0.02781335,
                   -0.04082439,  0.03324904,  0.03689038,  0.00002702,  0.02901776,
                    0.03288809, -0.04124146, -0.05165096, -0.02892164, -0.05043417,
                   -0.00997580, -0.00024128,  0.00307327, -0.02688246,  0.00031094,
                   -0.04755402, -0.04872100]
    # LUMO orbital expansion results obtained from Fortran code:
    lumo_result = [ 0.01019334, -0.03985858, -0.04740753,  0.08131806, -0.02842964,
                   -0.06900307,  0.06019631,  0.01863115, -0.01397425,  0.00146082,
                   -0.00519832, -0.01734235,  0.04867160,  0.00004154,  0.01719274,
                    0.06485919, -0.01288644, -0.01262115, -0.00838637, -0.00849781,
                    0.00009801, -0.00910972, -0.01612369,  0.00810636,  0.00485805,
                   -0.02029558, -0.00936400]
    # check homo expansion
    test = orbtool.orbitals_exp(8)
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    test = rorbtool.orbitals_exp(8)
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    # check lumo expansion
    test = orbtool.orbitals_exp(9)
    assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    test = rorbtool.orbitals_exp(9)
    assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    # check homo & lumo expansion with list of orbital indices
    test = orbtool.orbitals_exp([8, 9])
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    test = rorbtool.orbitals_exp([8, 9])
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with tuple of orbital indices
    test = orbtool.orbitals_exp((8, 9))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    test = rorbtool.orbitals_exp((8, 9))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with array of orbital indices
    test = orbtool.orbitals_exp(np.array([8, 9]))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    test = rorbtool.orbitals_exp(np.array([8, 9]))
    assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
