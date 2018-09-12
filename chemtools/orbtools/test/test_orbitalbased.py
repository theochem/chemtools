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
"""Test chemtools.toolbox.orbitalbased."""

from numpy.testing import assert_raises
import numpy as np
from horton import IOData
from chemtools.utils.wrappers import HortonMolecule
from chemtools import context
from chemtools.utils import CubeGen
from chemtools.orbtools.orbitalbased import OrbitalLocalTool


def test_orbital_tool_ch4_uhf_ccpvdz():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    mol = HortonMolecule.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalLocalTool(mol, cube.points)

    # density results obtained from Fortran code:
    result = [0.00003304, 0.00053319, 0.00019292, 0.00111552, 0.00679461,
              0.00153604, 0.00015922, 0.00030448, 0.00003973, 0.00045413,
              0.00754940, 0.00043585, 0.01189345, 120.661406, 0.00488532,
              0.00085596, 0.00715178, 0.00084528, 0.00015549, 0.00192313,
              0.00004713, 0.00034775, 0.00541748, 0.00042815, 0.00003358,
              0.00103735, 0.00021200]
    # check density array
    test = orbtool.density
    np.testing.assert_array_almost_equal(test, result, decimal=6)

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
    np.testing.assert_array_almost_equal(test, result, decimal=6)

    # Weizsacker KE results obtained from Fortran code:
    result = [0.00002617, 0.00033546, 0.00013043, 0.00079549, 0.00421069,
              0.00111306, 0.00010605, 0.00016847, 0.00002982, 0.00026458,
              0.00485089, 0.00024972, 0.00756715, 0.16450352, 0.00285088,
              0.00058608, 0.00457311, 0.00057792, 0.00010346, 0.00138352,
              0.00003417, 0.00019521, 0.00311071, 0.00025077, 0.00002639,
              0.00073611, 0.00014453]
    # check Weizsacker kinetic energy density
    test = orbtool.weizsacker_kinetic_energy_density
    np.testing.assert_array_almost_equal(test, result, decimal=6)

    # Thomas-Fermi KE results obtained from Fortran code:
    result = [0.00000010, 0.00001007, 0.00000185,    0.00003445, 0.00069986,
              0.00005871, 0.00000134, 0.00000396,    0.00000013, 0.00000770,
              0.00083417, 0.00000719, 0.00177930, 8459.58828066, 0.00040385,
              0.00002216, 0.00076224, 0.00002170,    0.00000129, 0.00008539,
              0.00000018, 0.00000494, 0.00047980,    0.00000698, 0.00000010,
              0.00003052, 0.00000216]
    # check Thomas-Fermi kinetic energy density
    test = orbtool.thomas_fermi_kinetic_energy_density
    np.testing.assert_allclose(test, result, rtol=1e-08, atol=1e-08)

    # Positive Definite KE results obtained from Fortran code:
    result = [0.00002941, 0.00036577, 0.00013191, 0.00082047, 0.00530631,
              0.00113415, 0.00010757, 0.00020285, 0.00003336, 0.00030256,
              0.00574122, 0.00028868, 0.00824769, 5.43470706, 0.00392851,
              0.00061677, 0.00549151, 0.00060841, 0.00010513, 0.00140232,
              0.00003700, 0.00023018, 0.00437588, 0.00028569, 0.00002990,
              0.00076019, 0.00014593]
    # check positive definite kinetic energy density
    test = orbtool.kinetic_energy_density
    np.testing.assert_array_almost_equal(test, result, decimal=6)

    # HOMO orbital expansion results obtained from Fortran code:
    homo_result = [-0.00249675, -0.01320839, -0.00722774, -0.01116065, -0.05149882,
                   -0.02020972, -0.00320152, -0.00673363, -0.00218699,  0.00155883,
                    0.00418715, -0.00495368, -0.01006866,  0.00003386, -0.00863501,
                   -0.00162681,  0.01993861,  0.00698478,  0.00422124,  0.01551245,
                    0.00234245,  0.00843082,  0.04407199,  0.01013363,  0.00238077,
                    0.01567561,  0.00627160]
    # LUMO orbital expansion results obtained from Fortran code:
    lumo_result = [-0.01533933, -0.03700800, -0.03240375, -0.04711904, -0.04017300,
                   -0.05128246, -0.02993671, -0.02987038, -0.01693414, -0.03617181,
                   -0.04009498, -0.03576572, -0.04418628, -0.88010359, -0.03575185,
                   -0.04420089, -0.03958405, -0.04396912, -0.02985073, -0.05428580,
                   -0.01779872, -0.03181989, -0.03822247, -0.03474619, -0.01563934,
                   -0.04577599, -0.03352035]
    # check homo expansion
    test = orbtool.orbitals_exp(5)
    np.testing.assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    # check homo expansion (beta should equal alpha)
    test = orbtool.orbitals_exp(np.array([5]), spin='beta')
    np.testing.assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    # check lumo expansion
    test = orbtool.orbitals_exp(6)
    np.testing.assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    # check lumo expansion (beta should equal alpha)
    test = orbtool.orbitals_exp(np.array([6]), spin='beta')
    np.testing.assert_array_almost_equal(test[:, 0], lumo_result, decimal=6)
    # check homo & lumo expansion with list of orbital indices
    test = orbtool.orbitals_exp([5, 6])
    np.testing.assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    np.testing.assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with tuple of orbital indices
    test = orbtool.orbitals_exp((5, 6))
    np.testing.assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    np.testing.assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)
    # check homo & lumo expansion with array of orbital indices
    test = orbtool.orbitals_exp(np.array([5, 6]))
    np.testing.assert_array_almost_equal(test[:, 0], homo_result, decimal=6)
    np.testing.assert_array_almost_equal(test[:, 1], lumo_result, decimal=6)

    # check spin chemical potential
    # spin chemical potential result obtained from manually ecaluating the formula:
    result = [-0.1606881912, -0.1606881912]
    test = orbtool.spin_chemical_potential(25000.0)
    np.testing.assert_array_almost_equal(test, result, decimal=6)

    # check temperature dependent density at 25000K
    # Density results obtained from Fortran code:
    result = [0.00003919, 0.00058389, 0.00025735, 0.00121919, 0.00679814, 0.00166896,   0.00021239,
              0.00032721, 0.00004793, 0.00049444, 0.00755612, 0.00047382, 0.01189356, 120.67872015,
              0.00489612, 0.00093781, 0.00715898, 0.00092623, 0.00020734, 0.00207881,   0.00005839,
              0.00037579, 0.00542211, 0.00046599, 0.00003973, 0.00113522, 0.00028249]

    # check density array
    test = orbtool.temperature_dependent_density(25000.0)
    np.testing.assert_array_almost_equal(test, result, decimal=6)
    # check density array only using spin_chemical_potential
    test = orbtool.temperature_dependent_density(25000.0,
                                                 spin_chemical_potential=[-0.160688, -0.160688])
    np.testing.assert_array_almost_equal(test, result, decimal=6)

    # check KeyError
    assert_raises(KeyError, orbtool.orbitals_exp, np.array([9]), spin='error')


def test_orbital_tool_h2o_b3lyp_sto3g():
    file_path = context.get_fn('test/water_b3lyp_sto3g.fchk')
    mol = HortonMolecule.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalLocalTool(mol, cube.points)

    # mep results obtained from Fortran code:
    expected = np.array([-0.01239766, -0.02982537, -0.02201149,   -0.01787292, -0.05682143,
                         -0.02503563, -0.00405942, -0.00818772,   -0.00502268,  0.00321181,
                         -0.03320573, -0.02788605,  0.02741914, 1290.21135500, -0.03319778,
                          0.01428660,  0.10127092,  0.01518299,    0.01530548,  0.00197975,
                         -0.00894206,  0.04330806,  0.03441681,   -0.00203017,  0.02272626,
                          0.03730846,  0.01463959])

    test = orbtool.mep(mol.coordinates, mol.pseudo_numbers)

    np.testing.assert_array_almost_equal(test, expected, decimal=6)

    # check spin chemical potential
    # spin chemical potential result obtained from manually ecaluating the formula
    result = [0.10821228040, 0.10821228040]
    test = orbtool.spin_chemical_potential(25000.0)
    np.testing.assert_almost_equal(test, result, decimal=6)


def test_orbital_tool_elf_h2o_dimer():
    file_path = context.get_fn('test/h2o_dimer_pbe_sto3g.fchk')
    # load fchk
    mol = HortonMolecule.from_file(file_path)
    # Check against elf created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    elf_cube_path = context.get_fn('test/h2o_dimer_pbe_sto3g-elf.cube')
    elf = IOData.from_file(elf_cube_path)
    result = elf.cube_data.flatten()
    # Build OrbitalLocal tool
    cube = CubeGen.from_cube(elf_cube_path)
    orbtool = OrbitalLocalTool(mol, cube.points)
    test = orbtool.elf

    np.testing.assert_equal(test.shape, result.shape)
    np.testing.assert_array_almost_equal(test, result, decimal=5)


def test_localip_ch4_uhf_ccpvdz_alpha():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    mol = HortonMolecule.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalLocalTool(mol, cube.points)

    # local ip obtained with a Mathematica notebook:
    expected = [-0.583314, -0.587023, -0.565793,
                -0.582555, -0.623935, -0.581330,
                -0.565933, -0.595541, -0.579717,
                -0.588923, -0.620614, -0.589438,
                -0.611751, -10.88100, -0.631027,
                -0.583850, -0.621805, -0.583896,
                -0.566106, -0.580425, -0.575995,
                -0.592862, -0.629875, -0.589605,
                -0.583253, -0.582762, -0.565703]

    # compute the local ionization potential
    test = orbtool.local_ip

    np.testing.assert_almost_equal(expected, test, decimal=4)


def test_localip_ch4_uhf_ccpvdz_both():
    file_path = context.get_fn('test/ch4_uhf_ccpvdz.fchk')
    # load fchk
    mol = HortonMolecule.from_file(file_path)

    # creating cube file:
    ori = np.array([-3.000000, -3.000000, -3.000000])
    ax = np.array([[ 3.000000,  0.000000,  0.000000],
                   [ 0.000000,  3.000000,  0.000000],
                   [ 0.000000,  0.000000,  3.000000]])
    sh = np.array([3, 3, 3])
    cube = CubeGen(mol.numbers, mol.pseudo_numbers, mol.coordinates, ori, ax, sh)

    # initialize OrbitalLocalTool:
    orbtool = OrbitalLocalTool(mol, cube.points)

    # local ip obtained with a Mathematica notebook:
    expected = [-0.583314, -0.587023, -0.565793,
                -0.582555, -0.623935, -0.581330,
                -0.565933, -0.595541, -0.579717,
                -0.588923, -0.620614, -0.589438,
                -0.611751, -10.88100, -0.631027,
                -0.583850, -0.621805, -0.583896,
                -0.566106, -0.580425, -0.575995,
                -0.592862, -0.629875, -0.589605,
                -0.583253, -0.582762, -0.565703]
    # compute the local ionization potential
    test = orbtool.local_ip

    np.testing.assert_almost_equal(expected, test, decimal=4)
