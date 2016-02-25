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
import numpy as np, tempfile, shutil
from horton import IOData
from chemtools import *
from contextlib import contextmanager

@contextmanager
def tmpdir(name):
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)

def test_analyze_ch4_fchk_linear():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # IP = -E(HOMO) & EA = E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), 1.93295185E-01, -4.019868797400735E+01
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path, model='linear')
    np.testing.assert_almost_equal(desp.globaltool.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(11.), energy - ea, decimal=8)
    # Check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.globaltool.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.electron_affinity, ea, decimal=8)
    # Check chemical-potential, chemical-hardness & hyper-hardness
    print desp.globaltool.mu, desp.globaltool.chemical_potential
    np.testing.assert_equal(desp.globaltool.mu, None)
    np.testing.assert_equal(desp.globaltool.chemical_potential, None)
    np.testing.assert_equal(desp.globaltool.eta, None)
    np.testing.assert_equal(desp.globaltool.chemical_hardness, None)
    np.testing.assert_equal(desp.globaltool.hyper_hardness(2), None)
    np.testing.assert_equal(desp.globaltool.hyper_hardness(3), None)
    np.testing.assert_equal(desp.globaltool.hyper_hardness(4), None)
    # Check mu+, mu-, mu0
    np.testing.assert_almost_equal(desp.globaltool.mu_plus, -ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.mu_minus, -ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.mu_zero, -0.5 * (ip + ea), decimal=8)
    # Derivatives of E(N)
    np.testing.assert_equal(desp.globaltool.energy_derivative(10, 3), None)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.5, 4), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.0, 3), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10.4, 2), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(11, 1), -ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.0, 1), -ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.7, 1), -ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10.5, 1), -ea, decimal=8)
    # Check softness & hyper-softness
    # Check N_max and related descriptors


def test_analyze_ch4_fchk_quadratic():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # IP = -E(HOMO) & EA = E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), 1.93295185E-01, -4.019868797400735E+01
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path, model='quadratic')
    np.testing.assert_almost_equal(desp.globaltool.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(11.), energy - ea, decimal=8)
    # Check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.globaltool.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.electron_affinity, ea, decimal=8)
    # Check chemical-potential, chemical-hardness & hyper-hardness
    mu, eta = -0.5 * (ip + ea), ip - ea
    np.testing.assert_almost_equal(desp.globaltool.mu, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_potential, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.eta, eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_hardness, eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(2), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(3), 0.0, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(4), 0.0, decimal=8)
    # Check softness & hyper-softness
    np.testing.assert_almost_equal(desp.globaltool.softness, 1.0/eta, decimal=8)
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(2), 0.0, decimal=8)
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(3), 0.0, decimal=8)
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(4), 0.0, decimal=8)
    # Check N_max and related descriptors
    np.testing.assert_almost_equal(desp.globaltool.n0, 10, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.n_max, 10 - mu/eta, decimal=8)
    # Check Electrophilicity
    value = 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.globaltool.electrophilicity, value, decimal=8)
    value = (ip + ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.globaltool.electrophilicity, value, decimal=8)
    # Check Nucleofugality
    value = (ip - 3 * ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.globaltool.nucleofugality, value, decimal=8)
    value = (mu + eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.globaltool.nucleofugality, value, decimal=8)
    value = - ea + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.globaltool.nucleofugality, value, decimal=8)
    # Check Electrofugality
    value = (3 * ip - ea)**2 / (8 * (ip - ea))
    np.testing.assert_almost_equal(desp.globaltool.electrofugality, value, decimal=8)
    value = (mu - eta)**2 / (2 * eta)
    np.testing.assert_almost_equal(desp.globaltool.electrofugality, value, decimal=8)
    value = ip + 0.5 * mu * mu / eta
    np.testing.assert_almost_equal(desp.globaltool.electrofugality, value, decimal=8)


def test_analyze_ch4_fchk_exponential():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # IP = -E(HOMO) & EA = E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), 1.93295185E-01, -4.019868797400735E+01
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path, model='exponential')
    np.testing.assert_almost_equal(desp.globaltool.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(11.), energy - ea, decimal=8)
    # Check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.globaltool.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.electron_affinity, ea, decimal=8)
    # Check chemical-potential, chemical-hardness & hyper-hardness
    a, g, b = 0.30010587313, 1.03307732519, -40.4987938471
    mu, eta = -a * g, a * g**2
    print 'A, gamma, B:', a, g, desp.globaltool._B
    np.testing.assert_almost_equal(desp.globaltool.mu, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_potential, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.eta, eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_hardness, eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(2), a * (-g)**3, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(3), a * (-g)**4, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(4), a * (-g)**5, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(10), a * (-g)**11, decimal=8)
    # Check softness & hyper-softness
    np.testing.assert_almost_equal(desp.globaltool.softness, 1.0/eta, decimal=8)
    # value = 1. / (a**2 * g**3)
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(2), value, decimal=8)
    # value = -4.0 / (a**3 * g**4)
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(3), value, decimal=8)
    # Check N_max and related descriptors


def test_analyze_ch4_fchk_rational():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # IP = -E(HOMO) & EA = E(LUMO)
    ip, ea, energy = -(-5.43101269E-01), 1.93295185E-01, -4.019868797400735E+01
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path, model='rational')
    np.testing.assert_almost_equal(desp.globaltool.energy(10.), energy, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(9.), energy + ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy(11.), energy - ea, decimal=8)
    # Check ionization potential and electron affinity
    np.testing.assert_almost_equal(desp.globaltool.ip, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ionization_potential, ip, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.ea, ea, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.electron_affinity, ea, decimal=8)
    # Check chemical-potential & chemical-hardness
    a0, a1, b1 = -40.958945888, 5.16779063198, -0.126664950953
    mu = (a1 - a0 * b1) / (1 + b1 * 10)**2
    eta = 2 * b1 * (a1 - a0 * b1) / (1 + b1 * 10)**3
    np.testing.assert_almost_equal(desp.globaltool.mu, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_potential, mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.eta, eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.chemical_hardness, eta, decimal=8)
    # Check derivative of E(N)
    value = lambda N, n: math.factorial(n) * (a1 - a0 * b1) * b1**(n - 1) / (1 + b1 * N)**(n + 1)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10, 1), value(10, 1), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10, 1), mu, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10, 2), value(10, 2), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10, 2), eta, decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.6, 1), value(9.6, 1), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(9.1, 2), value(9.1, 2), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(8.5, 1), value(8.5, 1), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10.3, 1), value(10.3, 1), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10.8, 2), value(10.8, 2), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(11.5, 1), value(11.5, 1), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.energy_derivative(10.2, 3), value(10.2, 3), decimal=8)
    # Check hyper-hardness
    value = lambda n: math.factorial(n + 1) * (a1 - a0 * b1) * b1**n / (1 + b1 * 10)**(n + 2)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(2), value(2), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(3), value(3), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(4), value(4), decimal=8)
    np.testing.assert_almost_equal(desp.globaltool.hyper_hardness(5), value(5), decimal=8)
    # Check softness & hyper-softness
    np.testing.assert_almost_equal(desp.globaltool.softness, 1.0 / eta, decimal=8)
    # value =
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(2), value, decimal=8)
    # value =
    # np.testing.assert_almost_equal(desp.globaltool.hyper_softness(3), value, decimal=8)
    # Check N_max and related descriptors

def test_analyze_nci_h2o_dimer_wfn():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g.wfn')
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path)
    # testing against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    dens_cube1_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g-dens.cube')
    cube = CubeGen.from_cube(dens_cube1_path)
    # testing against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    grad_cube1_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g-grad.cube')
    dmol1 = IOData.from_file(dens_cube1_path)
    gmol1 = IOData.from_file(grad_cube1_path)

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_h2o_dimer_fchk') as dn:
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g')
        desp.compute_nci(cube2,cube=cube)
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        mol2 = IOData.from_file(cube2)

        assert abs(dmol1.coordinates - mol2.coordinates).max() < 1e-4
        assert (dmol1.numbers == mol2.numbers).all()
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        assert abs(ugrid1.grid_rvecs - ugrid2.grid_rvecs).max() < 1e-4
        assert (ugrid1.shape == ugrid2.shape).all()
        assert abs((dmol1.cube_data - mol2.cube_data)/dmol1.cube_data).max() < 1e-4
        assert abs(dmol1.pseudo_numbers - mol2.pseudo_numbers).max() < 1e-4

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)

        assert abs(gmol1.coordinates - mol2.coordinates).max() < 1e-4
        assert (gmol1.numbers == mol2.numbers).all()
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        assert abs(ugrid1.grid_rvecs - ugrid2.grid_rvecs).max() < 1e-4
        assert (ugrid1.shape == ugrid2.shape).all()
        assert abs((gmol1.cube_data - mol2.cube_data)/gmol1.cube_data).max() < 1e-4
        assert abs(gmol1.pseudo_numbers - mol2.pseudo_numbers).max() < 1e-4

def test_analyze_nci_h2o_dimer_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    file_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g.fchk')
    # Build conceptual DFT descriptor tool
    desp = Analyze_1File(file_path)
    # testing against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    dens_cube1_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g-dens.cube')
    cube = CubeGen.from_cube(dens_cube1_path)
    # testing against .cube files created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    grad_cube1_path = os.path.join(path, 'data/test/h2o_dimer_pbe_sto3g-grad.cube')
    dmol1 = IOData.from_file(dens_cube1_path)
    gmol1 = IOData.from_file(grad_cube1_path)

    with tmpdir('chemtools.analysis.test.test_base.test_analyze_nci_h2o_dimer_fchk') as dn:
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g')
        desp.compute_nci(cube2,cube=cube)
        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-dens.cube')
        mol2 = IOData.from_file(cube2)

        assert abs(dmol1.coordinates - mol2.coordinates).max() < 1e-4
        assert (dmol1.numbers == mol2.numbers).all()
        ugrid1 = dmol1.grid
        ugrid2 = mol2.grid
        assert abs(ugrid1.grid_rvecs - ugrid2.grid_rvecs).max() < 1e-4
        assert (ugrid1.shape == ugrid2.shape).all()
        assert abs((dmol1.cube_data - mol2.cube_data)/dmol1.cube_data).max() < 1e-4
        assert abs(dmol1.pseudo_numbers - mol2.pseudo_numbers).max() < 1e-4

        cube2 = '%s/%s' % (dn, 'h2o_dimer_pbe_sto3g-grad.cube')
        mol2 = IOData.from_file(cube2)

        assert abs(gmol1.coordinates - mol2.coordinates).max() < 1e-4
        assert (gmol1.numbers == mol2.numbers).all()
        ugrid1 = gmol1.grid
        ugrid2 = mol2.grid
        assert abs(ugrid1.grid_rvecs - ugrid2.grid_rvecs).max() < 1e-4
        assert (ugrid1.shape == ugrid2.shape).all()
        assert abs((gmol1.cube_data - mol2.cube_data)/gmol1.cube_data).max() < 1e-4
        assert abs(gmol1.pseudo_numbers - mol2.pseudo_numbers).max() < 1e-4
