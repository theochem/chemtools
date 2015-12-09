#!/usr/bin/env python


import os
from horton import IOData
from chemtools import *


def test_analyze_ch4_fchk():
    # Temporary trick to find the data files
    path = os.path.abspath(os.path.dirname(__file__)).rsplit('/', 3)[0]
    path_file = os.path.join(path, 'data/test/ch4_uhf_ccpvdz.fchk')
    # Get CH4 molecule
    mol = IOData.from_file(path_file)
    # Conceptual DFT analysis
    anal = Analyze(mol)
    # IP = -E(HOMO)
    assert abs(anal.glob.ip - 0.543101269) < 1.e-6
    assert abs(anal.glob.ionization_potential - 0.543101269) < 1.e-6
    # EA = -E(LUMO)
    assert abs(anal.glob.ea - (-0.193295185)) < 1.e-6
    assert abs(anal.glob.electron_affinity - (-0.193295185)) < 1.e-6
    # Electrophilicity
    value = (-0.5 * (0.543101269 - 0.193295185))**2
    value /= (2 * (0.543101269  + 0.193295185 ))
    assert abs(anal.glob.electrophilicity - value) < 1.e-6
    # Fukui Function arrays has the same number of points as the grid
    assert anal.local.ff_plus.shape == (len(anal.grid.points), 1)
    assert anal.local.ff_minus.shape == (len(anal.grid.points), 1)
    assert anal.local.ff_central.shape == (len(anal.grid.points), 1)
    assert anal.local.dual_descriptor.shape == (len(anal.grid.points), 1)
    # Fukui Function integrates to 1.0
    assert abs(anal.grid.integrate(anal.local.ff_plus) - 1.0) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.ff_minus) - 1.0) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.ff_central) - 1.0) < 1.e-4
    # Dual Descriptor integrates to 0.0
    assert abs(anal.grid.integrate(anal.local.dual_descriptor)) < 1.e-4
    # Local Descriptors integrate to the value of global descriptor
    # Local ionization potential
    assert abs(anal.grid.integrate(anal.local.ip_plus) - 0.543101269) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.ip_minus) - 0.543101269) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.ip_central) - 0.543101269) < 1.e-4
    # Local electrophilicity
    value = (-0.5 * (0.543101269 - 0.193295185))**2
    value /= (2 * (0.543101269  + 0.193295185 ))
    assert abs(anal.grid.integrate(anal.local.electrophilicity_plus) - value) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.electrophilicity_minus) - value) < 1.e-4
    assert abs(anal.grid.integrate(anal.local.electrophilicity_central) - value) < 1.e-4
