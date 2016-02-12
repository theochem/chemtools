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

import numpy as np
from chemtools import *


def test_local():
    # Empty instance of local descriptors
    loc = QuadraticLocalTool()
    assert loc.ff_plus is None
    assert loc.ff_minus is None
    assert loc.ff_zero is None
    assert loc.dual_descriptor is None
    # Basic instance of local descriptors w/ FF+
    loc = QuadraticLocalTool(ff_plus=np.arange(10))
    assert loc.ff_minus is None
    assert loc.ff_zero is None
    assert loc.dual_descriptor is None
    assert (abs(loc.ff_plus - np.arange(10)) < 1.e-6).all()
    # Basic instance of local descriptors w/ FF-
    loc = QuadraticLocalTool(ff_minus=np.arange(10))
    assert loc.ff_plus is None
    assert loc.ff_zero is None
    assert loc.dual_descriptor is None
    assert (abs(loc.ff_minus - np.arange(10)) < 1.e-6).all()


def test_local_H_fake_FF():
    # H atom: IP=13.59843401, EA=0.754195
    glob = QuadraticGlobalTool(0, -0.754195, 13.59843401)
    loc = QuadraticLocalTool(ff_plus=np.array([ 1.1, 2.3, -3.7]),
                              ff_minus=np.array([-4.2, 5.0, 9.4]),
                              global_instance=glob)
    # check the fukui functions
    assert (abs(loc.ff_plus - np.array([ 1.1, 2.3, -3.7])) < 1.e-6).all()
    assert (abs(loc.ff_minus - np.array([-4.2, 5.0, 9.4])) < 1.e-6).all()
    assert (abs(loc.ff_zero - np.array([-1.55, 3.65, 2.85])) < 1.e-6).all()
    assert (abs(loc.dual_descriptor - np.array([5.3, -2.7, -13.1])) < 1.e-6).all()
    # check local ionization potential
    expected = 13.59843401 * np.array([ 1.1, 2.3, -3.7])
    assert (abs(loc.ip_plus - expected) < 1.e-6).all()
    assert (abs(loc.ionization_potential_plus - expected) < 1.e-6).all()
    expected = 13.59843401 * np.array([-4.2, 5.0, 9.4])
    assert (abs(loc.ip_minus - expected) < 1.e-6).all()
    assert (abs(loc.ionization_potential_minus - expected) < 1.e-6).all()
    expected = 13.59843401 * np.array([-1.55, 3.65, 2.85])
    assert (abs(loc.ip_zero - expected) < 1.e-6).all()
    assert (abs(loc.ionization_potential_zero - expected) < 1.e-6).all()
    # check local chemical hardness
    global_chemical_hardness = 13.59843401 - 0.754195
    expected = global_chemical_hardness * np.array([ 1.1, 2.3, -3.7])
    assert (abs(loc.eta_plus - expected) < 1.e-6).all()
    assert (abs(loc.chemical_hardness_plus - expected) < 1.e-6).all()
    expected = global_chemical_hardness * np.array([-4.2, 5.0, 9.4])
    assert (abs(loc.eta_minus - expected) < 1.e-6).all()
    assert (abs(loc.chemical_hardness_minus - expected) < 1.e-6).all()
    expected = global_chemical_hardness * np.array([-1.55, 3.65, 2.85])
    assert (abs(loc.eta_zero - expected) < 1.e-6).all()
    assert (abs(loc.chemical_hardness_zero - expected) < 1.e-6).all()
    # check local nucleofugality
    global_nucleofugality  = (13.59843401 - 3 * 0.754195)**2
    global_nucleofugality /= (8 * (13.59843401 - 0.754195))
    expected = global_nucleofugality * np.array([ 1.1, 2.3, -3.7])
    assert (abs(loc.nucleofugality_plus - expected) < 1.e-6).all()
    expected = global_nucleofugality * np.array([-4.2, 5.0, 9.4])
    assert (abs(loc.nucleofugality_minus - expected) < 1.e-6).all()
    expected = global_nucleofugality * np.array([-1.55, 3.65, 2.85])
    assert (abs(loc.nucleofugality_zero - expected) < 1.e-6).all()
