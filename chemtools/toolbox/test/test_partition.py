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
# pragma pylint: disable=invalid-name
"""Test chemtools.toolbox.motbased.MOTBasedTool.compute_populations"""


import numpy as np
from numpy.testing import assert_raises
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path

from chemtools.toolbox.motbased import MOTBasedTool


def test_populations_mulliken_h2o():
    # Test against Gaussian
    with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as fname:
        mot = MOTBasedTool.from_file(str(fname))
    mulliken = np.array([-4.32227787E-01, 2.16114060E-01, 2.16113727E-01])
    assert np.allclose(mulliken, mot.compute_charges(), atol=1e-6)
    assert_raises(ValueError, mot.compute_charges, "bad type")


def test_populations_mulliken_h2o_cation():
    # Test against Gaussian
    with path("chemtools.data", "h2o_q+1_ub3lyp_ccpvtz.fchk") as fname:
        mot = MOTBasedTool.from_file(str(fname))
    mulliken = np.array([3.49417097E-01, 3.25291762E-01, 3.25291141E-01])
    assert np.allclose(mulliken, mot.compute_charges(), atol=1e-6)


def test_populations_mulliken_h2o_anion():
    # Test against Gaussian
    with path("chemtools.data", "h2o_q-1_ub3lyp_ccpvtz.fchk") as fname:
        mot = MOTBasedTool.from_file(str(fname))
    mulliken = np.array([-2.64833827E-01, -3.67583325E-01, -3.67582849E-01])
    assert np.allclose(mulliken, mot.compute_charges(), atol=1e-6)


def test_compute_bond_order():
    """Test MOTBasedTool.compute_bond_order.

    Results are compared with those generated from Gaussian. The system is H2O UB3LYP/aug-cc-pVDZ,
    singlet, and at zero net charge. The following the its coordinates (au):

    O 0.0159484498, 0.0170042791, 0.0238579956
    H -0.772778442, 0.561446550, 1.57501231
    H 1.29850109, 1.26951236, -0.309113326

    """
    with path("chemtools.data.examples", "h2o.fchk") as fname:
        mot = MOTBasedTool.from_file(str(fname))

    bond_order = np.array(
        [[0.0, 1.059127, 1.059127], [1.059127, 0.0, -0.008082], [1.059127, -0.008082, 0.0]]
    )

    assert np.allclose(bond_order, mot.compute_bond_orders(), atol=1e-6)

    assert_raises(ValueError, mot.compute_bond_orders, "bad type")
