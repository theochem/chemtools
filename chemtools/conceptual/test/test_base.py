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
"""Test chemtools.conceptual.base Module."""

import numpy as np
from numpy.testing import assert_raises
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool


def test_global_base_raises():
    # check invalid N0
    assert_raises(ValueError, BaseGlobalTool, -3.0, 1.0)
    assert_raises(ValueError, BaseGlobalTool, 0.0, 0.1)
    # check invalid Nmax
    assert_raises(ValueError, BaseGlobalTool, 8.0, -1.0)
    assert_raises(ValueError, BaseGlobalTool, 8.0, -0.1)


def test_local_base_raises():
    # fake density arrays
    d0 = np.array([1.0, 3.0, 5.0, 2.0, 7.0])
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, 3.0, 2.0, 8.0])
    # check invalid N
    assert_raises(ValueError, BaseLocalTool, {-1.: d0, 0.0: dp, 1.: dm}, -1.)
    assert_raises(ValueError, BaseLocalTool, {0.5: d0, -1.: dp, -1.5: dm}, 0.5)
    # check invalid density arrays
    assert_raises(ValueError, BaseLocalTool, {10.: -1 * d0, 11.: dp, 9.: dm}, 10)
    assert_raises(ValueError, BaseLocalTool, {8.5: d0, 9.5: -1 * dp, 7.5: dm}, 8.5)
    assert_raises(ValueError, BaseLocalTool, {9.1: d0, 10.1: dp, 8.1: -1 * dm}, 9.1)
    dp = np.array([0.5, -4.5, 6.0, 1.0, 5.0])
    assert_raises(ValueError, BaseLocalTool, {5.5: d0, 6.5: dp, 4.5: dm}, 5.5)
    dp = np.array([0.5, 4.5, 6.0, 1.0, 5.0])
    dm = np.array([1.0, 4.0, -3.0, -2.0, -8.0])
    assert_raises(ValueError, BaseLocalTool, {5.5: d0, 6.5: dp, 4.5: dm}, 5.5)
