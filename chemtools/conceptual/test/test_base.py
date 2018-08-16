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


from numpy.testing import assert_raises
from chemtools.conceptual.base import BaseGlobalTool, BaseLocalTool, BaseCondensedTool


def test_global_base_raises():
    # check invalid N0
    assert_raises(ValueError, BaseGlobalTool, -.0, 1.0)
    assert_raises(ValueError, BaseGlobalTool, 0.0, 0.1)
    # check invalid Nmax
    assert_raises(ValueError, BaseGlobalTool, 8.0, -1.0)
    assert_raises(ValueError, BaseGlobalTool, 8.0, -0.1)


def test_local_base_raises():
    # check invalid N0 & Nmax
    assert_raises(ValueError, BaseLocalTool, -1., 2.0)
    assert_raises(ValueError, BaseLocalTool, 2., -1.5)
    assert_raises(ValueError, BaseLocalTool, 2., -3.0)


def test_condensed_base_raises():
    # check invalid N0 & Nmax
    assert_raises(ValueError, BaseCondensedTool, 0., 2.0)
    assert_raises(ValueError, BaseCondensedTool, -1., 2.0)
    assert_raises(ValueError, BaseCondensedTool, 2., -1.5)
    assert_raises(ValueError, BaseCondensedTool, 2., -3.0)
    # build model & check properties/methods
    model = BaseCondensedTool(1.0, 1.5, None)
    assert_raises(NotImplementedError, model.population, 3.0)
    assert_raises(NotImplementedError, model.population, 1.0)
    assert_raises(NotImplementedError, model.population_derivative, 1.0)
    assert_raises(NotImplementedError, model.population_derivative, 1.3, 1)
    assert_raises(NotImplementedError, model.population_derivative, 2.4, 2)
