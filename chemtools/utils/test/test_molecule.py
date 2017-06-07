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
from chemtools import BaseMolecule


def test_molecule_check_raises():
    # check invalid coordinates argument
    assert_raises(TypeError, BaseMolecule, [1., 2., 3.], np.array([6]))
    assert_raises(TypeError, BaseMolecule, np.array([1., 2., 3.]), np.array([6]))
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2.]]), np.array([6]))
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2., 3.], [4., 5., 6.]]), np.array([6]))
    # check invalid numbers argument
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2., 3.]]), 6)
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2., 3.]]), [6])
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2., 3.]]), np.array([[6]]))
    assert_raises(TypeError, BaseMolecule, np.array([[1., 2., 3.]]), np.array([6, 8]))
    # check invalid coordinates & numbers argument simultaneously
    assert_raises(TypeError, BaseMolecule, np.array([1., 2., 3.]), [7])
    assert_raises(TypeError, BaseMolecule, [1., 2., 3.], [6])
