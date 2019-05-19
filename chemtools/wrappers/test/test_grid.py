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
"""Test chemtools.wrappers.grid."""


import numpy as np
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path

from numpy.testing import assert_raises, assert_allclose

from chemtools.wrappers.beckegrid import BeckeGrid
from chemtools.wrappers.molecule import Molecule


def test_wrapper_grid_raises():
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fpath:
        assert_raises(TypeError, BeckeGrid.from_molecule, fpath, 'exp:1e-5:20:40:50')
        assert_raises(ValueError, BeckeGrid.from_file, fpath, 'ex:1e-5:20:40:50')
        assert_raises(ValueError, BeckeGrid.from_file, fpath, 'exp:1e-5:-20:40:50')
        assert_raises(ValueError, BeckeGrid.from_file, fpath, 'exp:1e-5:20:40:10')
        assert_raises(ValueError, BeckeGrid.from_file, fpath, 'pow:1e-5:20:40:50')
        assert_raises(ValueError, BeckeGrid.from_file, fpath, 'veryfin')
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fpath:
        grid = BeckeGrid.from_file(fpath)
    assert_raises(ValueError, grid.integrate, np.array([[1., 2., 3.]]))


def test_wrapper_grid_ch4():
    with path('chemtools.data', 'ch4_uhf_ccpvdz.fchk') as fpath:
        mol = Molecule.from_file(fpath)
    grid = BeckeGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'exp:1e-5:25:80:230')
    # check grid basics
    assert_allclose(mol.coordinates, grid.coordinates, rtol=0., atol=1.e-7)
    assert_allclose(mol.numbers, grid.numbers, rtol=0., atol=1.e-7)
    assert_allclose(mol.pseudo_numbers, grid.pseudo_numbers, rtol=0., atol=1.e-7)
    # check integrate
    assert_allclose(10., grid.integrate(mol.compute_density(grid.points)), rtol=0., atol=1.e-4)


def test_wrapper_grid_from_file_o2():
    with path('chemtools.data', 'o2_uhf.wfn') as fpath:
        grid = BeckeGrid.from_file(fpath, 'veryfine')
        mol = Molecule.from_file(fpath)
    # check integrate
    assert_allclose(16., grid.integrate(mol.compute_density(grid.points)), rtol=0., atol=1.e-4)
