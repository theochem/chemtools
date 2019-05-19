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

from unittest import TestCase

import numpy as np
from horton import BeckeMolGrid
from importlib_resources import path

from chemtools.wrappers.beckegrid import BeckeGrid
from chemtools.wrappers.molecule import Molecule




class TestGrid(TestCase):
    """Grid test class."""

    def setUp(self):
        """Set up function for all tests."""
        with path('chemtools.data', 'water_b3lyp_sto3g.fchk') as file_path:
            self.mol = Molecule.from_file(file_path)
        self.grid_mol = BeckeGrid(self.mol.coordinates, self.mol.numbers, self.mol.pseudo_numbers)

    def test_string_init(self):
        """Test Grid object initialization."""
        grid = BeckeGrid(self.mol.coordinates, self.mol.numbers, self.mol.pseudo_numbers,
                         specification='exp:1e-5:20:40:50')
        # assert grid._grid_type is None
        assert grid.rrad == 40
        assert grid.rpoint == 50

    def test_init(self):
        """Check initial values and attributs."""
        assert isinstance(self.grid_mol, BeckeGrid)
        assert np.allclose(self.grid_mol._coordinates, self.mol.coordinates)
        assert np.allclose(self.grid_mol._numbers, self.mol.numbers)
        assert np.allclose(self.grid_mol._pseudo_n, self.mol.pseudo_numbers)
        assert self.grid_mol._grid_type == 'medium'
        assert self.grid_mol._custom_type == [None, None, None, None, None]
        assert self.grid_mol.grid_type == 'medium'
        assert self.grid_mol._k == 3
        assert self.grid_mol._random_rotate is False
        assert self.grid_mol._mode == "discard"

    def test_error_type(self):
        """Check invalid string format."""
        with self.assertRaises(ValueError):
            self.grid_mol.grid_type = 'linear:1e-5:2e1:120'
        with self.assertRaises(ValueError):
            self.grid_mol.grid_type = 'coarses'
        with self.assertRaises(ValueError):
            self.grid_mol.grid_type = 'lin:1e-5:2e1:120:110'
        # fail to converge str(xxx) to float
        with self.assertRaises(ValueError):
            self.grid_mol.grid_type = 'linear:xxxx:2e1:120:110'
        with self.assertRaises(TypeError):
            self.grid_mol.rrad = 120.
        with self.assertRaises(TypeError):
            self.grid_mol.rpoint = 120.
        with self.assertRaises(ValueError):
            self.grid_mol.rpoint = 119
        with self.assertRaises(TypeError):
            self.grid_mol.rrange = ('start', 'end')

    def test_switching_format(self):
        """Change grid format between two styles."""
        # set grid with string
        ref_grid = 'linear:1e-5:2e1:120:110'
        self.grid_mol.grid_type = ref_grid
        # check return grid information
        contents = self.grid_mol.grid_type.split(':')
        assert str(contents[0]) == 'linear'
        assert float(contents[1]) == 1e-5
        assert float(contents[2]) == 20
        assert int(contents[3]) == 120
        assert int(contents[4]) == 110
        # assert self.grid_mol._grid_type is None
        assert self.grid_mol.rname == 'linear'
        assert np.allclose(self.grid_mol.rrange, (1e-5, 20))
        assert self.grid_mol.rrad == 120
        assert self.grid_mol.rpoint == 110

        # switch back to test grid
        self.grid_mol.grid_type = 'fine'
        assert self.grid_mol.grid_type == 'fine'
        self.grid_mol.rname = 'exp'
        self.grid_mol.rrange = [0.0005, 20]
        # test not have enough info for generate a grid
        with self.assertRaises(ValueError):
            grid = self.grid_mol.grid
        self.grid_mol.rrad = 40
        self.grid_mol.rpoint = 50

        # compare two grids
        grid = self.grid_mol.grid
        ref_grid = BeckeGrid.from_molecule(self.mol, 'exp:5e-4:2e1:40:50').grid
        assert np.array_equal(grid.points, ref_grid.points)
        assert np.array_equal(grid.weights, ref_grid.weights)
        # assert self.grid_mol._grid_type is None

    def test_form_grid(self):
        """Generate the same grid."""
        gene_grid = self.grid_mol.grid
        ref_grid = BeckeMolGrid(
            self.mol.coordinates,
            self.mol.numbers,
            self.mol.pseudo_number,
            random_rotate=False)
        assert np.array_equal(gene_grid.points, ref_grid.points)
        assert np.array_equal(gene_grid.weights, ref_grid.weights)

        # generate a grid from str and compare
        self.grid_mol.grid_type = 'linear:1e-5:2e1:50:50'
        gene_grid = self.grid_mol.grid
        ref_grid = BeckeMolGrid(
            self.mol.coordinates,
            self.mol.numbers,
            self.mol.pseudo_number,
            agspec='linear:1e-5:2e1:50:50',
            random_rotate=False)
        assert np.array_equal(gene_grid.points, ref_grid.points)
        assert np.array_equal(gene_grid.weights, ref_grid.weights)

        # change back to preset grid
        self.grid_mol.grid_type = 'coarse'
        gene_grid = self.grid_mol.grid
        ref_grid = BeckeMolGrid(
            self.mol.coordinates,
            self.mol.numbers,
            self.mol.pseudo_number,
            agspec='coarse',
            random_rotate=False)
        assert np.array_equal(gene_grid.points, ref_grid.points)
        assert np.array_equal(gene_grid.weights, ref_grid.weights)

        # specify by each elements
        self.grid_mol.rname = 'exp'
        self.grid_mol.rrange = [1e-4, 18]
        self.grid_mol.rrad = 30
        self.grid_mol.rpoint = 38
        gene_grid = self.grid_mol.grid
        ref_grid = BeckeMolGrid(
            self.mol.coordinates,
            self.mol.numbers,
            self.mol.pseudo_number,
            agspec='exp:1e-4:18:30:38',
            random_rotate=False)
        assert np.array_equal(gene_grid.points, ref_grid.points)
        assert np.array_equal(gene_grid.weights, ref_grid.weights)

    def test_grid_from_molecule(self):
        """Test from_molecule generate the same grid."""
        ref_grid = BeckeGrid.from_molecule(self.mol).grid
        gene_grid = self.grid_mol.grid
        assert np.array_equal(gene_grid.points, ref_grid.points)
        assert np.array_equal(gene_grid.weights, ref_grid.weights)

        # set to wrong mode
        with self.assertRaises(ValueError):
            BeckeGrid.from_molecule(self.mol, mode='random')
        with self.assertRaises(TypeError):
            BeckeGrid.from_molecule('random')
