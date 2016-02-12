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
'''Analyze Output Module.'''


import numpy as np
from horton import BeckeMolGrid
from chemtools.tool.globaltool import QuadraticGlobalTool
from chemtools.tool.localtool import QuadraticLocalTool


class Analyze(object):
    '''
    Conceptual DFT Ananlysis of molecule.
    '''
    def __init__(self, molecule, approx='FMO', model='quadratic', grid=None):
        '''
        Parameters
        ----------
        mol :
            instance of ``horton.IOData``.
        grid :
            instance of ``horton.BeckeMolGrid``.
        '''
        self._mol = molecule
        if grid is None:
            # Set up the default grid
            self._grid = BeckeMolGrid(self._mol.coordinates, self._mol.numbers,
                                      self._mol.pseudo_numbers, 'fine',
                                      random_rotate=False, mode='keep')
        else:
            # Check the consistency of grid and molecule
            assert (abs(grid.numbers - self._mol.numbers) < 1.e-6).all()
            assert (abs(grid.coordinates - self._mol.coordinates) < 1.e-6).all()
            self._grid = grid

        # Get energy and density of HOMO and LUMO
        homo_energy = self._mol.exp_alpha.homo_energy
        lumo_energy = self._mol.exp_alpha.lumo_energy
        homo_index = self._mol.exp_alpha.get_homo_index()
        lumo_index = self._mol.exp_alpha.get_lumo_index()
        homo_density = self._mol.obasis.compute_grid_orbitals_exp(
            self._mol.exp_alpha, self._grid.points, np.array([homo_index]))**2
        lumo_density = self._mol.obasis.compute_grid_orbitals_exp(
            self._mol.exp_alpha, self._grid.points, np.array([lumo_index]))**2

        # Make instances of the ConceptualTools
        # self._glob = QuadraticGlobalTool(0, ip=-homo_energy, ea=-lumo_energy)
        self._glob = QuadraticGlobalTool(0, lumo_energy, -homo_energy)
        self._local = QuadraticLocalTool(ff_plus=homo_density,
                                          ff_minus=lumo_density,
                                          global_instance=self._glob)

    @property
    def mol(self):
        ''''Instance of ``IOData`` class.'''
        return self._mol

    @property
    def grid(self):
        ''''Instance of ``BeckeMolGrid`` class.'''
        return self._grid

    @property
    def glob(self):
        '''Instance of ``GlobalConceptualTool`` class.'''
        return self._glob

    @property
    def local(self):
        '''Instance of ``LocalConceptualTool`` class.'''
        return self._local
