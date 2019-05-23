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
"""Wrapper of Part Module."""


import numpy as np

from horton import BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes

from chemtools.wrappers.molecule import Molecule


__all__ = ['DensPart']


class DensPart(object):

    def __init__(self, coordinates, numbers, pseudo_numbers, density, grid, scheme="h", **kwargs):

        wpart = wpart_schemes[scheme]
        # make proatom database
        if scheme.lower() not in ["mbis", "b"]:
            if "proatomdb" not in kwargs.keys():
                proatomdb = ProAtomDB.from_refatoms(numbers)
            kwargs["proatomdb"] = proatomdb
        # partition
        self.part = wpart(coordinates, numbers, pseudo_numbers, grid, density, **kwargs)
        self.part.do_all()

        self.grid = grid
        self.density = density
        self.coordines = coordinates
        self.numbers = numbers
        self.pseudo_numbers = pseudo_numbers
        self.charges = self.part['charges']

    @classmethod
    def from_molecule(cls, mol, scheme=None, grid=None, spin="ab", **kwargs):
        if grid is None:
            grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                agspec="fine", random_rotate=False, mode='keep')
        else:
            check_molecule_grid(mol, grid)
        # compute molecular electron density
        dens = mol.compute_density(grid.points, spin=spin)
        if mol.pesudo_numbers is None:
            pesudo_numbers = mol.numbers
        else:
            pesudo_numbers = mol.pesudo_numbers
        return cls(mol.coordinates, mol.numbers, pesudo_numbers, dens, grid, scheme, **kwargs)

    @classmethod
    def from_file(cls, fname, scheme=None, grid=None, spin="ab", **kwargs):
        """Initialize class given a file.

        Parameters
        ----------
        fname : str
            Path to molecule's files.
        grid
        spin

        Returns
        -------

        """
        # load molecule & make/check grid
        mol = Molecule.from_file(fname)
        print('MOL = ', mol)
        return cls.from_molecule(mol, scheme=scheme, grid=grid, spin=spin, **kwargs)

    def condense_to_atoms(self, property):
        condensed = np.zeros(self.part.natom)
        for index in range(self.part.natom):
            at_grid = self.part.get_grid(index)
            at_weight = self.part.cache.load("at_weights", index)
            wcor = self.part.get_wcor(index)
            local_prop = self.part.to_atomic_grid(index, property)
            condensed[index] = at_grid.integrate(at_weight, local_prop, wcor)
        return condensed


def check_molecule_grid(molecule, grid):
    """

    Parameters
    ----------
    molecule
    grid

    Returns
    -------

    """
    if not np.max(abs(grid.coordinates - molecule.coordinates)) < 1.e-6:
        raise ValueError("Argument molecule & grid should have the same coordinates.")
    if not np.max(abs(grid.numbers - molecule.numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same numbers.")
    if not np.max(abs(grid.pseudo_numbers - molecule.pseudo_numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same pseudo_numbers.")
