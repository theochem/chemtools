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

from horton import ProAtomDB
from horton.scripts.wpart import wpart_schemes

from chemtools.wrappers.molecule import Molecule
from chemtools.wrappers.grid import MolecularGrid


__all__ = ['DensPart']


class DensPart(object):
    """Density-based Atoms-in-Molecules Partitioning Class."""

    def __init__(self, coordinates, numbers, pseudo_numbers, density, grid, scheme="h", **kwargs):
        """Initialize class.

        Parameters
        ----------
        coordinates : np.ndarray, shape=(M, 3)
            Cartesian coordinates of `M` atoms in the molecule.
        numbers : np.ndarray, shape=(M,)
            Atomic number of `M` atoms in the molecule.
        pseudo_numbers : np.ndarray, shape=(M,)
            Pseudo-number of `M` atoms in the molecule.
        density : np.ndarray, shape=(N,)
            Total density to be partitioned.
        grid : BeckeMolGrid
            Instance of BeckeMolGrid numerical integration grid.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.

        """
        wpart = wpart_schemes[scheme]
        # make proatom database
        if scheme.lower() not in ["mbis", "b"]:
            if "proatomdb" not in kwargs.keys() or kwargs['proatomdb'] is None:
                proatomdb = ProAtomDB.from_refatoms(numbers)
                kwargs["proatomdb"] = proatomdb
            kwargs["local"] = False
        # partition
        self.part = wpart(coordinates, numbers, pseudo_numbers, grid, density, **kwargs)
        self.part.do_charges()

        self.grid = grid
        self.density = density
        self.coordines = coordinates
        self.numbers = numbers
        self.pseudo_numbers = pseudo_numbers
        self.charges = self.part['charges']

    @classmethod
    def from_molecule(cls, mol, scheme=None, grid=None, spin="ab", **kwargs):
        """Initialize class given a Molecule instance.

        Parameters
        ----------
        mol : Molecule
            Instance of Molecule class.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.
        grid : MolecularGrid
            Instance of MolecularGrid numerical integration grid.
        spin : str
           Type of occupied spin orbitals; choose either "a" (for alpha), "b" (for beta),
           and "ab" (for alpha + beta).

        """
        if grid is None:
            grid = MolecularGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                 specs="fine", rotate=False, k=3)
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
        """Initialize class given a wave-function file.

        Parameters
        ----------
        fname : str
            Path to molecule's files.
        scheme : str
            Type of atoms-in-molecule partitioning scheme.
        grid : MolecularGrid
            Instance of MolecularGrid integration grid.
        spin : str
           Type of occupied spin orbitals; choose either "a" (for alpha), "b" (for beta),
           and "ab" (for alpha + beta).

        """
        mol = Molecule.from_file(fname)
        return cls.from_molecule(mol, scheme=scheme, grid=grid, spin=spin, **kwargs)

    def condense_to_atoms(self, value, w_power=1):
        condensed = np.zeros(self.part.natom)
        for index in range(self.part.natom):
            at_grid = self.part.get_grid(index)
            at_weight = self.part.cache.load("at_weights", index)
            # wcor = self.part.get_wcor(index)
            local_prop = self.part.to_atomic_grid(index, value)
            condensed[index] = at_grid.integrate(at_weight**w_power * local_prop)
        return condensed

    def condense_to_fragments(self, value, fragments=None, w_power=1):
        if fragments is None:
            fragments = [[index] for index in range(self.part.natom)]
        else:
            # check fragments
            segments = sorted([item for frag in fragments for item in frag])
            segments = np.array(segments)
            if segments.size != self.numbers.size:
                raise ValueError("Items in Fragments should uniquely represent all atoms.")
        condensed = np.zeros(len(fragments))
        for index, frag in enumerate(fragments):
            weight = np.zeros(self.grid.points.shape[0])
            for item in frag:
                weight += self.part.cache.load("at_weights", item)
            share = self.grid.integrate(weight ** w_power, value)
            condensed[index] = share
        return condensed


def check_molecule_grid(mol, grid):
    """Check that molecule and grid represent the same system.

    Parameters
    ----------
    mol : Molecule
        Instance of Molecule class.
    grid : MolecularGrid
        Instance of MolecularGrid numerical integration grid.

    """
    if not np.max(abs(grid.centers - mol.coordinates)) < 1.e-6:
        raise ValueError("Argument molecule & grid should have the same coordinates/centers.")
    if not np.max(abs(grid.numbers - mol.numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same numbers.")
    if not np.max(abs(grid.pseudo_numbers - mol.pseudo_numbers)) < 1.e-6:
        raise ValueError("Arguments molecule & grid should have the same pseudo_numbers.")
