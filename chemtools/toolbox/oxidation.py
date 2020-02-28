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
"""Module for Oxidation State."""


from chemtools.wrappers.molecule import Molecule


class EOS(object):
    def __init__(self, molecule, denspart, grid):
        self.molecule = molecule
        self.denspart = denspart
        self.grid = grid

    @classmethod
    def from_molecule(cls, molecule, denspart, grid):
        """Initialize class from `Molecule` object.

        Parameters
        ----------
        molecule : `Molecule`
            Instance of `Molecular` class.
        denspart : `DensPart`
            Instance of `DensPart` class.
        grid : instance of `MolecularGrid`, optional
            Molecular numerical integration grid.

        """
        return cls(molecule, denspart, grid)

    @classmethod
    def from_file(cls, fname, denspart, grid):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's fname.
        denspart : `DensPart`
            Instance of `DensPart` class.
        grid : instance of `MolecularGrid`, optional
            Molecular numerical integration grid.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, denspart, grid)

    def compute_fragment_overlap(self, index, fragment=None):
        # compute MO overlap matrix for the atom specified with index
        pass

    def compute_oxidation_state(self, fragments=None):
        # if fragments=None, EOS is computed for each atom
        pass

    def compute_effective_orbitals(self):
        pass
