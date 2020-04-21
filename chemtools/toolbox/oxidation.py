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
from chemtools.wrappers.grid import MolecularGrid
from chemtools.wrappers.part import DensPart

import numpy as np
import itertools


class EOS(object):
    def __init__(self, molecule, part):
        if not isinstance(part, DensPart):
            raise TypeError('Argument part should be an instance of DensPart class.')

        self.molecule = molecule
        self.part = part
        self._frags = None
        self._reliability = None

    @classmethod
    def from_molecule(cls, molecule, scheme, grid=None, proatomdb=None):
        """Initialize class from `Molecule` object.

        Parameters
        ----------
        molecule : `Molecule`
            Instance of `Molecular` class.
        part : `DensPart`
            Instance of `DensPart` class.
        grid : `MolecularGrid`
            Molecular numerical integration grid.

        """
        part = DensPart.from_molecule(molecule, scheme=scheme, grid=grid, local=False,
                                      proatomdb=proatomdb)
        return cls(molecule, part)

    @classmethod
    def from_file(cls, fname, scheme, grid=None, proatomdb=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's fname.
        part : `DensPart`
            Instance of `DensPart` class.
        grid : `MolecularGrid`
            Molecular numerical integration grid.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, scheme=scheme, grid=grid, proatomdb=proatomdb)

    @property
    def reliability(self):
        """float : Reliability index of oxidation state assignment."""
        return self._reliability

    def compute_fragment_overlap(self, fragments=None, spin='ab'):
        # compute MO overlap matrix for fragments

        if spin == 'a':
            ne = int(self.molecule.mo.nelectrons[0])
        elif spin == 'b':
            ne = int(self.molecule.mo.nelectrons[1])
        else:
            raise NotImplementedError('Not clear what to do here!')

        # defining fragments/atoms
        if fragments is None:
            self._frags = [[index] for index in range(len(self.molecule.numbers))]
        else:
            self._frags = fragments

        orbitals = self.molecule.compute_molecular_orbital(self.part.grid.points, spin=spin).T

        # generating qij array for each atom/fragment
        arr = np.zeros((len(self._frags), ne, ne))

        for i, j in itertools.combinations_with_replacement(range(ne), 2):
            qij = self.part.condense_to_fragments(orbitals[i] * orbitals[j], self._frags, w_power=2)
            for x in range(len(self._frags)):
                arr[x][i][j] = qij[x]
                arr[x][j][i] = qij[x]

        return arr

    def compute_fragment_occupation(self, fragments=None, spin='ab'):
        # computes effective orbitals occupation for each fragment passed

        # compute fragment overlap matrix
        arr = self.compute_fragment_overlap(fragments, spin=spin)

        # diagonalize overlap matrix
        _, s, _ = np.linalg.svd(arr)

        return s

    def compute_oxidation_state(self, fragments=None):

        # compute oxidation state for fragments

        nalpha = int(self.molecule.mo.nelectrons[0])
        nbeta = int(self.molecule.mo.nelectrons[1])

        # TODO: avoid repeated calculation for restricted case
        s_a = self.compute_fragment_occupation(fragments, spin='a')
        s_b = self.compute_fragment_occupation(fragments, spin='b')

        sorted_alpha = sorted([(s, i) for i, row in enumerate(s_a) for s in row], reverse=True)
        sorted_beta = sorted([(s, i) for i, row in enumerate(s_b) for s in row], reverse=True)

        occs_a = [item[1] for item in sorted_alpha]
        occs_b = [item[1] for item in sorted_beta]

        occs_frag_a = [occs_a[:nalpha].count(index) for index in range(len(self._frags))]
        occs_frag_b = [occs_b[:nbeta].count(index) for index in range(len(self._frags))]

        z_frag = [sum([self.molecule.numbers[index] for index in frag]) for frag in self._frags]
        oxidation = np.array(z_frag) - np.array(occs_frag_a) - np.array(occs_frag_b)

        # Reliability index
        lumo_a = nalpha
        self._reliability = 100.00
        if len(self._frags) != 1:
            while sorted_alpha[nalpha-1][1] == sorted_alpha[lumo_a][1]:
                lumo_a += 1

            self._reliability = 100 * min(1, sorted_alpha[nalpha-1][0] - sorted_alpha[lumo_a][0] + 0.5)

        return oxidation

    def compute_effective_orbital(self):
        pass
