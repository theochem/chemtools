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
# pragma pylint: disable=invalid-name
"""Orbital-Based Local Tools."""


from chemtools.utils.utils import doc_inherit
from chemtools.utils.cube import UniformGrid
from chemtools.outputs.vmd import print_vmd_script_isosurface
from chemtools.wrappers.molecule import Molecule


__all__ = ['MOTBasedTool']


class MOTBasedTool(object):
    """Molecular Orbital Theory Based Descriptive Tools."""

    def __init__(self, molecule):
        r"""Initialize class from Molecule instance.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class.

        """
        self._molecule = molecule

    @classmethod
    def from_molecule(cls, molecule):
        """Initialize class from Molecule instance.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class.

        """
        return cls(molecule)

    @classmethod
    def from_file(cls, fname):
        """Initialize class from wave-function file.

        Parameters
        ----------
        fname : str
            Path to molecule's wave-function file.

        """
        molecule = Molecule.from_file(fname, wavefunction=True)
        return cls.from_molecule(molecule)

    @property
    @doc_inherit(Molecule, 'coordinates')
    def coordinates(self):
        return self._molecule.coordinates

    @property
    @doc_inherit(Molecule, 'numbers')
    def numbers(self):
        return self._molecule.numbers

    @property
    @doc_inherit(Molecule, 'nbasis')
    def nbasis(self):
        return self._molecule.nbasis

    @property
    @doc_inherit(Molecule, 'nelectrons')
    def nelectrons(self):
        return self._molecule.nelectrons

    @property
    @doc_inherit(Molecule, 'homo_index')
    def homo_index(self):
        return self._molecule.homo_index

    @property
    @doc_inherit(Molecule, 'lumo_index')
    def lumo_index(self):
        return self._molecule.lumo_index

    @property
    @doc_inherit(Molecule, 'homo_energy')
    def homo_energy(self):
        return self._molecule.homo_energy

    @property
    @doc_inherit(Molecule, 'lumo_energy')
    def lumo_energy(self):
        return self._molecule.lumo_energy

    @property
    @doc_inherit(Molecule, 'orbital_occupation')
    def orbital_occupation(self):
        return self._molecule.orbital_occupation

    @property
    @doc_inherit(Molecule, 'orbital_energy')
    def orbital_energy(self):
        return self._molecule.orbital_energy

    @property
    @doc_inherit(Molecule, 'orbital_coefficient')
    def orbital_coefficient(self):
        return self._molecule.orbital_coefficient

    @doc_inherit(Molecule, 'compute_orbital_overlap')
    def compute_orbital_overlap(self):
        return self._molecule.compute_orbital_overlap()

    @doc_inherit(Molecule, 'compute_density_matrix')
    def compute_density_matrix(self, spin='ab'):
        return self._molecule.compute_density_matrix(spin=spin)

    @doc_inherit(Molecule, 'compute_molecular_orbital')
    def compute_molecular_orbital(self, points, spin='ab', index=None):
        return self._molecule.compute_molecular_orbital(points, spin, index)

    @doc_inherit(Molecule, 'compute_density')
    def compute_density(self, points, spin='ab', index=None):
        return self._molecule.compute_density(points, spin, index)

    def compute_populations(self):
        raise NotImplementedError

    def generate_scripts(self, fname, spin='a', index=None, isosurf=0.05, grid=None):

        if spin not in ['a', 'b']:
            raise ValueError('Argument spin can only be "a" or "b".')
        if index is not None and not isinstance(index, int):
            raise ValueError('Argument index is either None or be an integer for visualization.')
        if grid is None:
            grid = UniformGrid.from_molecule(self._molecule, spacing=0.2, extension=7.0, rotate=True)
        elif not isinstance(grid, UniformGrid):
            raise ValueError('Argument grid should be a UniformGrid to generate cube files.')

        if index is None:
            spin_index = {'a': 0, 'b': 1}
            index = range(1, self.homo_index[spin_index[spin]] + 1)
        else:
            index = [index]
        for mo_index in index:
            vmdname = fname + '_mo{0}.vmd'.format(mo_index)
            cubname = fname + '_mo{0}.cube'.format(mo_index)
            mo_value = self.compute_molecular_orbital(grid.points, spin=spin, index=mo_index)
            grid.generate_cube(cubname, mo_value)
            print_vmd_script_isosurface(vmdname, cubname, isosurf=isosurf, negative=True,
                                        material='BlownGlass')
