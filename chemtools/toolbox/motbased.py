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
from chemtools.orbstools.partition import OrbitalPartitionTools
from chemtools.outputs.vmd import print_vmd_script_isosurface
from chemtools.wrappers.molecule import Molecule, MolecularOrbitals


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
        molecule = Molecule.from_file(fname)
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
    @doc_inherit(MolecularOrbitals, 'nelectrons')
    def nelectrons(self):
        return self._molecule.mo.nelectrons

    @property
    @doc_inherit(MolecularOrbitals, 'homo_index')
    def homo_index(self):
        return self._molecule.mo.homo_index

    @property
    @doc_inherit(MolecularOrbitals, 'lumo_index')
    def lumo_index(self):
        return self._molecule.mo.lumo_index

    @property
    @doc_inherit(MolecularOrbitals, 'homo_energy')
    def homo_energy(self):
        return self._molecule.mo.homo_energy

    @property
    @doc_inherit(MolecularOrbitals, 'lumo_energy')
    def lumo_energy(self):
        return self._molecule.mo.lumo_energy

    @property
    @doc_inherit(MolecularOrbitals, 'occupation')
    def orbital_occupation(self):
        return self._molecule.mo.occupation

    @property
    @doc_inherit(MolecularOrbitals, 'energy')
    def orbital_energy(self):
        return self._molecule.mo.energy

    @property
    @doc_inherit(MolecularOrbitals, 'coefficient')
    def orbital_coefficient(self):
        return self._molecule.mo.coefficient

    @doc_inherit(Molecule, 'compute_density_matrix')
    def compute_density_matrix(self, spin='ab'):
        return self._molecule.compute_density_matrix(spin=spin)

    @doc_inherit(Molecule, 'compute_molecular_orbital')
    def compute_orbital_expression(self, points, spin='ab', index=None):
        return self._molecule.compute_molecular_orbital(points, spin, index)

    @doc_inherit(Molecule, 'compute_density')
    def compute_density(self, points, spin='ab', index=None):
        return self._molecule.compute_density(points, spin, index)

    def compute_charges(self, scheme="mulliken"):
        """Return the partial charges at each atom using the given population analysis method.

        Parameters
        ----------
        scheme : {"lowdin", "mulliken"}
            Type of population analysis.
            Default is Mulliken population analysis.

        Returns
        -------
        populations : np.ndarray(N,)
            Number of electrons in each atom according the population analysis.

        Raises
        ------
        ValueError
            If scheme is not "wiberg-mayer".

        """
        coeff_ab_mo_alpha, coeff_ab_mo_beta = self._molecule.mo.coefficient
        occupations_alpha, occupations_beta = self._molecule.mo.occupation
        olp_ab_ab = self._molecule.ao.compute_overlap()
        atomic_charges = self._molecule.numbers
        num_atoms = len(self._molecule.numbers)
        ab_atom_indices = self._molecule._ind_basis_center

        orbpart_alpha = OrbitalPartitionTools(
            coeff_ab_mo_alpha, occupations_alpha, olp_ab_ab, num_atoms, ab_atom_indices
        )
        orbpart_beta = OrbitalPartitionTools(
            coeff_ab_mo_beta, occupations_beta, olp_ab_ab, num_atoms, ab_atom_indices
        )
        if scheme == "mulliken":
            pop = orbpart_alpha.mulliken_populations()
            pop += orbpart_beta.mulliken_populations()
        elif scheme == "lowdin":
            pop = orbpart_alpha.lowdin_populations()
            pop += orbpart_beta.lowdin_populations()
        else:
            raise ValueError("`scheme` must be one of 'mulliken' or 'lowdin'.")

        return atomic_charges - pop

    def compute_bond_orders(self, scheme="wiberg-mayer"):
        """Return the bond order for each pair of atoms.

        Parameters
        ----------
        scheme : "wiberg-mayer"
            Type of the bond order.
            Default is the Wiberg-Mayer bond order.

        Returns
        -------
        bond_orders : np.ndarray(N, N)
            Bond order for each atom pair.

        Raises
        ------
        ValueError
            If scheme is not "wiberg-mayer".

        """
        coeff_ab_mo_alpha, coeff_ab_mo_beta = self.orbital_coefficient
        occupations_alpha, occupations_beta = self.orbital_occupation
        olp_ab_ab = self._molecule.ao.compute_overlap()
        num_atoms = len(self._molecule.numbers)
        ab_atom_indices = self._molecule._ind_basis_center

        orbpart_alpha = OrbitalPartitionTools(
            coeff_ab_mo_alpha, occupations_alpha, olp_ab_ab, num_atoms, ab_atom_indices
        )
        orbpart_beta = OrbitalPartitionTools(
            coeff_ab_mo_beta, occupations_beta, olp_ab_ab, num_atoms, ab_atom_indices
        )

        if scheme == "wiberg-mayer":
            bond_order = orbpart_alpha.bond_order_wiberg_mayer_unrestricted
            bond_order += orbpart_beta.bond_order_wiberg_mayer_unrestricted
        else:
            raise ValueError("Bond order scheme must 'wiberg-mayer'.")

        return bond_order

    def generate_scripts(self, fname, spin='a', index=None, isosurf=0.05, grid=None):
        """Generate VMD script(s) and cube file(s) to visualize MO iso-surface of given orbitals.

        Parameters
        ----------
        fname : str
            A string representing the path to a fname of generated files.
            The VMD script and cube file will be named fname_mo{index}.vmd and
            fname_mo{index}.cube, respectively.
        spin : str, optional
           The type of occupied spin orbitals. Choose either 'a' or 'b'.
        index : int, optional
           Integer representing the index of spin orbital to visualize. Spin orbitals are each
           indexed from 1 to :attr:`nbasis`. If None, files for visualizing all orbitals are
           generated.
        isosurf : float, optional
            Value of MO iso-surface used in VMD script.
        grid : UniformGrid, optional
           Instance of UniformGrid used for computation and generating cube file(s).
           If None, a cubic grid is constructed from molecule with spacing=0.2 & extension=5.0.

        """
        if spin not in ['a', 'b']:
            raise ValueError('Argument spin can only be "a" or "b".')
        if index is not None and not isinstance(index, int):
            raise ValueError('Argument index is either None or an integer for visualization. '
                             'Given index={0}'.format(index))
        if grid is None:
            grid = UniformGrid.from_molecule(self._molecule, spacing=0.2, extension=5.0, rotate=True)
        elif not isinstance(grid, UniformGrid):
            raise ValueError('Argument grid should be a UniformGrid to generate cube files.')

        if index is None:
            spin_index = {'a': 0, 'b': 1}
            index = range(1, self._molecule.mo.homo_index[spin_index[spin]] + 1)
        else:
            index = [index]

        for mo_index in index:
            vmdname = fname + '_mo{0}.vmd'.format(mo_index)
            cubname = fname + '_mo{0}.cube'.format(mo_index)
            mo_value = self.compute_orbital_expression(grid.points, spin=spin, index=mo_index)
            grid.generate_cube(cubname, mo_value)
            print_vmd_script_isosurface(vmdname, cubname, isosurf=isosurf, negative=True,
                                        material='BlownGlass')
