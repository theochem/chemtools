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
"""Orbital-Based Population Analysis."""


import numpy as np

from chemtools.utils.utils import doc_inherit
from chemtools.orbstools.partition import OrbitalPartitionTools
from chemtools.wrappers.molecule import Molecule, MolecularOrbitals


class OrbPart(object):
    """Orbital Partitioning or Population Analysis Class."""

    def __init__(self, mo, ao, numbers, scheme='mulliken'):
        r"""Initialize class from Molecule instance.

        Parameters
        ----------
        mo : MolecularOrbitals
            An instance of `MolecularOrbitals` class.
        ao : AtomicOrbitals
            An instance of `AtomicOrbitals` class.
        numbers : np.ndarray
            Atomic numbers of orbital centers.
        scheme : str
            Type of population analysis scheme.

        """
        self._mo = mo
        self._ao = ao
        self._numbers = numbers
        self._scheme = scheme
        # compute weights
        if scheme == 'mulliken':
            weight = self._compute_weights_mulliken()
        elif scheme == 'lowdin':
            raise NotImplementedError
        else:
            raise ValueError("`scheme` must be one of 'mulliken' or 'lowdin'.")
        # compute atomic populations
        dm = self._mo.compute_dm()
        olp = self._ao.compute_overlap()
        self._pops = np.zeros(len(self._numbers))
        for i in range(len(self._numbers)):
            self._pops[i] = np.sum(olp * weight[i] * dm)
        self._charges = self._numbers - self._pops

    @classmethod
    def from_molecule(cls, molecule, scheme='mulliken'):
        """Initialize class from Molecule instance.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class.
        scheme : str
            Type of population analysis scheme.

        """
        return cls(molecule.mo, molecule.ao, molecule.numbers, scheme=scheme)

    @classmethod
    def from_file(cls, fname, scheme='mulliken'):
        """Initialize class from wave-function file.

        Parameters
        ----------
        fname : str
            Path to molecule's wave-function file.
        scheme : str
            Type of population analysis scheme.

        """
        return cls.from_molecule(Molecule.from_file(fname), scheme=scheme)

    @property
    @doc_inherit(Molecule, 'numbers')
    def numbers(self):
        return self._numbers

    @property
    def scheme(self):
        return self._scheme

    @property
    def charges(self):
        """"""
        return self._charges

    @property
    def populations(self):
        """np.ndarray(N,) : Number of electrons in each atom."""
        return self._pops

    def _compute_weights_mulliken(self):
        """Return the weights defined based on Mulliken population analysis."""
        ao_center_index = self._ao.center_index
        nbasis = self._ao.nbasis
        w = np.zeros((len(self._numbers), nbasis, nbasis))
        for index, arr in enumerate(w):
            for i in range(nbasis):
                for j in range(nbasis):
                    ci, cj = ao_center_index[i], ao_center_index[j]
                    if ci == cj == index:
                        arr[i, j] = 1.0
                    elif ci == index or cj == index:
                        arr[i, j] = 0.5
        return w

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
        coeff_ab_mo_alpha, coeff_ab_mo_beta = self._mo.coefficient
        occupations_alpha, occupations_beta = self._mo.occupation
        olp_ab_ab = self._ao.compute_overlap()
        num_atoms = len(self._numbers)
        ab_atom_indices = self._ao.center_index

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

    # def generate_scripts(self, fname, spin='a', index=None, isosurf=0.05, grid=None):
    #     """Generate VMD script(s) and cube file(s) to visualize MO iso-surface of given orbitals.
    #
    #     Parameters
    #     ----------
    #     fname : str
    #         A string representing the path to a fname of generated files.
    #         The VMD script and cube file will be named fname_mo{index}.vmd and
    #         fname_mo{index}.cube, respectively.
    #     spin : str, optional
    #        The type of occupied spin orbitals. Choose either 'a' or 'b'.
    #     index : int, optional
    #        Integer representing the index of spin orbital to visualize. Spin orbitals are each
    #        indexed from 1 to :attr:`nbasis`. If None, files for visualizing all orbitals are
    #        generated.
    #     isosurf : float, optional
    #         Value of MO iso-surface used in VMD script.
    #     grid : UniformGrid, optional
    #        Instance of UniformGrid used for computation and generating cube file(s).
    #        If None, a cubic grid is constructed from molecule with spacing=0.2 & extension=5.0.
    #
    #     """
    #     if spin not in ['a', 'b']:
    #         raise ValueError('Argument spin can only be "a" or "b".')
    #     if index is not None and not isinstance(index, int):
    #         raise ValueError('Argument index is either None or an integer for visualization. '
    #                          'Given index={0}'.format(index))
    #     if grid is None:
    #         grid = UniformGrid.from_molecule(self._molecule, spacing=0.2, extension=5.0, rotate=True)
    #     elif not isinstance(grid, UniformGrid):
    #         raise ValueError('Argument grid should be a UniformGrid to generate cube files.')
    #
    #     if index is None:
    #         spin_index = {'a': 0, 'b': 1}
    #         index = range(1, self._mo.homo_index[spin_index[spin]] + 1)
    #     else:
    #         index = [index]
    #
    #     for mo_index in index:
    #         vmdname = fname + '_mo{0}.vmd'.format(mo_index)
    #         cubname = fname + '_mo{0}.cube'.format(mo_index)
    #         mo_value = self.compute_orbital_expression(grid.points, spin=spin, index=mo_index)
    #         grid.generate_cube(cubname, mo_value)
    #         print_vmd_script_isosurface(vmdname, cubname, isosurf=isosurf, negative=True,
    #                                     material='BlownGlass')
