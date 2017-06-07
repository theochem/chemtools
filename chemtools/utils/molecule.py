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
"""The Input-Output (IO) Module."""


import numpy as np

__all__ = ['BaseMolecule']


class BaseMolecule(object):
    """Molecule Class.

    This class serves as a template for the internal datatype used in ChemTools.
    Each of the methods are required to run different parts of ChemTools.
    If these methods are not defined, then a NotImplementedError will be raised to indicate that
    this method is not compatible with the given Molecule instance/class.
    """

    def __init__(self, coordinates, numbers):
        """
        Initialize class.

        Parameters
        ----------
        coordinates : ndarray
           The 2d-array containing the cartesian coordinates of atomic centers.
           It has a shape of (M, 3) where M is the number of atoms.
        numbers : ndarray
           The 1d-array containing the atomic number of atoms.
           It has a shape of (M,) where M is the number of atoms.
        """
        # bare minimum attributes to define a molecule
        if not (isinstance(coordinates, np.ndarray) and coordinates.ndim == 2):
            raise TypeError('Argument coordinates should be a 2d-array.')
        if not (isinstance(numbers, np.ndarray) and numbers.ndim == 1):
            raise TypeError('Argument numbers should be a 1d-array.')
        if coordinates.shape[0] != numbers.size:
            raise TypeError('Arguments coordinates and numbers should represent the same number '
                            'of atoms! {0} != {1}'.format(coordinates.shape[0], numbers.size))
        if coordinates.shape[1] != 3:
            raise TypeError('Argument coordinates should be a 2d-array with 3 columns.')

        self._coordinates = coordinates
        self._numbers = numbers

    @classmethod
    def from_file(cls, filename):
        """
        Initialize class given a file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        """
        raise NotImplementedError

    @property
    def coordinates(self):
        """Cartesian coordinates of atomic centers."""
        return self._coordinates

    @property
    def numbers(self):
        """Aomic number of atomic centers."""
        return self._numbers

    @property
    def natom(self):
        """Number of atoms."""
        return self.numbers.size

    @property
    def nbasis(self):
        """Number of basis functions."""
        raise NotImplementedError

    @property
    def nelectrons(self):
        """Number of alpha and beta electrons."""
        raise NotImplementedError

    @property
    def homo_index(self):
        """Index of alpha and beta HOMO orbital."""
        raise NotImplementedError

    @property
    def lumo_index(self):
        """Index of alpha and beta LUMO orbital."""
        raise NotImplementedError

    @property
    def homo_energy(self):
        """Energy of alpha and beta HOMO orbital."""
        raise NotImplementedError

    @property
    def lumo_energy(self):
        """Energy of alpha and beta LUMO orbital."""
        raise NotImplementedError

    @property
    def orbital_occupation(self):
        """Orbital occupation of alpha and beta electrons."""
        raise NotImplementedError

    @property
    def orbital_energy(self):
        """Orbital energy of alpha and beta electrons."""
        raise NotImplementedError

    @property
    def orbital_coefficient(self):
        """Orbital coefficient of alpha and beta electrons.

        The alpha and beta orbital coefficients are each storied in a 2d-array in which
        the columns represent the basis coefficients of each molecular orbital.
        """
        raise NotImplementedError

    def compute_density_matrix_array(self, spin='ab'):
        """
        Return the density matrix array for the specified spin orbitals.

        Parameters
        ----------
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        """
        return NotImplementedError

    def compute_density(self, points, spin='ab', orbital_index=None, output=None):
        """
        Return electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        return NotImplementedError

    def compute_gradient(self, points, spin='ab', orbital_index=None, output=None):
        """
        Return gradient of electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n, 3) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        raise NotImplementedError

    def compute_hessian(self, points, spin='ab', orbital_index=None, output=None):
        """
        Return hessian of electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n, 6) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        raise NotImplementedError

    def compute_esp(self, points, spin='ab', orbital_index=None, output=None):
        """
        Return the molecular electrostatic potential on the given points for the specified spin.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        raise NotImplementedError

    def compute_kinetic_energy_density(self, points, spin='ab', orbital_index=None, output=None):
        r"""
        Return positive definite kinetic energy density on the given points for the specified spin.

        Positive definite kinetic energy density is defined as,

        .. math::
           \tau \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.
           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        raise NotImplementedError
