# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2021 The ChemTools Development Team
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
"""Wrapper Module."""

import logging
import numpy as np
from iodata import load_one
from gbasis.wrappers import from_iodata
from gbasis.integrals.overlap import overlap_integral
from gbasis.evals.eval import evaluate_basis
from gbasis.evals.density import (evaluate_density, evaluate_density_gradient,
                                  evaluate_density_hessian, evaluate_posdef_kinetic_energy_density)
from gbasis.evals.electrostatic_potential import electrostatic_potential
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


__all__ = ["Molecule"]


class Molecule:
    """ Molecule class based on IOData, gbasis, and grid"""

    def __init__(self, iodata):
        """
        Initialize class.

        Parameters
        ----------
        iodata : iodata.iodata.IOData
           An instance of iodata.iodata.IOData object.
        """
        self._iodata = iodata
        self._coords = self._iodata.atcoords
        self._numbers = self._iodata.atnums
        self._pseudo_numbers = self._iodata.atcorenums
        self._charges = self._iodata.atcharges

        # get the Atomic Orbital
        if hasattr(self._iodata, "obasis"):
            try:
                self._ao = AtomicOrbitals.from_molecule(self._iodata)
            except AttributeError:
                self._ao = None
        else:
            self._ao = None

        # get the Molecular Orbital
        if hasattr(self._iodata, "mo"):
            self._mo = MolecularOrbitals.from_molecule(self._iodata)
        else:
            self._mo = None

    @classmethod
    def from_file(cls, fname):
        try:
            iodata = load_one(str(fname))
        except IOError as _:
            try:
                with path("chemtools.data.examples", str(fname)) as fname:
                    logging.info("Loading {0}".format(str(fname)))
                    iodata = load_one(str(fname))
            except IOError as error:
                logging.info(error)
        return cls(iodata)

    def __getattr__(self, attr):
        return getattr(self._iodata, attr, None)

    @property
    def coordinates(self):
        """Cartesian coordinates of atomic centres"""
        return self._coords

    @property
    def numbers(self):
        """Atomic numbers of atomic centres."""
        return self._numbers

    @property
    def pseudo_numbers(self):
        """Pseudo numbers of atomic centres."""
        return self._pseudo_numbers

    @property
    def charges(self):
        """Return a dict with charges"""
        return self._charges

    @property
    def nbasis(self):
        return self._iodata.mo.nbasis

    @property
    def ao(self):
        """Atomic orbital instance"""
        return self._ao

    @property
    def mo(self):
        """Molecular orbital instance."""
        return self._mo

    def compute_density_matrix(self, spin="ab", index=None):
        """Compute the density matrix array for the specified spin orbitals.

        Parameters
        ----------
        spin : str, optional
           The type of occupied spin orbitals. Options are "a" (for alpha), "b" (for beta), and
           "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        return self.mo.compute_dm(spin, index=index)

    def compute_molecular_orbital(self, points, spin='ab', index=None):
        self._check_arguments(points)
        if self.mo:
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_orbitals(points, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def compute_density(self, points, spin='ab', index=None):
        dm = self.compute_density_matrix(spin=spin, index=index)
        if self.mo:
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_density(dm, points, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def compute_gradient(self, points, spin='ab', index=None):
        r"""Return gradient of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        self._check_arguments(points)
        dm = self.compute_density_matrix(spin=spin, index=index)
        return self._ao.compute_gradient(dm, points)

    def compute_hessian(self, points, spin='ab', index=None):
        r"""Return hessian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        self._check_arguments(points)
        dm = self.compute_density_matrix(spin=spin, index=index)
        return self._ao.compute_hessian(dm, points)

    def compute_laplacian(self, points, spin='ab', index=None):
        r"""Return Laplacian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        hess = self.compute_hessian(points, spin=spin, index=index)
        return np.trace(hess, axis1=1, axis2=2)

    def compute_esp(self, points, spin='ab', index=None, charges=None):
        r"""Return molecular electrostatic potential.

        The molecular electrostatic potential at point :math:`\mathbf{r}` is caused by the
        electron density :math:`\rho` of the specified spin orbitals and set of point charges
        :math:`\{q_A\}_{A=1}^{N_\text{atoms}}` placed at the position of the nuclei. i.e,

        .. math::
           V \left(\mathbf{r}\right) =
             \sum_{A=1}^{N_\text{atoms}} \frac{q_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
             \int \frac{\rho \left(\mathbf{r}"\right)}{\rvert \mathbf{r}" - \mathbf{r} \lvert}
                  d\mathbf{r}"

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        charges : np.ndarray, optional
           Array with shape (n,) representing the point charges at the position of the nuclei.
           When ``None``, the pseudo numbers are used.
        """
        self._check_arguments(points)

        if charges is None:
            charges = self.pseudo_numbers
        elif not isinstance(charges, np.ndarray) or charges.shape != self.numbers.shape:
            raise ValueError("Argument charges should be a 1d-array "
                             "with {0} shape.".format(self.numbers.shape))

        if self.mo:
            dm = self.compute_density_matrix(spin=spin, index=index)
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_esp(dm, points, self.coordinates,
                                        charges, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def compute_ked(self, points, spin='ab', index=None):
        r"""Return positive definite or Lagrangian kinetic energy density.

        .. math::
           \tau_\text{PD} \left(\mathbf{r}\right) =
           \tfrac{1}{2} \sum_i^N n_i \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        self._check_arguments(points)
        if self.mo:
            dm = self.compute_density_matrix(spin=spin, index=index)
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_ked(dm, points, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def _check_arguments(self, points):
        """Check given arguments.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2D-array of floats!")
        if self._ao is None:
            raise AttributeError("Atomic Orbitals information is needed!")
        if self._mo is None:
            raise AttributeError("Molecular Orbitals information is needed!")

class AtomicOrbitals:
    """Gaussian Basis Function or Atomic Orbitals"""

    def __init__(self, basis, coord_type):
        self._basis = basis
        self._coord_type = coord_type

    @classmethod
    def from_molecule(cls, mol):
        """Initialize class given an instance of `Molecule`.

        Parameters
        ----------
        mol : Molecule
            An instance of `Molecule` class.

        """
        basis, coord_type = from_iodata(mol)
        return cls(basis, coord_type)

    @classmethod
    def from_file(cls, fname):
        """Initialize class given a file.

        Parameters
        ----------
        fname : str
            Path to molecule"s files.

        """
        return cls.from_molecule(Molecule.from_file(str(fname)))

    @property
    def center_index(self):
        # FIXME: following code is pretty hacky. it will be used for the orbital partitioning code
        # GOBasis object stores basis set to atom and angular momentum mapping
        # by shell and not by contraction. So it needs to be converted
        try:
            ind_shell_atom = np.array(self._basis.shell_map)
            ind_shell_orbtype = np.array(self._basis.shell_types)

            def num_contr(orbtype):
                """Return the number of contractions for the given orbital type.

                Parameters
                ----------
                orbtype : int
                    Horton's orbital type scheme.
                    Positive number corresponds to cartesian type and negative number corresponds to
                    spherical types.

                Returns
                -------
                num_contr : int
                    Number of contractions for the given orbital type.

                """
                if orbtype < 0:
                    return 1 + 2 * abs(orbtype)
                return 1 + orbtype + sum(i * (i + 1) / 2 for i in range(1, orbtype + 1))

            numcontr_per_shell = np.array([num_contr(i) for i in ind_shell_orbtype])
            # Duplicate each shell information by the number of contractions in shell
            ind_basis_center = np.repeat(ind_shell_atom, numcontr_per_shell)
            # self._ind_basis_orbtype = np.repeat(ind_shell_orbtype, numcontr_per_shell)
        except AttributeError:
            pass
        return ind_basis_center

    def compute_overlap(self):
        r"""Return overlap matrix :math:`\mathbf{S}` of atomic orbitals.

        .. math:: [\mathbf{S}]_ij = int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}

        """
        # check if the coordinate types are mixed
        check_coords = all(elmt == self._coord_type[0] for elmt in self._coord_type)
        if check_coords:
            return overlap_integral(self._basis, coord_type=self._coord_type)
        else:
            raise NotImplementedError("Mixed coordinate types are not supported yet.")


    def compute_orbitals(self, points, transform=None):
        """ Evaluate Orbitals

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        transform : np.ndarray(K_orbs, K_cont)
            Transformation matrix from contractions in the given coordinate system (e.g. AO) to
            linear combinations of contractions (e.g. MO).
            Transformation is applied to the left.
            Rows correspond to the linear combinationes (i.e. MO) and the columns correspond to the
            contractions (i.e. AO).
        """
        return evaluate_basis(self._basis, points, transform=transform,
                              coord_type=self._coord_type)

    def compute_density(self, dm, points, transform=None):
        """Return electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_density(dm, self._basis, points, transform=transform,
                                coord_type=self._coord_type)

    def compute_gradient(self, dm, points):
        """Return gradient of the electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_density_gradient(dm, self._basis, points,
                                         coord_type=self._coord_type)

    def compute_hessian(self, dm, points):
        """Return hessian of the electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_density_hessian(dm, self._basis, points,
                                        coord_type=self._coord_type)

    def compute_esp(self, dm, points, coordinates, charges, transform=None):
        """Return electrostatic potential evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """

        check_coords = all(elmt == self._coord_type[0] for elmt in self._coord_type)
        if check_coords:
            return electrostatic_potential(self._basis, dm, points,
                                           coordinates, charges, transform=transform,
                                           coord_type=self._coord_type)
        else:
            raise NotImplementedError("Mixed coordinate types are not supported yet.")

    def compute_ked(self, dm, points, transform=None):
        """Return positive definite kinetic energy density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_posdef_kinetic_energy_density(dm, self._basis, points,
                                                      transform=transform,
                                                      coord_type=self._coord_type)

class MolecularOrbitals:
    """Molecular orbital class. """

    def __init__(self, occs_a, occs_b, energy_a, energy_b, coeffs_a, coeffs_b):
        self._occs_a, self._occs_b = occs_a, occs_b
        self._energy_a, self._energy_b = energy_a, energy_b
        self._coeffs_a, self._coeffs_b = coeffs_a, coeffs_b

    @classmethod
    def from_molecule(cls, mol):
        """Initialize class given an instance of `Molecule`.

        Parameters
        ----------
        mol : Molecule
            An instance of `Molecule` class.

        """
        if hasattr(mol, "mo") and mol.mo is not None:
            occs_a, occs_b = mol.mo.occsa, mol.mo.occsb
            energy_a, energy_b = mol.mo.energiesa, mol.mo.energiesb
            coeffs_a, coeffs_b = mol.mo.coeffsa, mol.mo.coeffsb
            return cls(occs_a, occs_b, energy_a, energy_b, coeffs_a, coeffs_b)
        else:
            raise AttributeError('Passed mol is missing MO information')

    @classmethod
    def from_file(cls, fname):
        """initialize class given a file.

        parameters
        ----------
        fname : str
            path to molecule"s files.

        """
        return cls.from_molecule(Molecule.from_file(fname))

    @property
    def homo_index(self):
        """Index of alpha and beta HOMO orbital."""
        index_a = np.argwhere(self._occs_a == 0.)[0, 0]
        index_b = np.argwhere(self._occs_b == 0.)[0, 0]
        return index_a, index_b

    @property
    def lumo_index(self):
        """Index of alpha and beta LUMO orbital."""
        return self.homo_index[0] + 1, self.homo_index[1] + 1

    @property
    def homo_energy(self):
        """Energy of alpha and beta HOMO orbital."""
        return self._energy_a[self.homo_index[0] - 1], self._energy_b[self.homo_index[1] - 1]

    @property
    def lumo_energy(self):
        """Energy of alpha and beta LUMO orbital."""
        return self._energy_a[self.lumo_index[0] - 1], self._energy_b[self.lumo_index[1] - 1]

    @property
    def occupation(self):
        """Orbital occupation of alpha and beta electrons."""
        return self._occs_a, self._occs_b

    @property
    def energy(self):
        """Orbital energy of alpha and beta electrons."""
        return self._energy_a, self._energy_b

    @property
    def coefficient(self):
        """Orbital coefficient of alpha and beta electrons.

        The alpha and beta orbital coefficients are each storied in a 2d-array in which
        the columns represent the basis coefficients of each molecular orbital.
        """
        return self._coeffs_a, self._coeffs_b

    @property
    def nelectrons(self):
        """Number of alpha and beta electrons."""
        return np.sum(self._occs_a), np.sum(self._occs_b)

    def compute_dm(self, spin='ab', index=None):
        available_spins = ['a', 'b', 'ab']

        if spin == 'a':
            arr = np.dot(self._coeffs_a * self._occs_a, self._coeffs_a.T)
        elif spin == 'b':
            arr = np.dot(self._coeffs_b * self._occs_b, self._coeffs_b.T)
        elif spin == 'ab':
            arr = np.dot(self._coeffs_a * self._occs_a, self._coeffs_a.T) + \
                  np.dot(self._coeffs_b * self._occs_b, self._coeffs_b.T)
        else:
            raise ValueError(f'Spin must be one of the following: {available_spins}')

        if index is not None:
            # convert to numpy array
            index = np.asarray(index)
            # check
            if index.ndim == 0:
                index = index.reshape(1)
            if index.ndim >= 2:
                raise ValueError("Indices should be given as a one-dimensional numpy array.")
            index -= 1
            if np.any(index < 0):
                raise ValueError(
                    "Indices cannot be less than 1. Note that indices start from 1."
                )
            arr = arr[index[:, np.newaxis], index[np.newaxis, :]]
        return arr

    def compute_overlap(self):
        return NotImplementedError
