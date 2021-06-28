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
        self._dm = self._iodata.one_rdms['scf']

        if hasattr(self._iodata, "obasis"):
            self._ao = AtomicOrbitals.from_molecule(self._iodata)
        else:
            self._ao = None

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
    def ao(self):
        """Atomic orbital instance"""
        return self._ao

    @property
    def mo(self):
        """Molecular orbital instance."""
        return self._mo

    def compute_molecular_orbital(self, points):
        self._check_arguments(points)
        if self.mo:
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_orbitals(points, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def compute_density(self, points):
        if self.mo:
            tf = (self._iodata.mo.coeffs * self._iodata.mo.occs).dot(self._iodata.mo.coeffs.T)
            return self._ao.compute_density(self._dm, points, transform=tf)
        else:
            raise ValueError("mo attribute cannot be empty or None")

    def compute_gradient(self, points):
        r"""Return gradient of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        self._check_arguments(points)
        return self._ao.compute_gradient(self._dm, points)

    def compute_hessian(self, points):
        r"""Return hessian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        self._check_arguments(points)
        return self._ao.compute_hessian(self._dm, points)

    def compute_laplacian(self, points):
        r"""Return Laplacian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        hess = self.compute_hessian(points)
        return np.trace(hess, axis1=1, axis2=2)

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

    def compute_overlap(self):
        r"""Return overlap matrix :math:`\mathbf{S}` of atomic orbitals.

        .. math:: [\mathbf{S}]_ij = int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}

        """
        # check if the coordinate types are mixed
        check_coords = all(elmt == self._coord_type[0] for elmt in self._coord_type)
        if check_coords:
            return overlap_integral(self._basis, coord_type=self._coord_type)
        else:
            raise NotImplementedError("GBasis does not support mixed coordinate types yet.")


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
        return electrostatic_potential(self._basis, dm, points,
                                       coordinates, charges)

    def compute_ked(self, dm, points):
        """Return positive definite kinetic energy density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_posdef_kinetic_energy_density(
            dm, self._basis, points, coord_type=self._coord_type)

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

    def compute_overlap(self):
        return NotImplementedError
