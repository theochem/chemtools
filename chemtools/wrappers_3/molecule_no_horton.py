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

        self._iodata = iodata
        self._coords = self._iodata.atcoords
        self._numbers = self._iodata.atnums

        if hasattr(self._iodata, "obasis"):
            self._ao = AtomicOrbitals.from_molecule(self)
        else:
            self._ao = None

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

    # TODO: Check whether mol has mol.one_rdms['scf']?
    # Add docs where needed

    # Include basis as an argument for computing properties?

    # Include new functionality (e.g. Laplassian of density)?

    def __init__(self, basis, coord_type):
        self._basis = basis
        self._coord_type = coord_type

    @classmethod
    def from_molecule(cls, mol):
        basis, coord_type = from_iodata(mol)
        return cls(basis, coord_type)

    @classmethod
    def from_file(cls, fname):
        return cls.from_molecule(Molecule.from_file(str(fname)))

    def compute_overlap(self):
        r"""Return overlap matrix :math:`\mathbf{S}` of atomic orbitals.

        .. math:: [\mathbf{S}]_ij = int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}

        """
        return overlap_integral(self._basis, coord_type="spherical")

    def compute_orbitals(self, points):

        return evaluate_basis(self._basis, points, coord_type="spherical")

    def compute_density(self, dm, points):
        """Return electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return evaluate_density(dm, self._basis, points,
                                coord_type="spherical")

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
                                         coord_type="spherical")

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
                                        coord_type="spherical")

    def compute_esp(self, dm, points, coordinates, charges):
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
            dm, self._basis, points, coord_type="spherical")

class MolecularOrbitals:
    pass
