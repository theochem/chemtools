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
# pragma pylint: disable=wrong-import-position
"""Module for topological analysis of Quantum Chemistry Output Files."""


import numpy as np

from chemtools.wrappers.molecule import Molecule
from chemtools.utils.cube import UniformGrid
from chemtools.topology.critical import Topology


class TopologicalTool(Topology):
    """Topological Analysis of Electron Density."""

    def __init__(self, func_dens, func_grad, func_hess, points, coordinates=None):
        """Initialize class.

        Parameters
        ----------
        func_dens : callable
            The function to compute the scalar density value at a given point.
        func_grad : callable
            The function to compute the gradient of density at a given point.
        func_hess : callale
            The function to compute the hessian of density at a given point.
        points : np.ndarray
            Cartesian coordinates of points used for critical point search.
        coordinates : np.ndarray, optional
            Cartesian coordinates of atoms. If given, they are added for critical point search.

        """
        if points.ndim != 2 and points.shape[1] != 3:
            raise ValueError("Argument point should be a 2D-array with 3 columns.")
        if coordinates.ndim != 2 and coordinates.shape[1] != 3:
            raise ValueError("Argument coordinates should be a 2D-array with 3 columns.")
        super(TopologicalTool, self).__init__(func_dens, func_grad, func_hess, points, coordinates)
        self.find_critical_points()

    @classmethod
    def from_molecule(cls, molecule, spin="ab", index=None, points=None):
        """Initialize class from ``Molecule`` object.

        Parameters
        ----------
        molecule : instance of `Molecule` class.
            Instance of `Molecular` class.
        spin : str, optional
            The type of occupied spin orbitals; options are "a", "b" & "ab".
        index : int or Sequence of int, optional
            Sequence of integers representing the index of spin orbitals.
            If None, all occupied spin orbitals are included.
        points : np.ndarray, optional
            Cartesian coordinates of points used for critical point search.
            If None, cubic grid points with spacing=0.1 & extension=0.1 are used.

        """
        if points is None:
            points = UniformGrid.from_molecule(molecule, spacing=0.1, extension=0.1).points
        func_dens = cls._wrapper_compute_density(molecule, spin=spin, index=index)
        func_grad = cls._wrapper_compute_gradient(molecule, spin=spin, index=index)
        func_hess = cls._wrapper_compute_hessian(molecule, spin=spin, index=index)
        return cls(func_dens, func_grad, func_hess, points, molecule.coordinates)

    @classmethod
    def from_file(cls, fname, spin="ab", index=None, points=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        fname : str
            A string representing the path to a molecule's fname.
        spin : str, optional
            The type of occupied spin orbitals; options are "a", "b" & "ab".
        index : int or Sequence of int, optional
            Sequence of integers representing the index of spin orbitals.
            If None, all occupied spin orbitals are included.
        points : np.ndarray, optional
            Cartesian coordinates of points used for critical point search.
            If None, cubic grid points with spacing=0.1 & extension=0.1 are used.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, spin=spin, index=index, points=points)

    @staticmethod
    def _wrapper_compute_density(molecule, spin, index):
        def compute_density(point):
            if point.ndim == 1:
                point = point[np.newaxis, :]
            dens = molecule.compute_density(point, spin, index)
            return dens
        return compute_density

    @staticmethod
    def _wrapper_compute_gradient(molecule, spin, index):
        def compute_gradient(point):
            revert = False
            if point.ndim == 1:
                point = point[np.newaxis, :]
                revert = True
            grad = molecule.compute_gradient(point, spin, index)
            if revert:
                grad = grad.flatten()
            return grad
        return compute_gradient

    @staticmethod
    def _wrapper_compute_hessian(molecule, spin, index):
        def compute_hessian(point):
            if point.ndim == 1:
                point = point[np.newaxis, :]
            hess = molecule.compute_hessian(point, spin, index)[0]
            hess += hess.T
            hess[np.diag_indices(3)] /= 2.
            return hess
        return compute_hessian
