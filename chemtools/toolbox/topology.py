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
"""Module for topological analysis of scalar fields."""


import sys
import numpy as np

from chemtools.utils.cube import UniformGrid
from chemtools.topology.critical import Topology
from chemtools.outputs.vmd import print_vmd_script_topology

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule


class TopologicalTool(Topology):
    """Topological analysis of scalar functions."""

    def __init__(self, func_dens, func_grad, func_hess, points, coordinates=None):
        """Initialize class for topological analysis of arbitrary scalar function.

        Parameters
        ----------
        func_dens : callable
            The function to compute the scalar value at a given point.
        func_grad : callable
            The function to compute the gradient of scalar function at a given point.
        func_hess : callable
            The function to compute the hessian of scalar function at a given point.
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
        """Initialize class from `Molecule` object for topological analysis of electron density.

        Parameters
        ----------
        molecule : `Molecule`
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
        """Initialize class from wave-function file for topological analysis of electron density.

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

    def generate_scripts(self, fname, radius=0.2):
        """Generate VMD script to visualize critical points & gradient path.

        Parameters
        ----------
        fname : str
            The name of the VMD script file.
        radius : float, optional
            Radius of spheres representing the critical points.

        """
        cps = {'gray': [nna.coordinate for nna in self.nna],
               'blue': [bcp.coordinate for bcp in self.bcp],
               'green': [rcp.coordinate for rcp in self.rcp],
               'red': [ccp.coordinate for ccp in self.ccp]}
        print_vmd_script_topology(fname, cps, radius=radius)

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
            return hess
        return compute_hessian
