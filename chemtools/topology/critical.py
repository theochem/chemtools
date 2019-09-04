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
r"""Functionality for finding critical points of any scalar function."""


import warnings
import numpy as np

from scipy.spatial import cKDTree

from chemtools.topology.point import CriticalPoint


class Topology(object):
    """Topology class for searching critical points given scalar function."""

    def __init__(self, func, func_grad, func_hess, points, coords=None, n_neighbours=4):
        r"""Initialize Topology class instance.

        Parameters
        ----------
        func : callable[np.ndarray(N, 3) -> float]
            Method for computing the scalar function.
        func_grad : callable[np.ndarray(N, 3) -> np.ndarray(3,)]
            Method for computing the gradient vector of scalar function.
        func_hess : callable[np.ndarray(N, 3) -> np.ndarray(3, 3)]
            Method for computing the hessian matrix of scalar function.
        points : np.ndarray(M, 3)
            Cartesian coordinates of :math:`M` initial guess points.
        coords : np.ndarray(N, 3), optional
            Cartesian coordinates of :math:`N` atomic centers to use as additional initial guesses.
        n_neighbours: int, optional
            Number of polyhedron vertices.

        """
        if points.ndim != 2 and points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D-array with 3 columns!")
        if points.shape[0] < 4:
            raise ValueError("At least 4 points are needed for critical point search!")
        self._points = points

        if coords is not None and coords.ndim != 2 and coords.shape[1] != 3:
            raise ValueError("Argument coords should be a 2D-array with 3 columns.")
        if coords is not None:
            self._coords = coords
            self._kdtree = cKDTree(np.vstack((self._coords, self._points)))
        else:
            self._coords = None
            self._kdtree = cKDTree(points)

        self.func = func
        self.grad = func_grad
        self.hess = func_hess
        # pre-compute coordinates of polyhedron vertices
        self._neighbours = self._polyhedron_coordinates(n_neighbours)
        # dictionary for storing critical points using (rank, signature) as key
        self._cps = {}

    @property
    def cps(self):
        """Sequence of CriticalPoint instances."""
        return [cp for cp_class in self._cps.values() for cp in cp_class]

    @property
    def nna(self):
        """Sequence of CriticalPoint instances representing (non-)nuclear attractors."""
        return self._cps.get((3, -3), [])

    @property
    def bcp(self):
        """Sequence of CriticalPoint instances representing bond critical points."""
        return self._cps.get((3, -1), [])

    @property
    def rcp(self):
        """Sequence of CriticalPoint instances representing ring critical points."""
        return self._cps.get((3, 1), [])

    @property
    def ccp(self):
        """Sequence of CriticalPoint instances representing cage critical points."""
        return self._cps.get((3, 3), [])

    @property
    def poincare_hopf_equation(self):
        """bool: whether the Poincare–Hopf equation is satisfied."""
        return len(self.nna) - len(self.bcp) + len(self.rcp) - len(self.ccp) == 1

    def add_points(self, points):
        """Add extra points as seeds to search for critical points.

        Parameters
        ----------
        points : np.ndarray
            Cartesian coordinates of points.

        """
        self._kdtree = cKDTree(np.vstack((self._kdtree.data, points)))

    def find_critical_points(self):
        """Find and store the critical points."""
        neighs = np.zeros((4 * len(self._kdtree.data), 3))
        for index, point in enumerate(self._kdtree.data):
            # compute distance to 4 closest grid points
            dists, _ = self._kdtree.query(point, 4)
            # store coordinates of neighbouring polyhedron vertices surrounding the point
            neighs[4 * index: 4 * index + 4, :] = point + np.max(dists) * self._neighbours

        # compute the gradient norm of points & surrounding vertices
        points_norm = np.linalg.norm(self.grad(self._kdtree.data), axis=-1)
        neighs_norm = np.linalg.norm(self.grad(neighs), axis=-1)

        for index, point in enumerate(self._kdtree.data):
            point_norm = points_norm[index]
            neigh_norm = neighs_norm[4 * index: 4 * index + 4]
            # use central point as initial guess for critical point finding
            if index < len(self._coords) or np.all(point_norm < neigh_norm):
                try:
                    cp_point = self._root_vector_func(point.copy())
                except Exception as _:
                    continue
                # add critical point if it is new
                if not np.any([np.linalg.norm(cp_point - cp.coordinate) < 1.e-3 for cp in self.cps]):
                    dens = self.func(cp_point)
                    grad = self.grad(cp_point)
                    # skip critical point if its dens & grad are zero
                    if abs(dens) < 1.e-4 and np.all(abs(grad) < 1.e-4):
                        continue
                    # compute rank & signature of critical point
                    eigenvals, eigenvecs = np.linalg.eigh(self.hess(cp_point))
                    cp = CriticalPoint(cp_point, eigenvals, eigenvecs, 1e-4)
                    self._cps.setdefault((cp.rank[0], cp.signature[0]), []).append(cp)
        # check Poincare–Hopf equation
        if not self.poincare_hopf_equation:
            warnings.warn("Poincare–Hopf equation is not satisfied.", RuntimeWarning)

    def _root_vector_func(self, guess, maxiter=5000):
        """Find root of a multivariate function using Newton-Raphson method.

        Parameters
        ----------
        guess : np.ndarray(3,)
            Cartesian coordinates of initial guess.
        maxiter: int, optional
            Maximum number of iterations.

        """
        niter = 0
        norm = np.inf
        while niter < maxiter and norm > 1.0e-6:
            grad = self.grad(guess)
            hess = self.hess(guess)
            guess = guess - np.dot(np.linalg.inv(hess), grad[:, np.newaxis]).flatten()
            norm = np.linalg.norm(grad, axis=-1)
            niter += 1
        return guess

    @staticmethod
    def _polyhedron_coordinates(n_vertices):
        """Compute Cartesian coordinates of polyhedron vertices.

        Parameters
        ----------
        n_vertices : int
            Number of polyhedron vertices.

        """
        if n_vertices == 4:
            coords = np.array([[np.sqrt(8.0) / 3.0, 0.0, -1.0 / 3.0],
                               [-np.sqrt(2.0) / 3.0, np.sqrt(2.0 / 3.0), -1.0 / 3.0],
                               [-np.sqrt(2.0) / 3.0, -np.sqrt(2.0 / 3.0), -1.0 / 3.0],
                               [0.0, 0.0, 1.0]])
            coords *= np.sqrt(24)
        else:
            raise NotImplementedError("Number of vertices {} is not supported".format(n_vertices))
        return coords
