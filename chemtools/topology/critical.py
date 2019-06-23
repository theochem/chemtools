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

from scipy.spatial import KDTree

from chemtools.topology.eigenvalues import EigenValueTool


class CriticalPoint(EigenValueTool):
    """Critical Point data class."""

    def __init__(self, point, eigenvalues, eigenvectors, eps=1e-15):
        """Initialize CriticalPoint data class.

        Parameters
        ----------
        point : np.ndarray(N,)
            Coordinates of the critical point
        eigenvalues : np.ndarray
            eigenvalue of its hessian function
        eigenvectors : TYPE
            eigenfunction of its hessian function

        """
        super(CriticalPoint, self).__init__(eigenvalues[np.newaxis, :], eps=eps)
        self._point = point
        self._eigenvectors = eigenvectors

    @property
    def point(self):
        """Cartesian coordinates of critical point."""
        return self._point


class Topology(object):
    """Topology class for searching critical points given scalar function."""

    def __init__(
        self, value_func, grad_func, hess_func, points, coords=None, n_neighbours=4,
    ):
        """Initialize Topology class instance.

        Parameters
        ----------
        coors : np.ndarray(N, 3)
            Cartesian coordinates of all atoms
        value_func : Callable[np.ndarray(N, 3) -> float]
            Objective scalar function in space
        grad_func : Callable[np.ndarray(N, 3) -> np.ndarray(3,)]
            1st order derivative of the scarlar function
        hess_func : Callable[np.ndarray(N, 3) -> np.ndarray(3, 3)]
            2nd order derivative of the scarlar function
        points : np.ndarray(M, 3), optional
            An array of 3 dimension initial guess points. If not given,
            a meshgrid (min(x, y, z) - extra, max(x, y, z) + extra, step=0.2)
            will be generated by default.
        extra : int, optional
            Extra space for generating meshgrid. Used in above situation
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
            self._kdtree = KDTree(np.vstack((self._coords, self._points)))
        else:
            self._coords = []
            self._kdtree = KDTree(points)

        self.v_f = value_func
        self.g_f = grad_func
        self.h_f = hess_func
        # pre-compute coordinates of polyhedron vertices
        self._neighbours = self._polyhedron_coordinates(n_neighbours)

        # dictionary for storing critical points using (rank, signature) as key
        self._cps = {}

    @property
    def cps(self):
        """sequence of CriticalPoint instances."""
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

    def find_critical_points(self):
        """Find and store the critical points."""
        for index, point in enumerate(self._kdtree.data):
            # compute distance to 4 closest grid points
            dists, _ = self._kdtree.query(point, 4)
            # compute coordinates of neighbouring polyhedron vertices surrounding the point
            neigh = point + np.max(dists) * self._neighbours
            # compute the gradient norm of point & surrounding vertices
            point_norm = np.linalg.norm(self.g_f(point))
            neigh_norm = np.linalg.norm(self.g_f(neigh), axis=-1)
            # use central point as initial guess for critical point finding
            if index < len(self._coords) or np.all(point_norm < neigh_norm):
                try:
                    cp_point = self._root_vector_func(point.copy())
                except Exception as _:
                    continue
                # add critical point if it is new
                if not np.any([np.linalg.norm(cp_point - cp.point) < 1.e-3 for cp in self.cps]):
                    dens = self.v_f(cp_point)
                    grad = self.g_f(cp_point)
                    # skip critical point if its dens & grad are zero
                    if abs(dens) < 1.e-4 and np.all(abs(grad) < 1.e-4):
                        continue
                    # compute rank & signature of critical point
                    eigenvals, eigenvecs = np.linalg.eigh(self.h_f(cp_point))
                    cp = CriticalPoint(cp_point, eigenvals, eigenvecs, 1e-4)
                    self._cps.setdefault((cp.rank[0], cp.signature[0]), []).append(cp)
        # check Poincare–Hopf equation
        if not self.poincare_hopf_equation:
            warnings.warn("Poincare–Hopf equation is not satisfied.", RuntimeWarning)

    def _root_vector_func(self, guess, maxiter=5000):
        """Find root of a multivariate function using Newton-Raphson method.

        Parameters
        ----------
        guess : np.ndarray
            Cartesian coordinates of initial guess.
        maxiter: int, optional
            Maximum number of iterations.

        """
        niter = 0
        norm = np.inf
        while niter < maxiter and norm > 1.e-4:
            grad = self.g_f(guess)
            hess = self.h_f(guess)
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
            Number of vertices of polyhedron.

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
