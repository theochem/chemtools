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

from scipy.optimize import root
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
        self, coors, value_func, grad_func, hess_func, points=None, extra=5
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
        if coors.ndim != 2:
            raise ValueError("Input array need to be (N, 3) shape.")
        self.coors = coors
        self.v_f = value_func
        self.g_f = grad_func
        self.h_f = hess_func
        # num of the maximum equals to num of atoms
        if points.ndim != 2 and points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D-array with 3 columns!")
        if points.shape[0] < 4:
            raise ValueError("At least 4 points are needed for critical point search!")
        self._kdtree = KDTree(np.vstack((coors, points)))
        self._found_ct = np.zeros((0, 3))
        self._found_ct_type = np.zeros(0, dtype=int)
        self._nna = []
        self._bcp = []
        self._rcp = []
        self._ccp = []

    def add_points(self, points):
        """Add a point to exiting initial guess points.

        Parameters
        ----------
        points : np.ndarray(3,)
            a numpy 3d array
        """
        new_points = np.vstack((self._kdtree.data, points))
        self._kdtree = KDTree(new_points)

    @staticmethod
    def _construct_cage(point, length, n_points=4):
        """Construct points to encage given guess point.

        Parameters
        ----------
        point : np.ndarray(3,)
            initial guess point
        length : float
            multiplier for the lengh of cage side.
        n_points : int, default to 4
            number of points to form the cage

        Returns
        -------
        np.ndarray(4, 3)
            coordinates of 4 points forming a tetrahedral around
            initial guess point

        Raises
        ------
        NotImplementedError
            n_points only support 4 at this moment
        """
        if n_points == 4:
            constant = np.sqrt(24) * length
            p1 = point + constant * np.array([np.sqrt(8.0) / 3.0, 0.0, -1.0 / 3.0])
            p2 = point - constant * np.array([np.sqrt(2.0) / 3.0, -np.sqrt(2.0 / 3.0), 1.0 / 3.0])
            p3 = point - constant * np.array([np.sqrt(2.0) / 3.0, np.sqrt(2.0 / 3.0), 1.0 / 3.0])
            p4 = point + constant * np.array([0.0, 0.0, 1.0])
            points = np.vstack((p1, p2, p3, p4))
        else:
            raise NotImplementedError(
                "Given args n_point={} is not valid".format(n_points)
            )
        return points

    def find_critical_pts(self):
        """Start the critical point finding main function."""
        for index, init_point in enumerate(self._kdtree.data):
            length, _ = self._kdtree.query(init_point, 4)
            tetrahedral = self._construct_cage(init_point, np.max(length))
            g_values = self.g_f(tetrahedral)
            central_g = self.g_f(init_point)
            good_guess = np.all(np.linalg.norm(central_g) < np.linalg.norm(g_values, axis=-1))
            if index < len(self.coors) or good_guess:
                try:
                    cp_coord = self.converge_to_cp(init_point.copy())
                except Exception as _:
                    continue
                # add if a new CP
                if not np.any([np.linalg.norm(p - cp_coord) < 1.e-3 for p in self._found_ct]):
                    dens = self.v_f(cp_coord[np.newaxis, :])
                    grad = self.compute_grad(cp_coord)
                    # add if dens & grad are not zero
                    if abs(dens) < 1.e-4 and np.all(abs(grad) < 1.e-4):
                        continue
                    cp, sig = self._classify_critical_pt(cp_coord)
                    self._add_critical_point(cp, sig)
        if not self.poincare_hopf_equation:
            warnings.warn("Poincare–Hopf equation is not satisfied.", RuntimeWarning)

    def converge_to_cp(self, guess, maxiter=5000):
        niter = 0
        norm = np.inf
        while niter < maxiter and norm > 1.e-4:
            grad = self.g_f(guess)
            norm = np.linalg.norm(grad, axis=-1)
            hess = self.compute_hess(guess)
            guess = guess - np.dot(np.linalg.inv(hess), grad[:, np.newaxis]).flatten()
            niter += 1
        return guess

    def _add_critical_point(self, ct_pt, ct_type):
        """Add criticla point to instance.

        Parameters
        ----------
        ct_pt : np.ndarray(3,)
            critical point coordinates
        ct_type : int
            sum of +(-) eigenvalues of points
        """
        signature_dict = {
            -3: self._nna,
            3: self._ccp,
            -1: self._bcp,
            1: self._rcp,
        }
        signature_dict[ct_type].append(ct_pt)
        self._found_ct = np.vstack((self._found_ct, ct_pt.point))
        self._found_ct_type = np.append(self._found_ct_type, int(ct_type))

    @property
    def poincare_hopf_equation(self):
        """bool: whether the Poincare–Hopf equation is satisfied."""
        pre_hopf = len(self._nna) - len(self._bcp) + len(self._rcp) - len(self._ccp)
        return pre_hopf == 1

    def _classify_critical_pt(self, point, eigen_cutoff=1e-4):
        """Classify the type of given critical point.

        Parameters
        ----------
        point : np.ndarray(3,)
            Coordinates of given critical point
        eigen_cutoff : float, default to 1e-4
            The engenvalue cutoff incase too small value

        Returns
        -------
        tuple(CriticalPoint, int)
            CriticalPoint instance with all property of crit pt
            and the sum of sign of eigenvalues of given point
        """
        hess_crit = self.h_f(point)
        eigenvals, eigenvecs = np.linalg.eigh(hess_crit)
        cp = CriticalPoint(point, eigenvals, eigenvecs, eigen_cutoff)
        return cp, cp.signature[0]
