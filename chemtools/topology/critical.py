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
from scipy.spatial import cKDTree, KDTree
from scipy.spatial.distance import cdist

from chemtools.topology.point import CriticalPoint
from chemtools.topology.ode import bond_paths_rk45


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
        # dictionary for storing bond paths of each bond-critical point using the index
        #  of critical point as key and the items are the list: [ndarray(M_1, 3), ndarray(M_2, 3)]
        #  signifying the two different bond paths.
        self._bp = {}

    @property
    def bp(self):
        """Dictionary of bond paths for each bond-critical point. Keys are indices of `bcp`."""
        return self._bp

    @property
    def cps(self):
        """Sequence of CriticalPoint instances."""
        return [cp for cp_class in list(self._cps.values()) for cp in cp_class]

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
                    coord = self._root_vector_func(point.copy())
                except np.linalg.LinAlgError as _:
                    continue
                # add critical point if it is new
                if not np.any([np.linalg.norm(coord - cp.coordinate) < 1.e-3 for cp in self.cps]):
                    dens = self.func(coord)
                    grad = self.grad(coord)
                    # skip critical point if its dens & grad are zero
                    if abs(dens) < 1.e-4 and np.all(abs(grad) < 1.e-4):
                        continue
                    # compute rank & signature of critical point
                    eigenvals, eigenvecs = np.linalg.eigh(self.hess(coord))
                    cp = CriticalPoint(coord, eigenvals, eigenvecs, 1e-4)
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

    def newton_step(self, pts, alpha=1.0, use_log=False):
        def grad_func_log(dens_vals, grad_vals):
            return grad_vals / dens_vals[:, None]

        def hess_func_log(dens_vals, grad_vals, hess_vals):
            grad_mat = np.einsum("ni,nj->nij", grad_vals, grad_vals)
            return hess_vals / dens_vals[:, None, None] - grad_mat / dens_vals[:, None, None] ** 2.0

        dens = self.func(pts)
        grad = self.grad(pts)
        hess = self.hess(pts)
        if use_log:
            grad_log = grad_func_log(dens, grad)
            hess_log = hess_func_log(dens, grad, hess)
            solution = np.linalg.solve(hess_log, grad_log)
        else:
            solution = np.linalg.solve(hess, grad)

        return pts - alpha[:, None] * solution, np.linalg.norm(grad, axis=1), hess

    def remove_duplicates(self, pts_converged, decimal_uniq, rad_uniq):
        r"""
        Remove duplicates found in `pts_converged`.

        The procedure is as follows:
        1. Points are rounded to `decimal_uniq` decimal place and grouped together. If
           the distance between points is still less than `rad_uniq`, then step 2 is performed.
        2. Using a KDTree, each point a ball of radius `rad_uniq` is constructed and
            the point within the ball with the smallest gradient norm is selected as
            the "unique" point. This process is iterated until all points are within a ball.

        """
        # Round to `decimal_uniq` decimal place and remove any non-unique points.
        test, indices = np.unique(
            np.round(pts_converged, decimals=decimal_uniq),
            axis=0,
            return_index=True
        )
        pts_converged = test[np.argsort(indices)]

        # Get final gradients
        final_gradients = np.linalg.norm(self.grad(pts_converged), axis=1)

        # See if any duplicates were found then run the kdtree algorithm
        dist = cdist(pts_converged, pts_converged)
        dist[np.diag_indices(len(pts_converged))] = 1.0
        maxim = np.any(dist < rad_uniq, axis=1)
        if np.any(maxim):
            indices = set(range(len(pts_converged)))
            tree = KDTree(pts_converged)
            indices_include = []
            while len(indices) != 0:
                index = indices.pop()
                query = tree.query_ball_point(pts_converged[index], r=rad_uniq)
                gradient_ball = final_gradients[query]
                sort_indices = np.argsort(gradient_ball)
                indices_include.append(query[sort_indices[0]])  # Pick the point with smallest gradient
                indices = indices.difference(query)

            final_gradients = final_gradients[indices_include]
            pts_converged = pts_converged[indices_include]

        return pts_converged, final_gradients

    def find_critical_points_vectorized(
        self,
        atom_coords,
        xtol=1e-8,
        ftol=1e-8,
        ss_rad=0.1,
        # u_bnd_rad=1.5,
        init_norm_cutoff=0.5,
        dens_cutoff=1e-5,
        decimal_uniq=12,
        rad_uniq=0.01,
        maxiter=250,
        failure_use_scipy=True,
        use_log=False,
        verbose=True
    ):
        r"""
        Find and store the critical points with a vectorization approach.

        Parameters
        ----------
        atom_coords: ndarray(M, 3)
            Atomic coordinates, used to construct atomic grids over each center as
            initial points.
        xtol: float, optional
            Point converges if the difference between two iterations is less than `xtol`
            in all three dimensions. Also needs to satisfy `ftol` convergence.
        ftol: float, optional
            Point converges if the norm of the gradient is less than `ftol`. Also needs
            to satisfy `xtol` convergence.
        ss_rad: float, optional
            Points that are 2`ss_rad` away from the atomic coordinates are included as
            initial guesses.
        init_norm_cutoff: float, optional
            Points whose gradient norm is less than `init_norm_cutoff` are excluded when constructing
            the initial atomic grid.
        dens_cutoff: float, optional
            If the point diverges during the critical point finding process then it excludes it based
            on its density value being less than `dens_cutoff`.
        decimal_uniq: float, optional
            Rounds each critical point to this `decimal_uniq` decimal place. Used to condense points
            together that are close together.
        rad_uniq: float, optional
            Points that are `rad_uniq` apart, are grouped together and the one with the smallest
            gradient norm is picked. Used to condense points together that are close together
        use_log: float, optional
            If true, then uses the logarithm of the electorn density. Speeds-up convergence.
        verbose: float, optional
            If true, then prints the sum of the norm of the gradient per iteration and number of
            points that didn't converge.

        Warns
        -----
        RuntimeWarning
            - Poincare-Hopf equation isn't satisfied. This implies not all critical points were
              found.
            - The norm of the gradient of the final crtical points are not less than 1e-5.

        Notes
        -----
        - Newton-Ralphson equation is used to find the roots of the gradient.

        - The algorithm is as follows:
            1. Constructs an atomic grid over each atomic center.
            2. Points that are close to center and that have small gradients are
                included to be initial guesses. Points whose density values is small are not.
            3. While-loop is done until each point from the initial guess converged.
            4. Points whose density values are too small are removed from the iteration process.
            5. If points didn't converge then Scipy is used instead but points are reduced
               considerably since Scipy is slower as it doesn't accept vectorization.
            6. Points that are identical are merged together by first rounding to some decimal
               places and then a KDTree is used to group the points into small balls and the
               points with the smallest gradient are picked out of each ball.

        """
        # Construct set of points around each center
        # from grid.onedgrid import OneDGrid
        # from grid.atomgrid import AtomGrid
        # natoms = len(atom_coords)
        # atom_grids_list = []
        # rad_grid = np.arange(0.0, u_bnd_rad, ss_rad)
        # rgrid = OneDGrid(points=rad_grid, weights=np.empty(rad_grid.shape))
        # for i in range(natoms):
        #     atomic_grid = AtomGrid.from_preset(10, preset="fine", rgrid=rgrid, center=atom_coords[i])
        #     atom_grids_list.append(atomic_grid.points)
        # pts = np.vstack(atom_grids_list)

        # Use density to remove points greater than certain value (0.001 au)
        cut_off = self.func(self._points) > dens_cutoff
        pts = self._points[cut_off, :]

        # Include points that have small gradients and close to the centers
        close_centers = np.any(cdist(pts, atom_coords) < 2 * ss_rad, axis=1)
        norm_grad = np.linalg.norm(self.grad(pts), axis=1)
        pts_decis = close_centers | (norm_grad < init_norm_cutoff)
        pts = pts[pts_decis]
        norm_grad = norm_grad[pts_decis]

        if verbose:
            print(f"Start Critical Points Finding.")
            print(f"Number of initial points: {len(pts)}")

        # Solve critical-point finding via Newton-Ralpson with backtracing
        # Keep going until all points converge
        p1, p0, grad_norm0 = np.ones(pts.shape) * 1000.0, pts, norm_grad
        alpha = np.ones(len(p1))  # Initial Step-Sizes
        converged_pts = np.empty((0, 3), dtype=float)
        niter = 0

        while len(p1) != 0 and niter < maxiter:
            # Take Newton-Ralphson Step
            p1, grad_norm1, hess = self.newton_step(p0, alpha=alpha, use_log=use_log)

            if verbose:
                print(f"Iteration: {niter + 1},   "
                      f"Sum Gradient Norm of pts left: {(np.sum(grad_norm1))}   "
                      f"Number of pts left: {len(p1)}")
            # TODO: Add Armijo condition

            # Remove points whose density is very small and gradient
            dens1 = self.func(p1)
            include = (dens1 > dens_cutoff)
            p1 = p1[include, :]
            p0 = p0[include, :]
            alpha = alpha[include]
            grad_norm1 = grad_norm1[include]

            # Store points that converge in its 3 dimensions
            err = np.abs(p1 - p0)
            did_converge = np.all(err < xtol, axis=1) & (grad_norm1 < ftol)
            converged_pts = np.vstack((converged_pts, p1[did_converge, :]))

            # Remove points that did converge
            didnt_converge = np.bitwise_not(did_converge)
            p1 = p1[didnt_converge, :]
            alpha = alpha[didnt_converge]
            grad_norm1 = grad_norm1[didnt_converge]

            # Update variables for next iteration
            p0 = p1
            niter += 1

        if niter == maxiter:
            if failure_use_scipy:
                p1, _ = self.remove_duplicates(p1, decimal_uniq, rad_uniq)
                if verbose:
                    print(f"Convergence not achieved, attempt using Scipy."
                          f" Number of points didn't converge: {len(p1)}")

                # Try Scipy instead: Go through each pt individually.
                #  Note this is slower but convergence is guaranteed.
                other_converged_pts = []
                for x in p1:
                    sol = root(lambda pt: self.grad(np.array([pt]))[0], x0=x,
                               # jac=lambda pt: self.hess(np.array([pt]))[0],
                               method="df-sane")
                    assert sol.success, f"One of the roots did not converge using scipy {sol}"
                    other_converged_pts.append(sol.x)

                converged_pts = np.vstack((converged_pts, np.array(other_converged_pts, dtype=float)))
            else:
                raise RuntimeError(f"Maximum number of iterations {niter} == {maxiter} was reached.")

        # Remove duplicates from entire set of points.
        converged_pts, final_gradients = self.remove_duplicates(converged_pts, decimal_uniq, rad_uniq)

        if verbose:
            print(f"Final number of points found {len(converged_pts)}")

        # Check the gradient are all small
        if np.any(final_gradients > 1e-5):
            warnings.warn("Gradient of the final critical points are not all smaller than 1e-5.",
                          RuntimeWarning)

        # compute rank & signature of critical point
        hess = self.hess(converged_pts)
        eigs, evecs = np.linalg.eigh(hess)
        for i in range(len(converged_pts)):
            cp = CriticalPoint(converged_pts[i], eigs[i], evecs[i], 1e-4)
            self._cps.setdefault((cp.rank[0], cp.signature[0]), []).append(cp)

        # check Poincare–Hopf equation
        if not self.poincare_hopf_equation:
            warnings.warn("Poincare–Hopf equation is not satisfied.", RuntimeWarning)

    def find_bond_paths(self, eps_dir=1e-4, tol=1e-8, ss_0=1e-8, max_ss=0.25, maxiter=2000):
        r"""
        Find bond paths of each bond critical point.

        The critical points must be found first.
        RungeRunge-Kutta of order 4(5) with adaptive step-size is used to solve for
        bond path.

        Parameters
        ----------
        eps_dir : float, optional
            The displacement in the direction in the eigenvector of the Hessian at the bond
            critical point.
        tol : float, optional
            Tolerance for the adaptive step-size.
        ss_0 : float, optional
            The initial step-size of the ODE (RK45) solver.
        max_ss : float, optional
            Maximum step-size of the ODE (RK45) solver.
        maxiter : int, optional
            The maximum number of iterations of taking a step in the ODE solver.

        """
        if not self._cps:
            raise RuntimeError(f"Dictionary empty, critical points should be found first.")
        if eps_dir < 0.0:
            raise ValueError(f"eps_dir {eps_dir} must be positive.")

        all_pts = []
        for i, cp in enumerate(self.bcp):
            coordinate = cp.coordinate
            eigs = cp.eigenvalues
            evecs = cp._eigenvectors
            if np.any(eigs[0] < -1e-8):
                bond_path_dir = evecs[:, eigs[0] > 0.0].T

                for b in bond_path_dir:
                    dir1 = coordinate + eps_dir * b
                    dir2 = coordinate - eps_dir * b
                    all_pts.append(dir1)
                    all_pts.append(dir2)
        all_pts = np.array(all_pts)
        bond_paths = bond_paths_rk45(
            all_pts,
            self.func,
            self.grad,
            ss_0=ss_0,
            max_ss=max_ss,
            tol=tol,
            maxiter=maxiter
        )

        count = 0
        for i_cp in range(len(self.bcp)):
            self._bp[i_cp] = [np.array(bond_paths[count]), np.array(bond_paths[count + 1])]
            count += 2
