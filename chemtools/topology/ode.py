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
import numpy as np

__all__ = [
    "bond_paths_rk45",
]


def _get_normalized_gradient_func(grad_func):
    r"""Returns the normalized version of the function."""
    def norm_grad_func(x):
        gradient = grad_func(x)
        return gradient / np.linalg.norm(gradient, axis=1)[:, None]
    return norm_grad_func


def _RK45_step(pts, grad_func, step_size, grad0=None):
    r"""
    Runge-Kutta fourth and five-order step for the following ode system:

    .. math::
        \frac{d(r(s))}{ds} = \nabla \rho (r(s)),

    where :math:`\nabla \rho(x)` is the gradient of a function.

    Parameters
    ----------
    pts : ndarray(M, 3)
        Points to take a step in.
    grad_func: callable(ndarray(M,3), ndarray(M,3))
        Callable function that takes in points and obtains the gradient.
    step_size: float
        Stepsize for the step
    grad0: ndarray(M, 3)
        If the gradient is already computed on `pts`, then avoids re-computing it.

    Returns
    -------
    (ndarray(M,3), ndarray(M,3)):
        Returns fourth-order and fifth-order Runge-Kutta step, respectively.

    """
    # Fourth and Fifth-Order Runge-Kutta
    if grad0 is None:
        k1 = step_size * grad_func(pts)
    else:
        k1 = step_size * grad0
    k2 = step_size * grad_func(pts + 0.4 * k1)
    k3 = step_size * grad_func(pts + (3.0 / 32) * k1 + (9.0 / 32.0) * k2)
    k4 = step_size * grad_func(pts + (1932 / 2197) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3)
    k5 = step_size * grad_func(
        pts + (439 / 216) * k1 - 8.0 * k2 + (3680 / 513) * k3 - (845 / 4104) * k4
    )
    k6 = step_size * grad_func(
        pts
        - (8.0 / 27.0) * k1
        + 2.0 * k2
        - (3544 / 2565) * k3
        + (1859 / 4104) * k4
        - (11 / 40) * k5
    )

    # Get the fourth and five-order approximation
    y_four = pts + (25.0 / 216.0) * k1 + (1408 / 2565) * k3 + (2197 / 4101) * k4 - k5 / 5.0
    y_five = (
        pts +
        (16.0 / 135.0) * k1
        + (6656 / 12825) * k3
        + (28561 / 56430) * k4
        - (9.0 / 50.0) * k5
        + (2.0 / 55.0) * k6
    )
    return y_four, y_five


def bond_paths_rk45(
    initial_pts,
    dens_func,
    grad_func,
    ss_0=1e-7,
    tol=1e-7,
    max_ss=0.25,
    tol_conv=1e-10,
    maxiter=2000
):
    r"""
    Solves the following problem ODE using Runge-Kutta of order 4(5) with adaptive step-size

    .. math::
        \frac{d(r(t))}{dt} = \frac{\nabla \rho(r(t))}{|| \rho(r() ||}

    over a set of points.

    Parameters
    ----------
    initial_pts: ndarray(N, 3)
        Initial points to solve for steepest-ascent/backtracing.
    dens_func: callable(ndarray(N,3), ndarray(N))
        The electron density function.
    grad_func: callable(ndarray(N,3), ndarray(N,3))
        The gradient of the electron density.
    beta_spheres: ndarray(M,)
        The beta-sphere/trust-region radius of each atom. These are spheres
        centered at each maxima that reduces the convergence of each point.
    maximas: ndarray(M,3)
        The position of each atoms in the molecule.
    ss_0: float, optional
        The initial step-size of the ODE (RK45) solver.
    max_ss: float, optional
        Maximum step-size of the ODE (RK45) solver.
    tol_conv: float, optional
        Tolerance for the convergence of a point from a ODE step.
    tol: float, optional
        Tolerance for the adaptive step-size.
    maxiter: int, optional
        The maximum number of iterations of taking a step in the ODE solver.
    terminate_if_other_basin_found : bool
        If true, then if multiple basin values were found, then the ODE solver will exit.
        If false, then the ODE solver will run until all points enter one of the
        beta-sphere/trust-region.

    Returns
    -------
    ndarray(N,):
        Integer array that assigns each point to a basin/maxima/atom of the molecule.
        If value is negative one, then the point wasn't assigned to a basin.

    """
    norm_grad_func = _get_normalized_gradient_func(grad_func)

    numb_pts = initial_pts.shape[0]
    if isinstance(ss_0, float):
        ss = np.ones((numb_pts, 1)) * ss_0
    elif isinstance(ss_0, np.ndarray):
        if not ss_0.shape[1] == 1:
            raise ValueError(f"Steps-size {ss_0.shape} should have shape of the form (N, 1).")
        ss = ss_0.copy()
    else:
        raise TypeError(f"Step-size ss_0 {type(ss_0)} should have type float or array.")

    pts = initial_pts.copy()
    dens_vals0 = dens_func(initial_pts)

    not_found_indices = np.arange(numb_pts)
    niter = 0  # Number of iterations
    grad0 = norm_grad_func(pts)  # Avoids re-computing the gradient twice, used to check for NNA
    pts_curr = pts.copy()
    gradient_paths = [[] for _ in range(len(pts_curr))]
    while len(not_found_indices) != 0:
        if niter == maxiter:
            raise RuntimeError(
                f"Number of iterations reached maximum {niter}, "
                f"this may be because of a non-nuclear attractor (NNA) which may cause the ODE "
                f"to cycle between two points. Repeat this calculation by including the "
                f"non-nuclear attractor to the list of critical points."
            )
        y_four, y_five = _RK45_step(pts_curr, norm_grad_func, ss, grad0)

        # Update step-size
        ss = (tol * ss / (2.0 * np.linalg.norm(y_five - y_four, axis=1)[:, None])) ** 0.25
        ss[ss > max_ss] = max_ss

        # Get density values and if the density-values decreased, reduce step-size
        dens_vals1 = dens_func(y_five)
        indices = np.where(dens_vals1 <= dens_vals0)[0]
        if len(indices) != 0:
            y_five[indices, :] = pts_curr[indices, :]
            ss[indices] *= 0.25

        for i, global_index_pt in enumerate(not_found_indices):
            if i not in indices:
                gradient_paths[global_index_pt].append(y_five[i])

        # Check any points are within the beta-spheres and remove them if they converged.
        converged = np.where(np.all(np.abs(y_five - pts_curr) < tol_conv, axis=1))[0]
        if len(converged) != 0:
            not_found_indices = np.delete(not_found_indices, converged)
            y_five = np.delete(y_five, converged, axis=0)
            ss = np.delete(ss, converged)[:, None]
            dens_vals1 = np.delete(dens_vals1, converged)

        if len(y_five) == 0:
            break
        # Update next iteration
        grad0 = norm_grad_func(y_five)
        pts_curr = y_five.copy()
        dens_vals0 = dens_vals1
        niter += 1

    return gradient_paths
