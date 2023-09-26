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
from grid.cubic import UniformGrid
import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import root_scalar

from chemtools.topology.ode import find_basins_steepest_ascent_rk45, steepest_ascent_rk45, gradient_path

r"""
Utility functions that is common between the QTAIM algorithms.
"""

__all__ = [
    "solve_for_oas_points",
    "construct_radial_grids",
    "find_optimize_centers",
    "determine_beta_spheres_and_nna",
    "solve_intersection_of_ias_point"
]


def solve_for_isosurface_pt(
    l_bnd, u_bnd, maxima, cart_sphere_pt, density_func, iso_val, iso_err
):
    r"""
    Solves for the point on a ray that satisfies the isosurface value equation.
    .. math::
        f(r) := \rho(\textbf{A} + r \textbf{\theta}) - c,
    where A is the position of the atom, :math:`\theta` is the Cartesian coordinates of the
    point on the sphere, r is the radius, and c is the isosurface value.  The radius
    is solved using a root-finding algorithm over an interval that contains the isosurface
    value.
    Parameters
    ----------
    l_bnd: float
        The lower-bound on the radius for the root-solver. Needs to be less than the
        isosurface value.
    u_bnd: float
        The upper-bound on the radius for the root-solver. Needs to be greater than the
        isosurface value.
    maxima: ndarray(3,)
        The maximum of the atom.
    cart_sphere_pt: ndarray(3,)
        The Cartesian coordinates of the point on the sphere.
    density_func: callable(ndarray(M,3), ndarray(M,))
        The electron density function.
    iso_val: float
        The isosurface value.
    iso_err: float
        The xtol for the root-solver.
    Returns
    -------
    ndarray(3,):
        The point :math:`\textbf{A} + r \textbf{\theta}` that satisfies the isosurface value.
    """
    # Given a series of points based on a maxima defined by angles `cart_sphere_pt` with
    #  radial pts `rad_pts`.   The `index_iso` tells us where on these points to construct another
    #  refined grid from finding l_bnd and u_bnd.
    bounds_found = False
    while not bounds_found:
        dens_l_bnd = density_func(np.array([maxima + l_bnd * cart_sphere_pt]))
        dens_u_bnd = density_func(np.array([maxima + u_bnd * cart_sphere_pt]))
        if iso_val < dens_u_bnd or dens_l_bnd < iso_val:
            if iso_val < dens_u_bnd:
                l_bnd = u_bnd
                u_bnd += 1.5
            elif dens_l_bnd < iso_val:
                u_bnd = l_bnd
                l_bnd -= 1.5
            # raise ValueError(f"Radial grid {l_bnd, u_bnd} did not bound {dens_l_bnd, dens_u_bnd} "
            #                  f"the isosurface value {iso_val}. Use larger radial grid.")
        else:
            bounds_found = True

    # Use Root-finding algorithm to find the isosurface point.
    root_func = lambda t: density_func(np.array([maxima + t * cart_sphere_pt]))[0] - iso_val
    sol = root_scalar(root_func, method="toms748", bracket=(l_bnd, u_bnd), xtol=iso_err)
    assert sol.converged, f"Root function did not converge {sol}."
    bnd_pt = maxima + sol.root * cart_sphere_pt
    return bnd_pt


def _solve_root_newton_raphson(initial_guess, angular_pts, root_and_grad, xtol, maxiter=1000):
    not_converged = np.arange(len(initial_guess))
    niter = 0  # Number of iterations
    success = True
    pts0 = initial_guess.copy()
    pts_converged = np.zeros(len(initial_guess))
    while len(not_converged) != 0:
        if niter == maxiter:
            success = False
            break
        dens_pts, grad_pts = root_and_grad(pts0, angular_pts[not_converged])

        pts1 = pts0 - dens_pts / grad_pts
        indices_converged = np.where(np.abs(pts1 - pts0) <= xtol)[0]
        if len(indices_converged) != 0:
            pts_converged[not_converged[indices_converged]] = pts0[indices_converged]

        not_converged = np.delete(not_converged, indices_converged)
        pts0 = pts1
        pts0 = np.delete(pts0, indices_converged)
        niter += 1

    return pts_converged, success


def _solve_root_grid_movement(initial_guess, angular_pts, root_func, iso_err, xtol=1e-3, maxiter=10000):
    not_converged = np.arange(len(initial_guess))
    niter = 0  # Number of iterations
    success = True
    pts0 = initial_guess.copy()
    root_pts = root_func(pts0, angular_pts[not_converged])
    pts_converged = np.zeros(len(initial_guess))
    while len(not_converged) != 0:
        if niter == maxiter:
            success = False
            break

        # Take a step based on the sign of the root equation, and take the step-size to be minimum of xtol
        #       and root_pts, i.e. you want the step-size to be smaller when you're closer to the root.
        pts1 = pts0 + np.sign(root_pts) * np.minimum(xtol, root_pts * 10)
        root_pts = root_func(pts1, angular_pts[not_converged])

        indices_converged = np.where(np.abs(root_pts) < iso_err)[0]
        if len(indices_converged) != 0:
            pts_converged[not_converged[indices_converged]] = pts1[indices_converged]

        not_converged = np.delete(not_converged, indices_converged)
        pts0 = pts1
        pts0 = np.delete(pts0, indices_converged)
        root_pts = np.delete(root_pts, indices_converged)
        niter += 1

    return pts_converged, success


def solve_for_oas_points(
     maximas, maximas_to_do, oas, angular_pts, dens_func, iso_val, iso_err, r_func
):
    r"""
    For each index in outer-atomic surface (OAS) solves for the isovalue point along a ray.

    This is stored inside `r_func`. Solves for the point on a ray that satisfies the isosurface value equation.

    .. math::
        f(r) := \rho(\textbf{A} + r \textbf{\theta}) - c,

    where A is the position of the atom, :math:`\theta` is the Cartesian coordinates of the
    point on the sphere, r is the radius, and c is the isosurface value.  The radius
    is solved using a root-finding algorithm over an interval that contains the isosurface
    value.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        All maximas of the density
    maximas_to_do: list[int]
        List of indices to solve for the OAS for each maximas.
    oas: list[list]
        List of indices that correspond to angular points whose ray intersect the isosurface
        of the electron density.
    radial_grid: list[ndarray]
        List of radial grids (arrays on zero to infinity) correspond to each maxima.
    angular_pts: list[ndarray]
        The angular points on the sphere for each maxima.
    dens_func: callable()
        The electron density function.
    iso_val: float
        The isosurface value that is to be solved
    iso_err: float
        The isosurface error
    r_func: list[ndarray()]
        This holds the radial coordinate on the ray that intersects the OAS.

    """
    for i_maxima in maximas_to_do:
        print("ATOM ", i_maxima)
        maxima = maximas[i_maxima]
        ang_pts = angular_pts[i_maxima]
        initial_guess = np.array([0.2] * len(oas[i_maxima]))

        def root_and_grad(t, angular_pts_max):
            # Using logarithm increases accuracy and convergence as Newton-Ralphson has quadratic convergence.
            real_pts = maxima + t[:, None] * angular_pts_max
            dens_pts = dens_func(real_pts)
            return dens_pts - iso_val

        sol_x, success = _solve_root_grid_movement(
            initial_guess, ang_pts[oas[i_maxima]], root_and_grad, xtol=0.01, iso_err=iso_err
        )
        radial_results = sol_x
        sol_fun = dens_func(maxima + sol_x[:, None] * ang_pts[oas[i_maxima]])
        # The points that weren't successful, try again.
        if not success:
            # Get rid of the points that converge, and re-try with the points that didn't.
            print("Try solving the root equations for OAS again on individual points.")
            indices = np.where(np.abs(sol_fun) > iso_err)[0]

            for i_oas in indices:
                oas_pt = solve_for_isosurface_pt(
                    radial_results[i_oas] - 0.1, radial_results[i_oas] + 0.1, maxima, ang_pts[i_oas],
                    dens_func, iso_val, iso_err
                )
                radial_results[i_oas] = np.linalg.norm(oas_pt - maxima)
        r_func[i_maxima][oas[i_maxima]] = radial_results


def find_optimize_centers(centers, grad_func):
    maximas = np.array(
        [gradient_path(x, grad_func, t_span=(0, 10), method="BDF",
                       first_step=1e-9, max_step=1e-1) for x in centers],
        dtype=np.float64
    )
    print("New maximas: \n ", maximas)
    # Check for duplicates
    distance = cdist(maximas, maximas)
    distance[np.diag_indices(len(maximas))] = 1.0  # Set diagonal elements to one
    if np.any(distance < 1e-6):
        raise RuntimeError(f"Optimized maximas contains duplicates: \n {maximas}.")
    return maximas


def construct_radial_grids(pts1, maximas, min_pts=0.2, pad=5.0, ss0=0.23):
    r"""
    Construct radial grids around each maxima depending on its neighbors

    Parameters
    ----------
    pts: ndarray(M_1, 3)
        Coordinates of the points to construct radial grids on.
    maximas: ndarray(M, 3)
        Coordinates of the maximas.
    min_pts: Union[float, list]
        The minimum radial value on [0, \infty) to construct a radial grid. If it is a list,
        it should be of size `M`.
    pad: float
        Extra padding to add to make sure the radial grid covers the intersection
        with the inter-atomic and outer-atomic surfaces.
    ss0: float
        The step-size of the uniform radial grid.


    Returns
    -------
    list[ndarray]:
        List of radial grid of length `M_1`. The radial grid are uniform
        grids that start at `min` and end on the maximum distance to the fifth atom plus
        an extra padding with stepsize `ss0`.

    """
    # Construct a radial grid for each atom by taking distance to the closest five atoms.
    #  Added an extra padding in the case of carbon in CH4
    #  TODO: the upper-bound should depend on distance to isosurface value and distance
    #         between atoms
    if isinstance(min_pts, (float, np.float)):
        # If it is a float, convert it to a list.
        min_pts = [min_pts] * len(maximas)
    dist_maxs = cdist(pts1, maximas)
    sorted_dists = np.sort(dist_maxs, axis=1)
    distance_maximas = sorted_dists[:, min(5, maximas.shape[0] - 1)]
    distance_minimas = sorted_dists[:, 1] / 4.0
    radial_grid = [
        np.arange(min(min_pts[i], distance_minimas[i]), x + pad, ss0) for i, x in enumerate(distance_maximas)
    ]
    return radial_grid


def determine_beta_spheres_and_nna(
        beta_spheres, maximas, radial_grids, angular_pts, dens_func, grad_func, hess_func=None
):
    r"""

    Notes this assumes the initial beta-sphere is 0.01, and so the distance between maximas
    cannot be smaller than this.

    """
    numb_maximas = len(maximas)
    initial_beta_sph = 0.01
    if beta_spheres is None:
        beta_spheres = [initial_beta_sph] * numb_maximas
    # TODO: add the sphere angle trick.

    # Determine the beta-spheres
    i_maxima = 0
    while i_maxima != len(maximas):
        maxima = maximas[i_maxima]
        if beta_spheres[i_maxima] == initial_beta_sph:
            for rad_pt in radial_grids[i_maxima]:
                if rad_pt > initial_beta_sph:
                    # Determine the points on the sphere with this radius
                    pts = maxima + rad_pt * angular_pts
                    # You want here for the ODE to be accurate in-order to find potential NNA.
                    if hess_func is None:
                        basins, _ = find_basins_steepest_ascent_rk45(
                            pts, dens_func, grad_func, beta_spheres, maximas, ss_0=0.1, max_ss=0.2, tol=1e-9,
                            hess_func=hess_func
                        )
                    else:
                        basins, maximas = find_basins_steepest_ascent_rk45(
                            pts, dens_func, grad_func, beta_spheres, maximas, ss_0=0.2, max_ss=0.2, tol=1e-9,
                            hess_func=hess_func, check_for_nna=True
                        )
                    basins = np.array(basins, dtype=np.int)

                    which_nna = np.where(basins >= numb_maximas)[0]
                    if len(which_nna) != 0:
                        # TODO: Easy way to determin min_pt is to take minimum distance between NNA and other maximas.
                        # Copy a radial grid from the previous method
                        radial_grids += \
                            construct_radial_grids(maximas[numb_maximas:], maximas[:numb_maximas],
                                                   min_pts=0.01, pad=5.0, ss0=0.01)

                        print(maximas)
                        print(beta_spheres)
                        beta_spheres += [initial_beta_sph] * (len(maximas) - numb_maximas)
                        numb_maximas = len(maximas)
                        #input("Found NNA")

                    # If all the basins went to same maxima, then update radius
                    # If the first radial point didn't suffice(e.g. NNA), then set the beta-sphere to something small.
                    # else then break out of this for loop.
                    if np.all(basins == i_maxima):
                        optimal_rad = rad_pt
                        beta_spheres[i_maxima] = optimal_rad
                        print(beta_spheres)
                        print("Optimal radius is ", optimal_rad)
                    elif np.abs(rad_pt - radial_grids[i_maxima][0]) < 1e-8:
                        beta_spheres[i_maxima] = min(rad_pt / 4.0, 0.3)
                        break
                    else:
                        break
        print(f"i Maxima {i_maxima} and Final optimal radius {beta_spheres[i_maxima]}")
        i_maxima += 1
        # input("next maxima")
    return beta_spheres, maximas, radial_grids


def find_non_nuclear_attractors(maximas, dens_func, grad_func, hess_func):
    r"""
    Finds non-nuclear attractors up to two decimal places.

    Returns
    -------
    ndarray(M+K, 3)
        Returns the original `M` maximas from `maximas` and adds `K` new maximas that are
        non-nuclear attractors.

    """
    grid = UniformGrid.from_molecule(
        [6.0] * len(maximas), maximas, spacing=0.1, rotate=False, extension=2.0,
    )

    dens_vals = dens_func(grid.points)
    indices2 = np.where(dens_vals > 0.001)[0]
    grads_vals = grad_func(grid.points[indices2])
    print(grads_vals)
    indices = np.where(np.all(np.abs(grads_vals) < 0.01, axis=1))[0]
    print(indices2)
    print(indices)
    if len(indices) != 0:
        np.set_printoptions(threshold=np.inf)
        pts = steepest_ascent_rk45(
            grid.points[indices2][indices], dens_func, grad_func, tol_conv=1e-8
        )
        print(pts)
        dist_maxima = cdist(pts, maximas)
        print(dist_maxima)
        beta_sph = dist_maxima <= [0.1] * len(maximas)
        which_basins = np.where(beta_sph)

        pts = np.delete(pts, which_basins[0], axis=0)
        print("Pts ")
        if len(pts) != 0:
            pts = np.unique(np.round(pts, decimals=1), axis=0)
            # print(pts)
            nna_attractors = np.array(
                [gradient_path(x, grad_func, t_span=(0, 30), method="BDF",
                               first_step=1e-9, max_step=1e-1) for x in pts],
                dtype=np.float64
            )
            # Round to two decimal places mostly due to empirical evidence of convergence of these ODE.
            nna_attractors = np.unique(np.round(nna_attractors, 2), axis=0)
            print(nna_attractors)
            eigs = np.linalg.eigvalsh(hess_func(nna_attractors))
            print("eigs", eigs)
            which_is_nna = np.where(np.all(eigs < -1e-10, axis=1))[0]
            print(which_is_nna)

            maximas = np.vstack((maximas, nna_attractors[which_is_nna]))
            print(maximas)
    return maximas


def solve_intersection_of_ias_point(
    maximas, ias_indices, angular_pts,
        dens_func, grad_func, beta_spheres, bnd_err, ias_lengths, ss_0, max_ss, tol, hess_func=None
):
    r"""
    Solves the intersection of the ray to the inner-atomic surface.

    A point is associated to each ray based on `ias_indices`.  The basin value
    is assigned to each point, and the point is moved along the ray until it keeps switching basins.
     The process is further repeated with a smaller step-size until the distance between two
     points on the ray is less than `bnd_err`.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        Optimized centers of the electron density.
    ias_indices: ndarray(N, 6)
       Rows correspond to each ray that intersects the IAS.
       First index is which index of maxima it originates from, then second
       index is index of angular point/ray, third index is the lower bound radius
       and fourth index is the upper-bound radius, fifth index is the step-size.
       The fifth index holds which basin the IAS point switches to. Note it may
       not be the true basin value that it switches to.
       The sixth index holds which index of the `IAS` list it points to.
    angular_pts: list[ndarray]
        List of size `M` of the angular points over each maxima.
    dens_func:
        The density of electron density.
    grad_func:
        The gradient of the electron density.
    beta_spheres: ndarray(M,)
        Beta-spheres radius of each atom.
    bnd_err: float
        The error of the intersection of the IAS. When the distance to two consequent points
        on the ray crosses different basins is less than this error, then the midpoint is accepted
        as the final radius value.
    ias_lengths: list[int]
        List of length `M` atoms that contains the number of IAS points
    ss_0: float, optional
        The initial step-size of the ODE (RK45) solver.
    max_ss: float, optional
        Maximum step-size of the ODE (RK45) solver.
    tol: float, optional
        Tolerance for the adaptive step-size.

    Return
    -------
    List[ndarray()], list[list[int]]:
        The first list holds arrays that gives the radius value that intersects the IAS or OAS.
        The second list of size `M`, holds which ias crosses which other basin.

    """
    if not isinstance(ias_indices, np.ndarray):
        raise TypeError(f"The parameters to solve indices should be numpy array instead of {type(ias_indices)}.")
    r_func = [np.zeros((len(angular_pts[i]),), dtype=np.float64) for i in range(len(maximas))]
    basin_ias = [[-1] * ias_lengths[i] for i in range(len(maximas))]

    while len(ias_indices) != 0:
        # Construct New Points
        points = []
        for (i_maxima, i_ang, l_bnd, u_bnd, _, _, _) in ias_indices:
            # Take the midpoint of interval [l_bnd, u_bnd]
            ray = maximas[int(i_maxima)] + (l_bnd + u_bnd) / 2.0 * angular_pts[int(i_maxima)][int(i_ang), :]
            points.append(ray)
        points = np.vstack(points)

        # Solve for basins
        basins, _ = find_basins_steepest_ascent_rk45(
            points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0, hess_func=hess_func
        )
        # print("Basins", basins)

        # Refine the rays further
        # print("Average step-size", np.mean(ias_indices[:, 4]))
        converge_indices = []
        # Note it always assumes that `l_bnd` goes to `i_maxima` and `u_bnd` goes to `basin_switch`.
        for i, (i_maxima, i_ang, l_bnd, u_bnd, ss, basin_switch, i_ias) in enumerate(ias_indices):
            #print((i_maxima, i_ang, l_bnd, u_bnd, ss, basin_switch))
            if basins[i] == i_maxima:
                # Calculate new step-size between [(l_bnd + u_bnd) / 2,  u_bnd]
                new_l_bnd = (l_bnd + u_bnd) / 2
                new_u_bnd = u_bnd
            else:
                # Calculate new stepsize between [l_bnd, (l_bnd + u_bnd) / 2]
                new_l_bnd = l_bnd
                new_u_bnd = (l_bnd + u_bnd) / 2
                # Update the basin that it switches to.
                if basin_switch != basins[i]:
                    basin_switch = basins[i]
            new_ss = (new_u_bnd - new_l_bnd)
            # If new stepsize was less than bnd_err, then we converge and should stop.
            if new_ss <= bnd_err:
                r_func[int(i_maxima)][int(i_ang)] = (new_l_bnd + new_u_bnd) / 2
                basin_ias[int(i_maxima)][int(i_ias)] = int(basin_switch)
                converge_indices.append(i)  # Put in list to remove indices.
            else:
                # Update ias_indices for the next iteration for example l_bnd, u_bnd, step-size
                ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, new_ss, basin_switch, i_ias]

        # Remove converged indices
        # print("Convergence indices", converge_indices)
        ias_indices = np.delete(ias_indices, converge_indices, axis=0)

    # Solve for multiple intersections
    return r_func, basin_ias
