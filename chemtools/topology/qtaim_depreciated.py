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

from scipy.integrate import solve_ivp
from scipy.interpolate import LSQSphereBivariateSpline, SmoothSphereBivariateSpline
from scipy.optimize import root_scalar
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist

from scipy.sparse import lil_matrix

from grid.atomgrid import AtomGrid
from grid.cubic import UniformGrid, _HyperRectangleGrid
from grid.angular import AngularGrid
from grid.utils import convert_cart_to_sph

from chemtools.topology.ode import steepest_ascent_rk45
from chemtools.topology.surface import SurfaceQTAIM

import time


def construct_points_between_ias_and_oas(
    ias: list, oas: int, angular_pts: np.ndarray, r_func_max: np.ndarray, maxima: np.ndarray
):
    r"""
    Construct points between the inner atomic surface and outer atomic surface.

    This is done by constructed a convex hull between IAS and OAS, seperetely.
    Each point on the IAS, the two closest points are found on the OAS, then
    a triangle is constructed.  Seven points are constructed within this triangle
    and the Cartesian coordinates of the sphere centered at the maxima is solved
    for each of these seven points.

    Parameters
    -----------
    ias : List[int]
        List of integers of `angular_pts` that are part of the inner atomic surface (IAS).
    oas : List[int]
        List of integers of `angular_pts` that are part of the outer atomic surface (OAS).
    angular_pts : np.ndarray
        Angular Points around the maxima for which rays are propgated from.
    r_func_max : np.ndarray
        The radial component for each angular point in `angular_pts` that either gives
        the radial value that intersects the OAS or the IAS.
    maxima : np.ndarray
        Maxima of the basin.

    Returns
    -------
    ndarray(K * 7, 3)
        Cartesian coordinates of :math:`K` points on the sphere centered at `maxima` such that
        they correspond to the seven points constructed above, where :math:`K` is the number
        of points on the IAS of `maxima`.

    """
    # Take a convex hull of both IAS and OAS seperately.
    ias_pts = maxima + r_func_max[ias, None] * angular_pts[ias, :]
    oas_pts = maxima + r_func_max[oas, None] * angular_pts[oas, :]
    ias_hull = ConvexHull(ias_pts)
    oas_hull = ConvexHull(oas_pts)
    ias_bnd = ias_hull.points[ias_hull.vertices]
    oas_bnd = oas_hull.points[oas_hull.vertices]

    # Compute the distance matrix
    dist_mat = cdist(ias_bnd, oas_bnd)
    # for each point in say ias take the closest two points in oas.
    new_ang_pts = np.zeros((0, 3), dtype=np.float64)  # usually 7 points per ias boundary are added.
    for i_ias, pt_ias in enumerate(ias_bnd):
        # Get the two closest points on OAS to this IAS pt.
        two_indices = dist_mat[i_ias].argsort()[:2]
        pt1, pt2 = oas_bnd[two_indices[0]], oas_bnd[two_indices[1]]

        # Take the center and midpoint between each line of the triangle (pt_ias, pt1, pt2)
        midpoint = (pt1 + pt2 + pt_ias) / 3.0
        line_pt1 = (pt1 + pt_ias) / 2.0
        line_pt2 = (pt2 + pt_ias) / 2.0
        line_pt3 = (pt1 + pt2) / 2.0

        # The triangle with the center can be split into three polygons, take the center of each.
        poly_pt1 = (midpoint + line_pt1 + line_pt2 + pt_ias) / 4.0
        poly_pt2 = (midpoint + line_pt1 + line_pt3 + pt1) / 4.0
        poly_pt3 = (midpoint + line_pt2 + line_pt3 + pt2) / 4.0

        new_pts = np.array([midpoint, line_pt1, line_pt2, line_pt3, poly_pt1, poly_pt2, poly_pt3])
        # Solve for the Cartesian angular coordinates of these 7 points by solving
        #  r = m + t direction, where m is the maxima, direction has norm one, r is each of
        #  these points

        # print("new pts ", new_pts)
        # matplotlib.use("Qt5Agg")
        # fig = plt.figure()
        # ax = plt.axes(projection='3d')
        # p = ias_bnd
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
        # p = oas_bnd
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
        # p = new_pts
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="y")
        # ax.scatter(pt1[0], pt1[1], pt1[2], color="m", s=30)
        # ax.scatter(pt2[0], pt2[1], pt2[2], color="m", s=30)
        # plt.show()

        direction = new_pts - maxima
        # Delete points that are on the maxima.
        direction = np.delete(direction, np.all(np.abs(direction) < 1e-10, axis=1), axis=0)
        print("Direction ", direction)
        t = np.linalg.norm(direction, axis=1)
        direction = direction / t[:, None]
        # Delete directions that are the same
        direction = np.unique(direction, axis=0)
        new_ang_pts = np.vstack((new_ang_pts, direction))
    return new_ang_pts


def gradient_path(pt, grad_func, maximas=None, t_span=(0, 1000), method="LSODA", max_step=100,
                  t_inc=400, max_tries=10, first_step=1e-3, beta_spheres=-np.inf):
    # TODO: If the density value is low, normalized_gradient low and trying ODE did not move much,
    #  then an option is to turn max_step tp np.inf
    is_converged = False
    y0 = pt.copy()
    # print("PT ", pt)
    numb_times = 0
    def grad(t, x):
        return grad_func(np.array([x]))[0]

    while not is_converged and numb_times < max_tries:
        sol = solve_ivp(
            grad, #lambda t, x: grad_func(np.array([x]))[0].T,
            y0=y0,
            t_span=t_span,
            method=method,
            max_step=max_step,
            first_step=first_step,
            # vectorized=True,
        )
        # print(sol)
        assert sol["success"], "ODE was not successful."
        # TODO: Write in docs that it summes all local maximas are identified.
        # If it is close to a maxima or within any of the beta-spheres, then stop.
        if maximas is not None:
            last_y_val = sol["y"][:, -1]
            dist_maxima = np.linalg.norm(last_y_val - maximas, axis=1)
            if np.any(dist_maxima < 0.1) or np.any(dist_maxima <= beta_spheres):
                return sol["y"][:, -1]
        # if maximas not specified, then just look at if it converged.
        else:
            convergence = np.linalg.norm(sol["y"][:, -2] - sol["y"][:, -1])
            if convergence < 1e-1:
                return sol["y"][:, -1]
        # No convergence occured, so increaes t-span.
        # print(sol["y"][:, -1], t_span, "YE")
        t_span = (t_span[1], t_span[1] + t_inc)
        y0 = sol["y"][:, -1]
        numb_times += 1

    if numb_times == max_tries:
        raise RuntimeError(f"No convergence in normalized_gradient path pt {pt},"
                           f" solution {sol['y'][:, -1]}, t_span {t_span}")


def gradient_path_all_pts(
    pts, grad_func, beta_spheres, i_maxima, maximas=None, t_span=(0, 1000), method="LSODA", max_step=100,
    t_inc=500, first_step=1e-3, rtol=1e-4, atol=1e-7, is_watershed_pt=False
):
    # i_maxima is the index of the maxima where the ray originates from
    def grad(t, x):
        return grad_func(np.array([x]))[0]

    basins = np.zeros((len(pts),))
    for i_pt, pt in enumerate(pts):
        found_basin = False
        y0 = pt.copy()
        print("Pt to backtrace", pt)

        while not found_basin:
            sol = solve_ivp(
                grad,
                y0=y0,
                t_span=t_span,
                method=method,
                max_step=max_step,
                first_step=first_step,
                rtol=rtol,
                atol=atol,
            )
            assert sol["success"], "ODE was not successful."
            y_vals = sol["y"][:, -1]

            # See if any of the points converge to their beta-spheres.
            dist_maxima = cdist(np.array([y_vals]), maximas)
            beta_sph = dist_maxima <= beta_spheres
            # print("Dist maxima ", dist_maxima)
            # print("Beta-Sphere ", beta_sph, beta_spheres)
            if np.any(beta_sph):
                which_basin = np.where(beta_sph[0])
                assert len(which_basin[0]) == 1, "More than one basin was found"
                print("Which basin ", which_basin)
                found_basin = True
                basins[i_pt] = which_basin[0][0]

                # If it is a point that is guaranteed to be watershed point
                # then stop when the point that it isn't i_maxima is found.
                if is_watershed_pt and which_basin[0][0] != i_maxima:
                    return basins

            # Could not found basin, so update
            t_span = (t_span[1], t_span[1] + t_inc)
            y0 = y_vals
    return basins


def gradient_path_vectorized(pts, grad_func, maximas=None, t_span=(0, 1000), method="LSODA", max_step=100,
                  t_inc=400, max_tries=10, first_step=1e-3, beta_spheres=-np.inf, rtol=1e-4, atol=1e-7):
    y0 = np.ravel(pts, order="C")
    numb_pts = len(pts)
    print("Numb_pts ", numb_pts)
    numb_tries = 0

    def grad(t, pts, numb_rays):
        pts_arr = np.reshape(pts, (numb_rays, 3), order="C")
        pts_arr = pts_arr.copy()
        return np.ravel(grad_func(pts_arr), order="C")

    indices = np.arange(0, numb_pts)  # indices not converged
    basins = np.zeros((numb_pts), dtype=int)
    while len(indices) != 0:
        sol = solve_ivp(
            grad,
            y0=y0,
            t_span=t_span,
            method=method,
            max_step=max_step,
            first_step=first_step,
            args=(numb_pts,),
            rtol=rtol,
            atol=atol,
            vectorized=False,
        )
        assert sol["success"]
        # print(sol)
        y_vals = sol["y"][:, -1]
        print("Yvals", np.reshape(y_vals, (numb_pts, 3)))
        y_vals = np.reshape(y_vals, (numb_pts, 3), order="C")

        # See if any of the points converge to maximas or their beta-spheres.
        dist_maxima = cdist(y_vals, maximas)
        print("Dist Maxima", dist_maxima)
        print("Which less than 0.1 ", np.any(dist_maxima < 0.1, axis=1))
        conv_to_maxima = np.any(dist_maxima < 0.1, axis=1)

        which_basins = np.argmin(dist_maxima, axis=1)  # which maximas it is closest to
        print("Which basins it converged to ", which_basins)

        beta_sph = dist_maxima <= beta_spheres

        print("Dist maxima <= beta_sphere ", beta_sph)
        which_beta_basins = np.where(dist_maxima <= beta_spheres)
        print("which pts are within basin based on beta-sphere", which_beta_basins)
        conv_to_beta = np.any(beta_sph, axis=1)
        print("Conv to maxima", conv_to_maxima)
        print("Conv to bet asphere", conv_to_beta)
        which_converged = (conv_to_maxima | conv_to_beta)
        print("which converged ", which_converged)
        print(np.argmin(which_converged, axis=0))

        # Update which basins it converged to
        basins[indices[which_converged]] = which_basins[which_converged]
        if len(which_beta_basins[1]) != 0:
            # If the distance to beta-sphere where found, then replace it with those values.
            print("basins", basins)
            print("indices", indices)
            print("which converged", which_converged)
            print(which_beta_basins[1], conv_to_beta)
            basins[indices[conv_to_beta]] = which_beta_basins[1]
        print("Basins ", basins)

        # delete indices that converged
        indices = np.delete(indices, which_converged)
        y_vals = np.delete(y_vals, which_converged, axis=0)
        print("indices didn't converge: ", indices)

        # the rest are continued increasing the t_span accordingly.
        numb_pts = len(indices)
        y0 = np.ravel(y_vals, order="C")
        t_span = (t_span[1], t_span[1] + t_inc)
        numb_tries += 1

        if numb_tries == max_tries:
            raise RuntimeError(f"No convergence in normalized_gradient path pt {y0},"
                               f" solution {sol['y'][:, -1]}, t_span {t_span}")
        # input("sds")
    return basins


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
    dens_l_bnd = density_func(np.array([maxima + l_bnd * cart_sphere_pt]))
    dens_u_bnd = density_func(np.array([maxima + u_bnd * cart_sphere_pt]))
    if iso_val < dens_u_bnd or dens_l_bnd < iso_val:
        if iso_val < dens_u_bnd:
            u_bnd += 1.5
        elif dens_l_bnd < iso_val:
            l_bnd -= 1.5
        # raise ValueError(f"Radial grid {l_bnd, u_bnd} did not bound {dens_l_bnd, dens_u_bnd} "
        #                  f"the isosurface value {iso_val}. Use larger radial grid.")

    # Use Root-finding algorithm to find the isosurface point.
    root_func = lambda t: density_func(np.array([maxima + t * cart_sphere_pt]))[0] - iso_val
    sol = root_scalar(root_func, method="toms748", bracket=(l_bnd, u_bnd), xtol=iso_err)
    assert sol.converged, f"Root function did not converge {sol}."
    bnd_pt = maxima + sol.root * cart_sphere_pt
    return bnd_pt


def solve_for_basin_bnd_pt(
    dens_cutoff, i_maxima, maximas, radial, cart_sphere_pt, dens_func, grad_func, bnd_err,
    iso_val, beta_spheres
):
    # Construct the ray and compute its density values based on a maxima defined by angles
    #  `cart_sphere_pt` with radial pts `rad_pts`.    It goes through each point on the ray
    #   if the ray density value is greater than dens_cutoff, then it is likely this ray
    #   tends towards infinity and has no basin boundary.  If density value is larger, then
    #   it solves for the normalized_gradient path via solving normalized_gradient ode.  If this ode solution,
    #   is close to other basins, then we found when it switched basins.  Then we take
    #  the two points where it switches basin and compute the distance, if this distance
    #  is less than `bnd_err`, then we take the midpoint to be the boundary point on the ray
    #  that intersects the ias.  If not, then we construct a new ray with different l_bnd
    #  and u_bnd and reduce the step-size further and repeat this process.
    rad_pts = radial.copy()
    ss_ray = np.mean(np.diff(rad_pts))   # Stay with a coarse ray then refine further.
    index_iso = None  # Needed to refine if the ray tends towards infinity.
    bnd_pt = None  # Boundary or Isosurface Point
    is_ray_to_inf = False  # Does this ray instead go towards infinity

    found_watershed_on_ray = False
    is_watershed_pt = False
    basin_id = None
    counter = 0
    while not found_watershed_on_ray:
        ray = maximas[i_maxima] + rad_pts[:, None] * cart_sphere_pt
        ray_density = dens_func(ray)
        print("Start of Ray ", ray[0], " Cartesian pt of Sphere ", cart_sphere_pt, "Final Ray Pt: ",
              ray[-1])

        # Cut off ray points that are less than dens_cutoff.
        ray_cutoff = ray_density < dens_cutoff
        ray = np.delete(ray, ray_cutoff, axis=0)
        ray_density = np.delete(ray_density, ray_cutoff)

        # print("The Ray", ray)
        # print("The Ray Density ", ray_density)
        grad_norm = np.min(np.linalg.norm(grad_func(ray), axis=1))
        print("Range ", min(max(0.5 / grad_norm, 10), 50))

        import time
        start = time.time()
        # basins = gradient_path_all_pts(
        #     ray, grad_func, beta_spheres, i_maxima, maximas,
        #     t_span=(0, min(max(0.5 / grad_norm, 10), 50)),
        #     max_step=np.inf,
        #     first_step=1e-6,
        #     method="LSODA",
        #     rtol=bnd_err,
        #     atol=1e-6,
        #     is_watershed_pt=is_watershed_pt
        # )
        basins = steepest_ascent_rk45(
            ray, dens_func, grad_func, beta_spheres, maximas, tol=1e-6
        )
        final = time.time()
        print("Difference ", final - start)
        print("basins ", basins)

        # If they all converged to the same basins, then it is ray to infinity.
        if np.all(basins == i_maxima):
            if counter == 0:
                is_ray_to_inf = True
                index_iso = np.argsort(np.abs(ray_density - iso_val))[0]
                print("Ray to infinity with index ", index_iso)
                break
            else:
                raise ValueError("Went from intersecting basin boundary to classifying as ray to inf.")
        elif len(np.unique(basins)) == 1:
            raise ValueError(f"Went from intersecting basin boundary to classifying as ray to inf."
                             f"{basins}.")

        # if some converged to other basins then refine further
        # first, find which points it switched from basin 1 to basin 2.
        i_switch = np.argmin(basins == i_maxima)
        # print("i_switch", i_switch)
        dist = np.linalg.norm(ray[i_switch - 1, :] - ray[i_switch, :])
        # print("Dist ", dist, np.linalg.norm(rad_pts[i_switch - 1] - rad_pts[i_switch]))

        # If the distance between the two points is less than bnd_err, then stop else refine.
        if np.abs(dist - bnd_err) < 1e-8:
            # Take the midpoint to be the boundary point.
            found_watershed_on_ray = True
            bnd_pt = (ray[i_switch] + ray[i_switch - 1]) / 2.0
            basin_id = basins[i_switch] # basins starts at 1
            print("Found the Watershed point ", bnd_pt, basin_id)
        else:
            # Refine Grid Further
            l_bnd = np.linalg.norm(ray[i_switch - 1] - maximas[i_maxima]) if i_switch != 0 else rad_pts[0] - 1e-3
            u_bnd = np.linalg.norm(ray[i_switch] - maximas[i_maxima])
            ss_ray = max(ss_ray / 10.0, bnd_err)  # Decrease step-size.
            rad_pts = np.arange(l_bnd, u_bnd + ss_ray, ss_ray)
            print("Refine the ray further with l_bnd, u_bnd, ss: ", l_bnd, u_bnd, ss_ray, rad_pts[-1])
            is_watershed_pt = True  # Update that it is a watershed point
        counter += 1   # increment counter so that it doesn't check if entire ray goes to infity.
        # input("Refine further")
    return bnd_pt, is_ray_to_inf, index_iso, found_watershed_on_ray, basin_id


def _optimize_centers(centers, grad_func):
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


def determine_beta_spheres(beta_spheres, maximas, radial_grid, angular_pts, dens_func, grad_func):
    r"""

    Notes this assumes the initial beta-sphere is 0.01, and so the distance between maximas
    cannot be smaller than this.

    """
    numb_maximas = len(maximas)
    initial_beta_sph = 0.01
    if beta_spheres is None:
        beta_spheres = [initial_beta_sph] * numb_maximas
    # Determine the beta-spheres
    for i_maxima, maxima in enumerate(maximas):
        if beta_spheres[i_maxima] == initial_beta_sph:
            optimal_rad = -np.inf
            for rad_pt in radial_grid[i_maxima]:
                if rad_pt > initial_beta_sph:
                    # Determine the points on the sphere with this radius
                    pts = maxima + rad_pt * angular_pts
                    print(pts)

                    # basins = gradient_path_all_pts(
                    #     pts, grad_func, beta_spheres, i_maxima, maximas,
                    #     t_span=(0, 100), max_step=np.inf, method="LSODA",
                    #     first_step=1e-7
                    # )
                    basins = steepest_ascent_rk45(
                        pts, dens_func, grad_func#, beta_spheres, maximas
                    )
                    basins = np.array(basins, dtype=np.int)
                    # If all the basins went to same maxima, then update radius
                    # else then break out of this for loop.
                    if np.all(basins == i_maxima):
                        optimal_rad = rad_pt
                        beta_spheres[i_maxima] = optimal_rad
                        print(beta_spheres)
                        print("Optimal radius is ", optimal_rad)
                    else:
                        break
        print("optimal radius", optimal_rad)
        # input("next maxima")
    return beta_spheres


def qtaim_surface(angular, centers, dens_func, grad_func, iso_val=0.001,
                  dens_cutoff=1e-5, bnd_err=1e-4, iso_err=1e-6,
                  beta_spheres=None, optimize_centers=True, refine=False):
    r"""
    Find the outer atomic and inner atomic surface based on QTAIM.

    For each maxima, a sphere is determined based on `angular` and for each
     point on the angular/sphere, a ray is created based on the radial grid `rgrids`.
     The ray is then determines to either go to infinity and cross the isosurface of the
     electron density or the ray intersects the inner-atomic surface (IAS) of another basin.
     This is determined for each point on the sphere.

    Parameters
    ----------
    angular: List[int] or ndarray(N, 3)
        Either integer specifying the degree to construct angular/Lebedev grid around each maxima
        or array of points on the sphere in Cartesian coordinates.
    centers: ndarray(M,3)
        List of local maximas of the density.
    dens_func: Callable(ndarray(N,3) ->  ndarray(N,))
        The density function.
    grad_func: Callable(ndarray(N,3) -> ndarray(N,3))
        The normalized_gradient of the density function.
    iso_val: float
        The isosurface value of the outer atomic surface.
    dens_cutoff: float
        Points on the ray whose density is less than this cutoff are ignored.
    bnd_err: float
        This determines the accuracy of points on the inner atomic surface (IAS) by controlling
        the step-size of the ray that cross the IAS.
    iso_err: float
        The error associated to points on the OAS and how close they are to the isosurface value.
    beta_spheres : list[float]
        List of size `M` of radius of the sphere centered at each maxima. It avoids backtracing
        of points within the circle. If None is provided, then beta-sphere is determined
        computationally.
    optimize_centers: bool
        If true, then it will optimize the centers/maximas to get the exact local maximas.
    refine : (bool, int)
        If true, then additional points between the IAS and OAS are constructed, added and
        solved for whether it is on the IAS or OAS.

    Returns
    -------
    SurfaceQTAIM
        Class that contains the inner-atomic surface, outer-atomic surface for each maxima.

    Notes
    -----
    - It is possible for a Ray to intersect the zero-flux surface but this algorithm will
        classify it as a ray to infinity because the points on the other side of the basin have
        density values so small that the ode doesn't converge to the maxima of the other basin.
        In this scenario it might be worthwhile to have a denser radial grid with less points
        away from infinity or have a smaller density cut-off.  Alternative for the developer,
        is to implement highly accurate ode solver at the expense of computation time.

    """
    if not isinstance(refine, (bool, int)):
        raise TypeError(f"Refine {type(refine)} should be either boolean or integer.")
    if dens_cutoff > iso_val:
        raise ValueError(f"Density cutoff {dens_cutoff} is greater than isosurface val {iso_val}.")
    if beta_spheres is not None and len(centers) != len(beta_spheres):
        raise ValueError(
            f"Beta sphere length {len(beta_spheres)} should match the"
            f" number of centers {len(centers)}"
        )

    # Using centers, update to the maximas
    maximas = centers
    if optimize_centers:
        # Using ODE solver to refine the maximas further.
        maximas = _optimize_centers(maximas, grad_func)

    # Construct a radial grid for each atom by taking distance to the closest five atoms.
    #  Added an extra padding in the case of carbon in CH4
    #  TODO: the upper-bound should depend on distance to isosurface value and distance
    #         between atoms
    dist_maxs = cdist(maximas, maximas)
    distance_maximas = np.sort(dist_maxs, axis=1)[:, min(5, maximas.shape[0] - 1)]
    print(cdist(maximas, maximas))
    print(distance_maximas + 5.0)
    ss0 = 0.23
    radial_grid = [
        np.arange(0.2, x + 5.0, ss0) for x in distance_maximas
    ]
    input("Hello")

    numb_maximas = len(maximas)
    angular_pts = AngularGrid(degree=angular).points if isinstance(angular, int) else angular
    r, thetas, phis = convert_cart_to_sph(angular_pts).T
    numb_ang_pts = len(thetas)

    # Determine beta-spheres from a smaller angular grid
    #  Degree can't be too small or else the beta-radius is too large and a IAS poitn got
    #  classified as a OAS point
    ang_grid = AngularGrid(degree=10)
    # TODO: Do the spherical trick then do the beta-sphere
    if beta_spheres is None:
        beta_spheres = determine_beta_spheres(
            beta_spheres, maximas, radial_grid, ang_grid.points, dens_func, grad_func
        )
        beta_spheres = np.array(beta_spheres)
    # Check beta-spheres are not intersecting
    condition = dist_maxs <= beta_spheres[:, None] + beta_spheres
    condition[range(len(maximas)), range(len(maximas))] = False  # Diagonal always true
    if np.any(condition):
        raise ValueError(f"Beta-spheres {beta_spheres} overlap with one another.")

    r_func = [np.zeros((numb_ang_pts,), dtype=np.float64) for _ in range(numb_maximas)]
    oas = [[] for _ in range(numb_maximas)]  # outer atomic surface
    ias = [[] for _ in range(numb_maximas)]  # inner atomic surface.
    basin_ias = [[] for _ in range(numb_maximas)]  # basin ids for inner atomic surface.
    refined_ang = [] if refine else None
    maxima_to_do = range(0, numb_maximas) if type(refine) == type(True) else [refine]  # refining
    for i_maxima, maxima in enumerate(maximas):
        # Maximas aren't usually large, so doing this is okay. Quick fix to use refinement without
        #  re-writing this function into seperate functions.
        if i_maxima in maxima_to_do:
            print("Start: Maxima ", maxima)

            # First classify points as either watershed/IAS or isosurface/OAS
            # Each angular point would have a different radial grid associated with it.
            radial = radial_grid[i_maxima][radial_grid[i_maxima] >= beta_spheres[i_maxima]]
            ias_indices = []
            ias_basin = []
            ias_radius = []  # Radius to start at
            indices_to_classify = np.arange(len(angular_pts))
            for i_rad in range(0, len(radial)):  # Go through each radial shell
                # Construct points on the angular points that aren't classified yet.
                all_points = maxima + radial[i_rad, None] * angular_pts[indices_to_classify, :]

                start = time.time()
                # basins = steepest_ascent_rk45(
                #     all_points, dens_func, grad_func, beta_spheres, maximas, tol=1e-7, max_ss=0.5, ss_0=0.23
                # )
                basins = steepest_ascent_rk45(
                    all_points, dens_func, grad_func #, beta_spheres, maximas, tol=1e-7, max_ss=0.5, ss_0=0.23
                )
                final = time.time()
                print("Basins", basins)
                print("Difference ", final - start)

                # Get indices of where they went to a different basin
                basin_switch_ind_local = np.where(basins != i_maxima)[0]
                watershed_indices = indices_to_classify[basin_switch_ind_local]
                print("Global indices (Points) that needs to be refined", watershed_indices)

                # If some points went to a different a basin, then record them as IAS
                indices_to_classify = np.delete(indices_to_classify, basin_switch_ind_local)
                ias_indices += list(watershed_indices)
                ias_basin += list(basins[basin_switch_ind_local])
                ias_radius += [radial[i_rad - 1]] * len(basin_switch_ind_local)


            # Rest of the points are OAS
            oas_indices = indices_to_classify

            # Sort the IAS points
            indices = np.argsort(ias_indices)
            ias_indices = np.array(ias_indices, dtype=int)[indices]
            ias_basin = np.array(ias_basin, dtype=int)[indices]
            ias_radius = np.array(ias_radius)[indices]

            print("IAS indices ", ias_indices)
            print("IAS basins", ias_basin)
            print("OAS indices ", oas_indices)

            """
            #Useful for debugging the classification process
            old_points = maxima + angular_pts
            import matplotlib
            import matplotlib.pyplot as plt
            from mpl_toolkits import mplot3d
            matplotlib.use("Qt5Agg")
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            p = maximas
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
            p = old_points[ias_indices]
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
            p = old_points[oas_indices]
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
            p = old_points[[99, 100]]
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="m", s=40)
            p = np.array([
                           [3.33215211e+00, 3.63210261e+00, -6.14962715e-01],
                           [3.33214688e+00, -3.63213146e+00, 6.14961201e-01]])
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="y", s=50)
            plt.show()
            """

            # Solve for the watershed/IAS points.
            all_ss = np.array([ss0 / 10.0] * len(ias_indices))  # Step-size of each watershed point
            all_ss = np.fmax(all_ss, bnd_err)
            ias_radius = np.array(ias_radius) + all_ss  # increase by ss
            indices_to_solve = ias_indices.copy()
            while len(indices_to_solve) != 0:
                print("Indices to solve watershed ", indices_to_solve)
                print("Current step-size ", all_ss)
                # Construct points on the angular points that aren't classified yet.
                all_points = maxima + ias_radius[:, None] * angular_pts[indices_to_solve, :]

                # If the density values is less than the isosurface value, then the point
                #  should have been classified as a OAS
                dens_vals = dens_func(all_points)
                print("Density Values ", dens_vals)
                dens_small_ind = np.where(dens_vals < iso_val)[0]
                if len(dens_small_ind) != 0:
                    print("(Local) Indices where density is small ", dens_small_ind)
                    print("(Global) indices where density is small ", indices_to_solve[dens_small_ind])
                    print("Before ias_indices ", ias_indices)
                    # Remove (globally) from ias_indices, add to oas_indices
                    is_in = np.isin(ias_indices, indices_to_solve[dens_small_ind])
                    oas_indices = np.hstack((oas_indices, ias_indices[np.where(is_in)[0]]))
                    ias_indices = ias_indices[np.where(~is_in)[0]]
                    ias_basin = ias_basin[np.where(~is_in)[0]]
                    # Delete from local information
                    indices_to_solve = np.delete(indices_to_solve, dens_small_ind)
                    ias_radius = np.delete(ias_radius, dens_small_ind)
                    all_ss = np.delete(all_ss, dens_small_ind)
                    all_points = np.delete(all_points, dens_small_ind, axis=0)

                # Calculate basins of all of these points
                start = time.time()
                basins = steepest_ascent_rk45(
                    all_points, dens_func, grad_func, beta_spheres, maximas, tol=1e-7, max_ss=0.5, ss_0=0.23
                )
                final = time.time()
                print("Basins Assigned ", basins)
                print("Difference ", final - start)
                # Get indices of where they went to a different basin
                basin_switch_ind_local = np.where(basins != i_maxima)[0]
                print("Global indices (Points) that needs to be refined: ", indices_to_solve[basin_switch_ind_local])

                # If basins are different, make sure they match the correct basins
                if len(basin_switch_ind_local) != 0:
                    basins_vals = basins[basin_switch_ind_local]
                    print("Basin_vals that switched ", basins_vals)
                    print("Local Indices that switched ", basin_switch_ind_local)
                    print("IAS indices ", ias_indices)
                    print("IAS basins", ias_basin)
                    print("Actual indices ", np.where(np.in1d(ias_indices, indices_to_solve[basin_switch_ind_local]))[0])
                    original_basins = ias_basin[np.where(np.in1d(ias_indices, indices_to_solve[basin_switch_ind_local]))[0]]
                    print("Original Basins ", original_basins)
                    # if np.any(basins_vals != original_basins):
                    #     raise ValueError(f"Basin switched")

                # Check convergence of watershed points that switch to different basin
                watershed_conv_ind = np.where(all_ss <= bnd_err)[0]
                print(all_ss[basin_switch_ind_local], bnd_err, all_ss[basin_switch_ind_local] <= bnd_err)
                indices_conv = indices_to_solve[watershed_conv_ind]
                if len(indices_conv) != 0:
                    print("Global Indices that converged ", indices_conv)
                    # Get the boundary points:
                    radius_bnd_pts = (2.0 * ias_radius[watershed_conv_ind] + all_ss[watershed_conv_ind]) / 2.0
                    bnd_pts = maxima + radius_bnd_pts[:, None] * angular_pts[watershed_conv_ind]

                    # Store the result:
                    r_func[i_maxima][indices_conv] = np.linalg.norm(bnd_pts - maxima, axis=1)
                    [ias[i_maxima].append(x) for x in indices_conv]
                    original_basins = ias_basin[np.where(np.in1d(ias_indices, indices_to_solve[basin_switch_ind_local]))[0]]
                    [basin_ias[i_maxima].append(x) for x in original_basins]

                    # Delete to avoid for the next iteration
                    indices_to_solve = np.delete(indices_to_solve, watershed_conv_ind)
                    ias_radius = np.delete(ias_radius, watershed_conv_ind)
                    all_ss = np.delete(all_ss, watershed_conv_ind)
                    basins = np.delete(basins, watershed_conv_ind)


                # The ones that different converge, adjust its radius and step-size
                basin_same_ind_local = np.where(basins == i_maxima)[0]
                basin_switch_ind_local = np.where(basins != i_maxima)[0]
                # The ones that basins didn't switch, take a step with step-size
                #   if it reached upper-bound then a problem occured.
                ias_radius[basin_same_ind_local] += all_ss[basin_same_ind_local]
                # TODO: Add a upper-bound check
                # The ones that basins switched, take a step-back and adjust step-size
                #  adjust upper-bound to be the current point.
                ias_radius[basin_switch_ind_local] -= all_ss[basin_switch_ind_local]
                all_ss[basin_switch_ind_local] = np.fmax(all_ss[basin_switch_ind_local] / 10.0, bnd_err)

                print("\n")

            # Solve for the root of each OAS indices
            for i_oas in oas_indices:
                # Construct upper and lower bound of the isosurface equation
                ang_pt = angular_pts[i_oas]
                iso_eq = np.abs(dens_func(maxima + ang_pt * radial[:, None]) - iso_val)
                i_iso = np.argsort(iso_eq)[0]
                l_bnd = radial[i_iso - 1] if i_iso >= 0 else radial[i_iso] / 2.0
                u_bnd = radial[i_iso + 1] if i_iso + 1 < len(radial) else radial[i_iso] * 2.0
                # Solve for the isosurface point
                oas_pt = solve_for_isosurface_pt(
                    l_bnd, u_bnd, maxima, angular_pts[i_oas], dens_func, iso_val, iso_err
                )
                # print("Check isosurface pt", oas_pt, dens_func(np.array([oas_pt])))
                # Record them
                r_func[i_maxima][i_oas] = np.linalg.norm(oas_pt - maxima)
                oas[i_maxima].append(i_oas)

            if type(refine) == type(True) and refine:  # refine can be integer, so this ignores it.
                # Take convex hull between ias and oas and construct additional points in that region.
                #  `new_pts` is concatenated to angular grids and is in cartesian coordinates.
                print("IAS ", ias[i_maxima])
                print("OAS", oas[i_maxima])
                new_pts = construct_points_between_ias_and_oas(
                    ias[i_maxima], oas[i_maxima], angular_pts, r_func[i_maxima], maxima
                )
                print("new pts ", new_pts, np.linalg.norm(new_pts, axis=1))
                # Re-do this qtaim algortihm only on this center
                refined_qtaim = qtaim_surface(new_pts, maximas, dens_func,
                                              grad_func, iso_val, dens_cutoff,
                                              bnd_err, iso_err, beta_spheres=beta_spheres,
                                              optimize_centers=False, refine=i_maxima)
                print("Refined", refined_qtaim.ias, refined_qtaim.oas)
                # Update this basin's result from the refined, + numb_ang_pts: corrects indices
                ias[i_maxima] += [x + numb_ang_pts for x in refined_qtaim.ias[i_maxima]]
                oas[i_maxima] += [x + numb_ang_pts for x in refined_qtaim.oas[i_maxima]]
                basin_ias[i_maxima] += refined_qtaim.basins_ias[i_maxima]
                refined_ang.append(new_pts)
                print(refined_qtaim.r_func, r_func[i_maxima].shape)
                r_func[i_maxima] = np.hstack((r_func[i_maxima], refined_qtaim.r_func[i_maxima]))
                # input("Why")

    print("\n")
    return SurfaceQTAIM(r_func, angular, maximas, oas, ias, basin_ias, refined_ang)


