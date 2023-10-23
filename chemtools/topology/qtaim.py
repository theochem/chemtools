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
import warnings
from collections import OrderedDict
import itertools
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial.transform.rotation import Rotation as r

from grid.angular import AngularGrid

from chemtools.topology.surface import SurfaceQTAIM
from chemtools.topology.utils import (
    construct_radial_grids,
    determine_beta_spheres_and_nna,
    find_optimize_centers,
    solve_for_oas_points,
)
from chemtools.topology.ode import find_basins_steepest_ascent_rk45


__all__ = ["qtaim_surface_vectorize"]


def _classify_rays_as_ias_or_oas(
    maximas, maximas_to_do, all_points, all_basins, index_to_atom, numb_rays_to_atom, numb_rad_to_radial_shell
):
    r"""
    Classify all rays in a molecule as either crossing the outer or inner atomic surface.

    Also provides the interval limits [r_0, r_1] of each ray that crosses the IAS.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        Optimized centers of the electron density where the rays need to be classified. Doesn't need to be
        all of the centers.
    all_points: ndarray(N, 3)
        All points in all rays across each atom in a molecule.
    all_basins: ndarray(N, 3)
        All basin values that were assigned for each point in all rays
    index_to_atom:  list[int]
        Gives the indices that assings each point in `all_points` to each atom, i.e.
        [0, i_1, i_2, ..., N] implies `all_points[0: i_1]` corresponds to atom 1.
    numb_rays_to_atom: list[int]
        List of size `M`, that holds the number of angular points or rays in each atom.
        Used to assign which points in `all_points` corresponds to which ray.
    numb_rad_to_radial_shell: list[list[int]]
        For each atom, for each angular pt/ray, tells the number of radial point.
        Used to assign which points in `all_points` corresponds to which point in each ray.

    Returns
    -------
    ias, oas, ias_bnds, ias_basins: list[list[int]], list[list[int]], list[OrderedDict], list[list[int]]
        A list of size `M`, that holds a list of indices of which angular pt/ray corresponds to
        either intersecting the ias or oas. The final element is a list of size `M`, of
        ordered dictionary whose keys are the indices of ias and items are the lower
        and upper-bound of the radius where the intersection occurs somewhere inbetween.
        `ias_basins` contains which basin each ias pt in `ias` switches to, when it crosses the boundary.

        These are repeat again for searching for second intersection ias_2, ias_basins_2, ias_bnds_2, and
        similarly for the third intersection.

    """
    numb_maximas = len(maximas)
    oas = [[] for _ in range(numb_maximas)]  # outer atomic surface
    ias = [[] for _ in range(numb_maximas)]  # inner atomic surface for first intersection.

    # The points that all converge to the same point are OAS, and are isosurface points. Remove
    #  them
    #     First to get the maxima, you would use index_to_atom.
    #     To get each ray, assuming the number of pts in each ray is the same, you would do
    #            use numb_rad_to_atom
    #     Initially each ray has the same number of points, then one can re-shape both
    #     points and basins to make it easy to index each ray of each maxima, then classify
    #     it as either a IAS or OAS.   Here you can determine whether it crosses twice, and
    #     determine which ray requires special attention.
    #
    # print("Index to atom ", index_to_atom, all_points.shape)
    ias_bnds = [OrderedDict() for _ in range(0, numb_maximas)]              # Keys are Points index
    ias_basins = [OrderedDict() for _ in range(numb_maximas)]               # Keys are Points index


    ias_2 = [[] for _ in range(numb_maximas)]  # inner atomic surface for second intersection.
    ias_3 = [[] for _ in range(numb_maximas)]  # inner atomic surface for third intersection.
    ias_bnds_2 = [OrderedDict() for _ in range(0, numb_maximas)]              # Keys are Points index
    ias_basins_2 = [OrderedDict() for _ in range(numb_maximas)]               # Keys are Points index
    ias_bnds_3 = [OrderedDict() for _ in range(0, numb_maximas)]              # Keys are Points index
    ias_basins_3 = [OrderedDict() for _ in range(numb_maximas)]               # Keys are Points index
    np.set_printoptions(threshold=np.inf)
    for i_do, i_maxima in enumerate(maximas_to_do):
        # print("ATom i ", i_maxima)
        # print("Starting and Final index", index_to_atom[i_maxima], index_to_atom[i_maxima + 1])
        basins_a = all_basins[index_to_atom[i_do]:index_to_atom[i_do + 1]]  # Basin of atom
        points_a = all_points[index_to_atom[i_do]:index_to_atom[i_do + 1]]  # Points of atom
        numb_rad_pts = numb_rad_to_radial_shell[i_do]

        # print("Basins of atom ", basins_a)
        # print(index_to_atom[i_maxima], numb_rad_pts)
        # print("Number of angular points in this atom", numb_rays_to_atom[i_maxima])
        i_ray = 0
        for i_ang in range(numb_rays_to_atom[i_do]):
            # print("Angular pt j", i_ang)
            # print("Number of radial points in this angular pt ", numb_rad_pts[i_ang])

            # Get the basin of the ray
            basins_ray = basins_a[i_ray:i_ray + numb_rad_pts[i_ang]]
            # print("Basin of the ray ", basins_ray)

            # Classify basins as either OAS and IAS, if IAS, then count the number of
            #     intersections of the IAS. In addition, get the l_bnd, u_bnd of each intersection.
            # Groups basins i.e. [1, 1, 1, 0, 0, 0, 1, 1, 2, 2] -> [1, 0, 1, 2]
            group_by = [(k, list(g)) for k, g in itertools.groupby(basins_ray)]
            unique_basins = np.array([x[0] for x in group_by])
            # print(basins_ray == i_maxima)
            # print(group_by)

            # All pts in the ray got assigned to the same basin of the maxima
            #     or if the numb_rad_pts is zero, then this means that the ray has density value less than
            #     the isosurface value (starting from the beta-sphere determination), thus it is an oas point.
            if (len(unique_basins) == 1 and unique_basins[0] == i_maxima) or numb_rad_pts[i_ang] == 0:
                # This implies it is an OAS point, else then it is an IAS with a bad ray.
                # print("OAS Point")
                oas[i_maxima].append(i_ang)
            else:
                # The point is an IAS, determine the number of intersections.
                conv_to_atom = unique_basins == i_maxima
                numb_intersections = np.sum(conv_to_atom)
                l_bnd_pad = 0.0  # This is the case with a bad ray, loweres the l_bnd by this amount
                if numb_intersections == 0:
                    # This is IAS with a bad ray, would have to re-determine the l_bnd. This was an IAS point
                    # Whose ray at the boundary is assigned to different basins depending on the accuracy of ODE solver.
                    # This is IAS with a bad ray, would have to re-determine the l_bnd
                    l_bnd_pad = 0.1

                if 0 <= numb_intersections:
                    # print("IAS With one Intersection.")
                    # Determine lower and upper-bound Point on ray.
                    if group_by[0][1][0] == i_maxima:
                        # if the ray started with a basin that converged ot i_maxima, then take the upper bound
                        #   to be the when it started to switch to a different basin.
                        index_u_bnd = len(group_by[0][1])
                        index_l_bnd = index_u_bnd - 1
                    else:
                        # Here the ray is a bad ray in the sense that the start of the ray should have converged to
                        #  i_maxima but it dind't, and so take the u_bnd to be when it or sure converges to the
                        #  different maxima from i_maxima.
                        index_u_bnd = min(2, len(group_by[0][1]))
                        index_l_bnd = 0
                        if index_u_bnd == index_l_bnd:
                            raise RuntimeError(f"Algorithm Error .")
                    # Determine radius from the upper and lower bound.
                    r_ubnd = np.linalg.norm(points_a[i_ray + index_u_bnd] - maximas[i_maxima])
                    r_lbnd = np.linalg.norm(points_a[i_ray + index_l_bnd] - maximas[i_maxima])
                    # Update containers
                    ias_bnds[i_maxima][i_ang] = [max(0.1, r_lbnd - l_bnd_pad), r_ubnd]
                    ias[i_maxima].append(i_ang)
                    # Get the basins where it switches
                    ias_basins[i_maxima][i_ang] = unique_basins[np.argmax(unique_basins != i_maxima)]
                    # print("Radius Lower and Upper bound ", r_lbnd, r_ubnd)

                    # Gather information about other intersections
                    if numb_intersections > 1:
                        # print("Angular pt j", i_ang)
                        # print("Number of radial points in this angular pt ", numb_rad_pts[i_ang])
                        # print("IAS With Multiple Intersections")
                        # print(group_by)
                        # print(unique_basins)

                        # Figure out if the number intersections is two or three. Code only checks up to three.
                        index_ray = 0  # Keeps track of how many points to go further
                        more_than_three = False
                        for i, (basin, assigned_basin_vals) in enumerate(group_by):
                            if i != 0:
                                # Multiple intersections found, find correct intervals to search for intersections
                                if basin == i_maxima:
                                    if more_than_three:
                                        print("Angular pt j", i_ang)
                                        print("Number of radial points in this angular pt ", numb_rad_pts[i_ang])
                                        print(f"Unique basins {unique_basins}")
                                        print(f"Group by {group_by}")
                                        # raise RuntimeError(f"More than three intersections was found."
                                        #                    f" Code doesn't check.")
                                        warnings.warn(f"More than three intersections was found."
                                                           f" Code doesn't check.")

                                    # Add the second intersection
                                    ias_2[i_maxima].append(i_ang)
                                    ias_basins_2[i_maxima][i_ang] = group_by[i - 1][0]
                                    l_bnd = points_a[i_ray + index_ray - 1]
                                    u_bnd = points_a[i_ray + index_ray]
                                    r_ubnd = np.linalg.norm(u_bnd - maximas[i_maxima])
                                    r_lbnd = np.linalg.norm(l_bnd - maximas[i_maxima])
                                    ias_bnds_2[i_maxima][i_ang] = [r_lbnd, r_ubnd]

                                    # Check for the third intersection that occurs afterwards
                                    if i + 1 < len(group_by):
                                        ias_3[i_maxima].append(i_ang)
                                        ias_basins_3[i_maxima][i_ang] = group_by[i + 1][0]
                                        i_lbnd = i_ray + index_ray + len(assigned_basin_vals) - 1
                                        i_ubnd = i_ray + index_ray + len(assigned_basin_vals)
                                        l_bnd = points_a[i_lbnd]
                                        u_bnd = points_a[i_ubnd]
                                        r_ubnd = np.linalg.norm(u_bnd - maximas[i_maxima])
                                        r_lbnd = np.linalg.norm(l_bnd - maximas[i_maxima])
                                        ias_bnds_3[i_maxima][i_ang] = [r_lbnd, r_ubnd]
                                        more_than_three = True
                            index_ray += len(assigned_basin_vals)

                        # print(ias_2[i_maxima])
                        # print(ias_basins_2[i_maxima])
                        # print(ias_bnds_2[i_maxima])
                        # print(ias_3[i_maxima])
                        # print(ias_basins_3[i_maxima])
                        # print(ias_bnds_3[i_maxima])
                        # import matplotlib
                        # import matplotlib.pyplot as plt
                        # from mpl_toolkits import mplot3d
                        # matplotlib.use("Qt5Agg")
                        # fig = plt.figure()
                        # ax = plt.axes(projection='3d')
                        # p = points_a[i_ray: i_ray + numb_rad_pts[i_ang]]
                        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
                        # ax.scatter(maximas[:, 0], maximas[:, 1], maximas[:, 2], color="r", s=60)
                        # plt.show()
                    else:
                        # Store [l_bnd, u_bnd] inside ias_bnds_2[i_maxima] for searching to the IAS in case
                        #  the user wants to search for multiple intersections.
                        ias_bnds_2[i_maxima][i_ang] = [[r_ubnd, r_ubnd]]

            i_ray += numb_rad_pts[i_ang]
    return ias, oas, ias_bnds, ias_basins, ias_2, ias_bnds_2, ias_basins_2, ias_3, ias_bnds_3, ias_basins_3


def construct_all_points_of_rays_of_atoms(
    maximas, angular_pts, radial_grid, maximas_to_do, dens_func, iso_val
):
    r"""
    Construct all points of all rays across all molecules.
    """
    #  Need a way to track which points correspond to which maxima,
    #  Need a way to track which sets of points correspond to a ray
    #  index_to_atom = [0, i_1, i_1 + i_2, ..., \sum_j^M i_j] first index always zero and last
    #      always the number of points.
    index_to_atom = [0] * (len(maximas_to_do) + 1)  # First index is always zero
    NUMB_RAYS_TO_ATOM = [len(angular_pts[i]) for i in maximas_to_do]
    numb_rad_to_radial_shell = []  # List of List: Number of radius points per ray
    points = []
    print(NUMB_RAYS_TO_ATOM)
    for i_do, i in enumerate(maximas_to_do):  #range(0, numb_maximas):
        # Construct all points on the atomic grid around atom i
        radial_shells = np.einsum("i,jk->jik", radial_grid[i], angular_pts[i])
        print("Number of radial points", len(radial_grid[i]))
        rs = maximas[i, None, None, :] + radial_shells   # has shape (K, N, 3)
        rs = np.reshape(rs, (rs.shape[0] * rs.shape[1], 3)) # has shape (KN, 3)
        print("Total number of points ", rs.shape)

        # Record information what indices it corresponds to
        numb_rad_to_radial_shell.append([len(radial_grid[i])] * NUMB_RAYS_TO_ATOM[i_do])

        # First remove the density values that are less than isosurface values.
        density_vals = dens_func(rs)
        indices = np.where(density_vals < iso_val)[0]
        if len(indices) != 0:
            rs = np.delete(rs, indices, axis=0)
            # Convert from index I to (i) where i is the angular index and j is the radial.
            for k in indices:
                numb_rad_to_radial_shell[i_do][k // len(radial_grid[i])] -= 1
        
        index_to_atom[i_do + 1] = index_to_atom[i_do] + rs.shape[0]  # Add what index it is
        points.append(rs)
    points = np.vstack(points)  # has shape (Product_{i=1}^M K_i N_i, 3)
    print("Total number of points Over All Molecules ", points.shape)
    return points, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell


def _solve_intersection_of_ias_interval(
    maximas, ias_indices, angular_pts, dens_func, grad_func, beta_spheres, bnd_err, ss_0, max_ss, tol, ias_lengths
):
    r"""
    Solves the intersection of the ray to the inner-atomic surface.

    A radial grid is constructed over each ray based on `ias_indices`.  The basin value
    is assigned to each point, and the point where it swtiches basins is recorded.
     The process is further repeated with a smaller step-size until the distance between two
     points on the ray is less than `bnd_err`.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        Optimized centers of the electron density.
    ias_indices: ndarray(N, 5)
       Rows correspond to each ray that intersects the IAS.
       First index is which index of maxima it originates from, then second
       index is index of angular point/ray, third index is the lower bound radius and fourth
       index is the upper-bound radius, fifth index is the step-size.
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
    r_func = [np.zeros((len(angular_pts[i]),), dtype=np.float64) for i in range(len(maximas))]
    basin_ias = [[-1] * ias_lengths[i] for i in range(len(maximas))]
    while len(ias_indices) != 0:
        # Construct New Points
        points = []
        numb_pts_per_ray = []
        for (i_maxima, i_ang, l_bnd, u_bnd, ss, _, _) in ias_indices:
            # print((i_maxima, i_ang, l_bnd, u_bnd, ss))
            ray = (
                maximas[int(i_maxima)] +
                np.arange(l_bnd, u_bnd + ss, ss)[:, None] * angular_pts[int(i_maxima)][int(i_ang), :]
            )
            points.append(ray)
            numb_pts_per_ray.append(len(ray))
        points = np.vstack(points)

        # Solve for basins
        basins, _ = find_basins_steepest_ascent_rk45(
            points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0
        )
        # print("Basins", basins)

        # Refine the rays further
        index_basins = 0  # Index to iterate through basins
        converge_indices = []
        # print("Average step-size", np.mean(ias_indices[:, -1]))
        make_ode_solver_more_accurate = False
        for i, (i_maxima, i_ang, l_bnd, u_bnd, ss, basin_switch, i_ias) in enumerate(ias_indices):
            basins_ray = basins[index_basins:index_basins + numb_pts_per_ray[i]]
            # print((i_maxima, i_ang, l_bnd, u_bnd, ss))
            # print("Basins ", i, basins_ray)

            # Basins that switch index
            i_switch = np.argmax(basins_ray != i_maxima)
            # print("Index of switch ", i_switch)
            if i_switch == 0:
                print("Basins with bad ray: ", basins_ray, (i_maxima, i_ang, l_bnd, u_bnd, ss))
                # raise ValueError(f"This ray lost it's ability to be an IAS point, as all points converged to the same maxima. "
                #                  f"Fix this.")

                # Update ias_indices for the next iteration for example l_bnd, u_bnd, step-size
                # This lower and upper bound is chosen to guarantee that the IAS point will be found.
                new_l_bnd = l_bnd - 20.0 * ss
                new_u_bnd = u_bnd + 20.0 * ss
                new_ss = max(ss / 10.0, bnd_err)
                # Make ODE solver more accurate
                make_ode_solver_more_accurate = True
                # print("New step-size ", new_ss)
                ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, new_ss, basin_switch, i_ias]
            else:
                # If ss was less than bnd_err, then we converge and should stop.
                if ss <= bnd_err:
                    # Take midpoint to be the radius of intersection
                    radius_mid_pt = (2.0 * l_bnd + ss * i_switch) / 2.0
                    r_func[int(i_maxima)][int(i_ang)] = radius_mid_pt
                    basin_ias[int(i_maxima)][int(i_ias)] = basins_ray[i_switch]
                    converge_indices.append(i)  # Put in list to remove indices.
                else:
                    # Update ias_indices for the next iteration for example l_bnd, u_bnd, step-size
                    new_l_bnd = l_bnd + ss * (i_switch - 1)
                    new_u_bnd = l_bnd + ss * (i_switch)
                    new_ss = max(ss / 10.0, bnd_err)
                    # print("New step-size ", new_ss)
                    ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, new_ss, basin_switch, i_ias]

            # Update index for the next ray
            index_basins += numb_pts_per_ray[i]
            # print("COnvergence indices", converge_indices)
        if make_ode_solver_more_accurate:
            tol /= 2.0
            max_ss = min(0.1, max_ss)
            ss_0 /= 2.0

        # Remove converged indices
        ias_indices = np.delete(ias_indices, converge_indices, axis=0)

    # Solve for multiple intersections
    return r_func, basin_ias


def _solve_intersection_of_ias_point(
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
        # Take the initial step-size of ODE solver based on step-size of ss_0 and the interval of [l_bnd, u_bnd]
        step_size = np.minimum(ias_indices[:, 4][:, None], ss_0)
        basins, _ = find_basins_steepest_ascent_rk45(
            points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=step_size, hess_func=hess_func
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
                # decrease step-size by two to make it smaller than the interval bound
                new_step_size = np.minimum(new_ss / 2.0, ss_0)
                ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, new_step_size, basin_switch, i_ias]

        # Remove converged indices
        # print("Convergence indices", converge_indices)
        ias_indices = np.delete(ias_indices, converge_indices, axis=0)

    # Solve for multiple intersections
    return r_func, basin_ias


def qtaim_surface_vectorize(
    angular,
    centers,
    dens_func,
    grad_func,
    iso_val=0.001,
    bnd_err=1e-4,
    iso_err=1e-5,
    beta_spheres=None,
    beta_sphere_deg=27,
    ss_0=0.1,
    max_ss=0.25,
    tol=1e-7,
    optimize_centers=True,
    hess_func=None,
    find_multiple_intersections=False,
    maximas_to_do=None,
    min_pt_radial=0.1,
    padding_radial=3.0,
    ss_radial=0.24,
):
    r"""
    Parameters
    ----------
    angular: List[AngularGrid] or List[int]
        List of angular grids over each atom, or a list of their degrees.
    centers: ndarray(M, 3)
        Atomic coordinates.
    dens_func: callable(ndarray(N, 3)->ndarray(N,))
        The electron density function.
    grad_func: callable(ndarray(N, 3)->ndarray(N,3))
        The gradient of the electron density.
    iso_val: float
        Isosurface value of the outer-atomic surface.
    bnd_err: float
        The error of the points on the inner-atomic surface.
    iso_err: float
        The error in solving for the isosurface points on the outer-atomic surface.
    beta_spheres: (List[float] or None)
        The radius of confidence that points are assigned to the atom. Should have length `M`.
    beta_sphere_deg: int
        Integer specifying angular grid of degree `beta_sphere_deg` that is used to find the beta-sphere
        automatically, if `beta_spheres` isn't provided. Default value is 21.
    ss_0: float, optional
        The initial step-size of the ODE (RK45) solver.
    max_ss: float, optional
        Maximum step-size of the ODE (RK45) solver.
    tol: float, optional
        Tolerance for the adaptive step-size of the ODE solver.
    optimize_centers: bool
        If true, then the steepest-ascent is performed on the centers to find the local maximas.
    hess_func: callable(ndarray(N, 3)->ndarray(N, 3, 3))
        The Hessian of the electron density. If this is provided, then non-nuclear attractors will be found.
        Adds a default Lebedev/angular grid of degree fifty and 0.1 a.u. to the `beta-spheres` if it is
        provided.
    find_multiple_intersections: bool
        If true, then it searches for up to three intersections of the inter-atomic surface. This is a
        time-consuming process but produces more accurate surfaces.
    maximas_to_do: (None, list[int])
        List of indices of the `centers`/`maximas` to solve for the QTAIM basin surface.  If this is provided,
        then `angular` should also be of this length.
    min_pt_radial: float
        The minimum point to each radial grids that are constructed over each atom. Sometimes, two
        maximas are very close to one another and so it controls how the beta-spheres are determined. 
    padding_radial: float
        Adds a padding to the maximum of each radial grids that are constructed over each atom.
        The default maximum is taken based on the maximum distance between the closest five atoms.
        It affects the decision that the ray crosses the IAS or OAS.
    ss_radial: float, list[float]
        The step-size of the radial grids over each atom. It affects the decision that the ray crosses the IAS or OAS.
        Smaller step-size is able to capture more points that should be on the IAS.

    Returns
    -------
    SurfaceQTAIM:
        Object that holds all information regarding the surface of each atom.

    Notes
    -----
    The algorithm is as follows:
    1. Optimize the centers provided to obtain the local maximas.
    2. Determine the beta-spheres over all atoms.
    3. Using an angular grid and radial grid, construct all rays propagating across all atoms.
    4. Solve for each basin value for each point.
    5. Analyze the basin values and classify each ray as either an outer-atomic or inner-atomic
        surface point.
    6. For the inner-atomic rays, find the point of intersection to the surface boundary based on moving
        along the ray.
    7. For the points on the outer-atomic surface, solve for when the point has density equal to the isosurface value.
    8. Analyze points on the IAS, and if their density values are less than the outer-atomic surface, then
        re-classify them as a OAS point.

    """
    if len(angular) != len(centers):
        raise ValueError(f"Length of angular {len(angular)} should be the same as the"
                         f"number of centers {len(centers)}.")
    if beta_spheres is not None and len(centers) != len(beta_spheres):
        raise ValueError(
            f"Beta sphere length {len(beta_spheres)} should match the"
            f" number of centers {len(centers)}"
        )
    if not isinstance(beta_spheres, (type(None), np.ndarray)):
        raise TypeError(f"Beta_sphers {type(beta_spheres)} should be of numpy type.")
    if maximas_to_do is not None and not isinstance(maximas_to_do, list):
        raise TypeError(f"Maximas to do {type(maximas_to_do)} should be either None or a list of integers.")
    if maximas_to_do is not None and max(maximas_to_do) >= len(centers):
        raise ValueError(f"Length of maximas_to_do {len(maximas_to_do)} should be less then"
                         f" length of centers {len(centers)}.")
    if isinstance(ss_radial, float):
        ss_radial = [ss_radial] * len(centers)
    if len(ss_radial) != len(centers):
        raise ValueError(f"The step-size for each radial grid {len(ss_radial)} should be"
                         f" equal to the number of molecules {len(centers)}.")
    if maximas_to_do is None:
        maximas_to_do = np.arange(len(centers))

    # Using centers, update to the maximas
    maximas = centers
    if optimize_centers:
        # Using ODE solver to refine the maximas further.
        maximas = find_optimize_centers(centers, grad_func)

    # Construct a dense radial grid for each atom by taking distance to the closest five atoms.
    ss0 = 0.1
    # radial_grids = construct_radial_grids(maximas[maximas_to_do], maximas, min_pts=0.1, pad=padding_radial, ss0=ss0)
    radial_grids = construct_radial_grids(maximas, maximas, min_pts=0.1, pad=padding_radial, ss0=ss0)

    # Determine beta-spheres and non-nuclear attractors from a smaller angular grid
    #  Degree can't be too small or else the beta-radius is too large and IAS point got classified
    #  as OAS point.
    ang_grid = AngularGrid(degree=beta_sphere_deg, use_spherical=True)
    if beta_spheres is None:
        beta_spheres, maximas, radial_grids = determine_beta_spheres_and_nna(
            beta_spheres, maximas, radial_grids, ang_grid.points, dens_func, grad_func, hess_func
        )
        beta_spheres = np.array(beta_spheres)
        print(f"Final Beta-spheres {beta_spheres}")
    # Check beta-spheres are not intersecting
    dist_maxs = cdist(maximas, maximas)
    condition = dist_maxs <= beta_spheres[:, None] + beta_spheres
    condition[range(len(maximas)), range(len(maximas))] = False  # Diagonal always true
    if np.any(condition):
        raise ValueError(f"Beta-spheres {beta_spheres} overlap with one another.")
    # Rotate the beta-sphere and double check the radius is correct
    rot_mat = r.from_euler("z", 10, degrees=True).as_matrix()  # Generate rotation by 10 degrees in z-axis
    ang_pts = ang_grid.points.dot(rot_mat.T)                   # Rotate the angular points
    assert np.all(np.abs(np.linalg.norm(ang_pts, axis=1) - 1.0) < 1e-8)
    for i_maxima, radius in enumerate(beta_spheres):
        pts_at_rad = maximas[i_maxima] + radius * ang_pts
        basins, _ = find_basins_steepest_ascent_rk45(
            pts_at_rad, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0
        )
        #assert np.all(basins == i_maxima)
        if not np.all(basins == i_maxima):
            # Decrease the beta-sphere by the step-size if all of the beta-sphere didn't converge to the
            #   correct atom.
            beta_spheres[i_maxima] -= ss0

    # Construct a coarse radial grid for each atom starting at the beta-spheres.
    radial_grids_old = radial_grids
    radial_grids = []
    i_do = 0
    for i_atom in range(len(maximas)):
        if i_atom in maximas_to_do:
            radial_grids.append(
                np.arange(beta_spheres[i_atom], radial_grids_old[i_do][-1], ss_radial[i_atom])
            )
            i_do += 1
        else:
            radial_grids.append([])

    # Construct Angular Points
    angular_pts = []
    for i in range(len(maximas)):
        # If it is not provided, then use what's specified
        if i < len(angular):
            if i in maximas_to_do:
                ang = angular[i]
                if isinstance(ang, int):
                    angular_pts.append(AngularGrid(degree=ang, use_spherical=True).points)
                else:
                    angular_pts.append(ang.points)
            else:
                angular_pts.append([])
        else:
            # If it is a Non-nuclear attractor
            angular.append(99)
            angular_pts.append(AngularGrid(degree=99, use_spherical=True).points)

    # First step is to construct a grid that encloses all radial shells across all atoms
    points, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell = \
        construct_all_points_of_rays_of_atoms(
            maximas, angular_pts, radial_grids, maximas_to_do, dens_func, iso_val
    )
    # print("Index to atom ", index_to_atom)

    # Then assign basins values for all the points.
    import time
    start = time.time()
    basins, _ = find_basins_steepest_ascent_rk45(
        points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0,
        hess_func=hess_func, check_for_nna=True
    )
    final = time.time()
    # print("Basins", basins)
    # print("Length of basins ", len(basins))
    print("Difference ", final - start)

    # Using basin values, classify each ray as either IAS or OAS and if IAS, find the interval
    #   along the ray that intersects the IAS.
    ias, oas, ias_bnds, ias_basins, ias_2, ias_bnds_2, ias_basins_2, ias_3, ias_bnds_3, ias_basins_3 = \
        _classify_rays_as_ias_or_oas(
            maximas, maximas_to_do, points, basins, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell
    )

    print("Total number of two intersections found ", [len(x) for x in ias_2])
    print("Total number of three intersections found ", [len(x) for x in ias_3])

    # The IAS is just refining the ray, till you find the exact intersection with the surface.
    # Checks docs of `_solve_intersection_of_ias` for what ias_indices is.
    # i_maxima, i_ang, l_bnd, u_bnd, ss, basin_switch, i_ias
    ias_indices = np.array(list(
        itertools.chain.from_iterable(
            [[(i, i_ang, ias_bnds[i][i_ang][0], ias_bnds[i][i_ang][1], max(ss_radial[i] / 10.0, bnd_err), ias_basins[i][i_ang], i_ias)
              for i_ias, i_ang in enumerate(ias[i])] for i in maximas_to_do]
        )
    ))
    start = time.time()
    r_func, basin_ias = _solve_intersection_of_ias_point(
        maximas, ias_indices, angular_pts, dens_func, grad_func, beta_spheres, bnd_err,
         tol=tol, max_ss=max_ss, ss_0=ss_0, ias_lengths=[len(x) for x in ias], hess_func=hess_func
    )
    final = time.time()
    print("Time Difference for Solving IAS ", final - start)

    # Pts on IAS that should be OAS are solved here. TODO:, the _solve_intersection should handle it
    #  so it doesn't solve it again.
    for i_atom in maximas_to_do:
        pts = maximas[i_atom] + r_func[i_atom][ias[i_atom], None] * angular_pts[i_atom][ias[i_atom], :]
        dens_vals = dens_func(pts)
        # Decrease by the OAS surface error "iso_err"
        ias_indices = np.where(dens_vals - iso_err < iso_val)[0]
        if len(ias_indices) != 0:
            print(f"Atom {i_atom} points that are IAS should be OAS")
            oas[i_atom] += list(np.array(ias[i_atom], dtype=int)[ias_indices])
            oas[i_atom] = sorted(oas[i_atom])
            ias[i_atom] = [k for j, k in enumerate(ias[i_atom]) if j not in ias_indices]
            basin_ias[i_atom] = [k for j, k in enumerate(basin_ias[i_atom]) if j not in ias_indices]

    # Solve OAS Points and updates r_func
    # Update the radial grid step-size so that it samples more points, this shouldn't decrease computational complexity
    #  since the electron density is cheaper to compute with.
    start = time.time()
    solve_for_oas_points(maximas, maximas_to_do, oas, angular_pts, dens_func, iso_val, iso_err, r_func)
    final = time.time()
    print("Time Difference for Solving OAS", final - start)

    if find_multiple_intersections:
        raise NotImplementedError(f"Multiple intersections was not implemented yet.")

    return SurfaceQTAIM(r_func, angular, maximas, maximas_to_do, oas, ias, basin_ias, iso_val, beta_spheres)




def _solve_intersection_of_ias_point_hirshfeld(
    maximas, ias_indices, angular_pts, basin_assignment_callable, bnd_err, ias_lengths,
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
        # Take the initial step-size of ODE solver based on step-size of ss_0 and the interval of [l_bnd, u_bnd]
        basins = basin_assignment_callable(points)
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
                # decrease step-size by two to make it smaller than the interval bound
                ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, None, basin_switch, i_ias]

        # Remove converged indices
        # print("Convergence indices", converge_indices)
        ias_indices = np.delete(ias_indices, converge_indices, axis=0)

    # Solve for multiple intersections
    return r_func, basin_ias


def hirshfeld_surface_vectorize(
    angular,
    centers,
    dens_func,
    basin_assignment_callable,
    iso_val=0.001,
    bnd_err=1e-4,
    iso_err=1e-5,
    find_multiple_intersections=False,
    maximas_to_do=None,
    min_pt_radial=0.1,
    padding_radial=3.0,
    ss_radial=0.24,
):
    if len(angular) != len(centers):
        raise ValueError(f"Length of angular {len(angular)} should be the same as the"
                         f"number of centers {len(centers)}.")
    if maximas_to_do is not None and not isinstance(maximas_to_do, list):
        raise TypeError(f"Maximas to do {type(maximas_to_do)} should be either None or a list of integers.")
    if maximas_to_do is not None and max(maximas_to_do) >= len(centers):
        raise ValueError(f"Length of maximas_to_do {len(maximas_to_do)} should be less then"
                         f" length of centers {len(centers)}.")
    if isinstance(ss_radial, float):
        ss_radial = [ss_radial] * len(centers)
    if len(ss_radial) != len(centers):
        raise ValueError(f"The step-size for each radial grid {len(ss_radial)} should be"
                         f" equal to the number of molecules {len(centers)}.")
    if maximas_to_do is None:
        maximas_to_do = np.arange(len(centers))

    # Using centers, update to the maximas
    maximas = centers

    # Construct a dense radial grid for each atom by taking distance to the closest five atoms.
    ss0 = 0.1
    radial_grids = construct_radial_grids(maximas, maximas, min_pts=0.1, pad=padding_radial, ss0=ss0)

    # Construct a coarse radial grid for each atom starting at the beta-spheres.
    radial_grids_old = radial_grids
    radial_grids = []
    i_do = 0
    for i_atom in range(len(maximas)):
        if i_atom in maximas_to_do:
            radial_grids.append(
                np.arange(min_pt_radial, radial_grids_old[i_do][-1], ss_radial[i_atom])
            )
            i_do += 1
        else:
            radial_grids.append([])

    # Construct Angular Points
    angular_pts = []
    for i in range(len(maximas)):
        if i in maximas_to_do:
            ang = angular[i]
            if isinstance(ang, int):
                angular_pts.append(AngularGrid(degree=ang, use_spherical=True).points)
            else:
                angular_pts.append(ang.points)
        else:
            angular_pts.append([])

    # First step is to construct a grid that encloses all radial shells across all atoms
    points, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell = \
        construct_all_points_of_rays_of_atoms(
            maximas, angular_pts, radial_grids, maximas_to_do, dens_func, iso_val
    )

    # Then assign basins values for all the points.
    import time
    start = time.time()
    basins = basin_assignment_callable(points)
    final = time.time()
    print("Difference ", final - start)

    # Using basin values, classify each ray as either IAS or OAS and if IAS, find the interval
    #   along the ray that intersects the IAS.
    ias, oas, ias_bnds, ias_basins, ias_2, ias_bnds_2, ias_basins_2, ias_3, ias_bnds_3, ias_basins_3 = \
        _classify_rays_as_ias_or_oas(
            maximas, maximas_to_do, points, basins, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell
    )

    print("Total number of two intersections found ", [len(x) for x in ias_2])
    print("Total number of three intersections found ", [len(x) for x in ias_3])

    # The IAS is just refining the ray, till you find the exact intersection with the surface.
    # Checks docs of `_solve_intersection_of_ias` for what ias_indices is.
    # i_maxima, i_ang, l_bnd, u_bnd, ss, basin_switch, i_ias
    ias_indices = np.array(list(
        itertools.chain.from_iterable(
            [[(i, i_ang, ias_bnds[i][i_ang][0], ias_bnds[i][i_ang][1], max(ss_radial[i] / 10.0, bnd_err), ias_basins[i][i_ang], i_ias)
              for i_ias, i_ang in enumerate(ias[i])] for i in maximas_to_do]
        )
    ))
    start = time.time()
    r_func, basin_ias = _solve_intersection_of_ias_point_hirshfeld(
        maximas, ias_indices, angular_pts, basin_assignment_callable, bnd_err, ias_lengths=[len(x) for x in ias]
    )
    final = time.time()
    print("Time Difference for Solving IAS ", final - start)

    # Pts on IAS that should be OAS are solved here. TODO:, the _solve_intersection should handle it
    #  so it doesn't solve it again.
    for i_atom in maximas_to_do:
        pts = maximas[i_atom] + r_func[i_atom][ias[i_atom], None] * angular_pts[i_atom][ias[i_atom], :]
        dens_vals = dens_func(pts)
        # Decrease by the OAS surface error "iso_err"
        ias_indices = np.where(dens_vals - iso_err < iso_val)[0]
        if len(ias_indices) != 0:
            print(f"Atom {i_atom} points that are IAS should be OAS")
            oas[i_atom] += list(np.array(ias[i_atom], dtype=int)[ias_indices])
            oas[i_atom] = sorted(oas[i_atom])
            ias[i_atom] = [k for j, k in enumerate(ias[i_atom]) if j not in ias_indices]
            basin_ias[i_atom] = [k for j, k in enumerate(basin_ias[i_atom]) if j not in ias_indices]

    # Solve OAS Points and updates r_func
    # Update the radial grid step-size so that it samples more points, this shouldn't decrease computational complexity
    #  since the electron density is cheaper to compute with.
    start = time.time()
    solve_for_oas_points(maximas, maximas_to_do, oas, angular_pts, dens_func, iso_val, iso_err, r_func)
    final = time.time()
    print("Time Difference for Solving OAS", final - start)

    if find_multiple_intersections:
        raise NotImplementedError(f"Multiple intersections was not implemented yet.")

    return SurfaceQTAIM(r_func, angular, maximas, maximas_to_do, oas, ias, basin_ias, iso_val, None)
