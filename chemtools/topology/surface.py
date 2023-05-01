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
r"""
Data structure that holds the outer-atomic (OAS) and intra-atomic (IAS) surfaces.

Can be used for
- analyzing the IAS and OAS.
- integration over basins.
"""
from chemtools.topology.utils import solve_for_oas_points, solve_intersection_of_ias_point

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

from grid.angular import AngularGrid
from grid.atomgrid import AtomGrid
from grid.basegrid import Grid

__all__ = ["SurfaceQTAIM"]


class SurfaceQTAIM:
    def __init__(self, r_func, angular_degs, maximas, indices_maxima, oas, ias, basins_ias, iso_val, beta_spheres):
        self._r_func = r_func
        self._maximas = maximas
        self._indices_maxima = indices_maxima
        self._angular_degs = angular_degs
        self._oas = oas
        self._ias = ias
        self._basins_ias = basins_ias
        self._iso_val = iso_val
        self._beta_spheres = beta_spheres

    @property
    def beta_spheres(self):
        # List[Float]: Radius points that indicate all points converge to the maxima.
        return self._beta_spheres
    @property
    def iso_val(self):
        # Float: Isosurface value of the outer atomic surface
        return self._iso_val

    @property
    def r_func(self):
        # List[M, np.ndarray(N_i,)] ::math:`r_j(\theta_i, \phi_i)` for each M basins and N_i
        # angular points. Length matches length of `maximas`.
        return self._r_func

    @property
    def oas(self):
        # List[List[int]] : First list is over basins, second over indices of points of outeratomic
        # surface. Length matches length of `indices_maxima`
        return self._oas

    @property
    def ias(self):
        # List[List[int]] : First list is over basins, second over indices of points of interatomic
        # surface. Length matches length of `indices_maxima`
        return self._ias

    @property
    def maximas(self):
        # ndarray(M, 3) : The maxima of all basins.
        return self._maximas

    @property
    def indices_maxima(self):
        # list[int]: List of indices of `centers` that correspond to each row of `maximas`.
        return self._indices_maxima

    @indices_maxima.setter
    def indices_maxima(self, value):
        self._indices_maxima = value

    @property
    def angular_degs(self):
        # int or List[int] : List of angular degrees over each basin. Length matches length of `indices_maxima`
        return self._angular_degs

    @property
    def rad_grids(self):
        # List[OneDGrids] : List of M OneDGrids for integrating over radial component in [0, \inty).
        return self._rad_grids

    @property
    def basins_ias(self):
        return self._basins_ias

    def __add__(self, other):
        if np.abs(self.iso_val - other.iso_val) > 1e-8:
            raise RuntimeError(f"Isosurface value should match each other when combing surface objects.")
        object = SurfaceQTAIM(
            self.r_func, self.angular_degs, self.maximas, self.indices_maxima, self.oas, self.ias,
            self.basins_ias, self.iso_val, self.beta_spheres
        )
        object.indices_maxima = sorted(list(set(self.indices_maxima).union(set(other.indices_maxima))))
        for i_replace in other.indices_maxima:
            object.ias[i_replace] = other.ias[i_replace]
            object.oas[i_replace] = other.oas[i_replace]
            object.basins_ias[i_replace] = other.basins_ias[i_replace]
            object.angular_degs[i_replace] = other.angular_degs[i_replace]
            object.r_func[i_replace] = other.r_func[i_replace]
        return object

    def save(self, filename):
        save_dict = {
            "ias": np.array(self.ias, dtype=np.object),
            "oas": np.array(self.oas, dtype=np.object),
            "basin_ias": np.array(self.basins_ias, dtype=np.object),
            "maximas": np.array(self.maximas),
            "indices_maxima": np.array(self.indices_maxima),
            "angular_degs": np.array(self.angular_degs),
            "r_func": np.array(self.r_func, dtype=np.object),
            "iso_val": self.iso_val
        }
        np.savez(filename + ".npz", **save_dict, allow_pickle=True)

    def generate_angular_grid_of_basin(self, i_basin):
        # Note this doesn't include the extra angular points generated by refinement.
        deg = self.angular_degs
        deg = deg[i_basin] if isinstance(deg, (list, np.ndarray)) else deg
        return AngularGrid(degree=deg, use_spherical=True)

    def generate_angular_pts_of_basin(self, i_basin):
        angular_grid = self.generate_angular_grid_of_basin(i_basin)
        points = angular_grid.points
        return points

    def get_atom_grid_over_basin(self, i_basin, rgrid=None):
        # integrate over a basin.
        deg = self.angular_degs
        deg = deg[i_basin] if isinstance(deg, list) else deg
        rgrid = self.rad_grids if rgrid is None else rgrid
        # If rgrid is a list for every center, then obtain the right one. Else is it is OneDGrid.
        if isinstance(rgrid, list):
            if len(rgrid) > 1:
                rgrid = self.rad_grids[i_basin]
            else:
                rgrid = rgrid[0]
        atom_grid = AtomGrid(rgrid, degrees=[deg], center=self.maximas[i_basin], use_spherical=True)

        ias_indices_a = self.ias[i_basin]
        r_limits = self.r_func[i_basin][ias_indices_a]

        # Holds indices of each point on the angular grid, where the radial points should be zero afterwards.
        ias_indices, rad_indices = np.where(atom_grid.rgrid.points[None, :] > r_limits[:, None])
        start_indices = atom_grid.indices[rad_indices]  # Get the radial shell that includes the index rad_indices.
        ias_indices = np.array(ias_indices, dtype=int)
        start_indices = np.array(start_indices, dtype=int)
        indices_zero = np.array(start_indices + np.array(ias_indices_a, dtype=int)[ias_indices], dtype=int)
        atom_grid.weights[indices_zero] = 0.0

        indices_zero = np.where(atom_grid.weights == 0.0)[0]
        points = np.delete(atom_grid.points, indices_zero, axis=0)
        weights = np.delete(atom_grid.weights, indices_zero)

        return Grid(points, weights)

    def generate_pts_on_surface(self, i_basin, include_other_surfaces=False):
        r"""
        Generates points on the surface of an atom.

        Parameters
        ----------
        i_basin : int
            Which basin you want to generate the surface for.
        include_other_surfaces: bool
            If true, then it add IAS points from other basins (other than `i_basin`)
            that crosses the `i_basin`th basin. This adds more points to the surface.

        Returns
        -------
        ndarray(N, 3):
            3D-coordinates of each point on the surface.

        """
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        if len(self.r_func[i_basin]) != 0:
            points = self.maximas[i_basin] + self.r_func[i_basin][:, None] * sph_pts
        else:
            points = np.empty((0, 3), dtype=np.float64)
        if include_other_surfaces:
            for i in range(len(self.maximas)):
                if i != i_basin:
                    # If this basin crosses the boundary.
                    indices = np.where(np.abs(i_basin - np.array(self.basins_ias[i])) < 1e-10)[0]
                    if len(indices) != 0:
                        ias_indices = np.array(self.ias[i])[indices]
                        sph_pts = self.generate_angular_pts_of_basin(i)
                        new_pts = self.maximas[i] + self.r_func[i][ias_indices, None] * sph_pts[ias_indices]
                        points = np.vstack((points, new_pts))
        # TODO: Probably better to round it and check if it is unique,
        return np.unique(np.round(points, 16), axis=0)

    def get_ias_pts_of_basin(self, i_basin, include_other_surfaces=False):
        ias = self.ias[i_basin]
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        points = self.maximas[i_basin] + self.r_func[i_basin][ias, None] * sph_pts[ias]
        if include_other_surfaces:
            for i in range(len(self.maximas)):
                if i != i_basin and i in self.indices_maxima:
                    # If this basin crosses the boundary.
                    indices = np.where(np.abs(i_basin - np.array(self.basins_ias[i])) < 1e-10)[0]
                    if len(indices) != 0:
                        ias_indices = np.array(self.ias[i])[indices]
                        sph_pts = self.generate_angular_pts_of_basin(i)
                        new_pts = self.maximas[i] + self.r_func[i][ias_indices, None] * sph_pts[ias_indices]
                        points = np.vstack((points, new_pts))
        return np.unique(np.round(points, 16), axis=0)

    def get_oas_pts_of_basin(self, i_basin):
        oas = self.oas[i_basin]
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        return self.maximas[i_basin] + self.r_func[i_basin][oas, None] * sph_pts[oas]

    def construct_points_between_ias_and_oas(self, basin_ids, dens_func, grad_func, ss_0, max_ss, tol, iso_err):
        r"""
        Construct points between the inner atomic surface and outer atomic surface.

        This is done by constructed a convex hull between IAS and OAS, separately.
        Each point on the IAS, the two closest points are found on the OAS, then
        a triangle is constructed.  Seven points are constructed within this triangle
        and the Cartesian coordinates of the sphere centered at the maxima is solved
        for each of these seven points.

        Parameters
        -----------
        basin_ids: list[int]
            List of integers specifying the index of the basins/maximas for refinement is to be performed.
        iso_err: float
            The error in solving for the isosurface points on the outer-atomic surface.

        Returns
        -------
        ndarray(K * 7, 3)
            Cartesian coordinates of :math:`K` points on the sphere centered at `maxima` such that
            they correspond to the seven points constructed above, where :math:`K` is the number
            of points on the IAS of `maxima`.

        """
        ias_parameters = []  # parameters needed for solving IAS
        all_angular_pts = []
        indices_for_each_basin = [0]
        for i_basin in basin_ids:
            maxima = self.maximas[i_basin]
            ias = self.ias[i_basin]
            oas = self.oas[i_basin]
            r_func_max = self.r_func[i_basin]
            angular_pts = self.generate_angular_pts_of_basin(i_basin)
            # Take a convex hull of both IAS and OAS seperately.
            ias_pts = maxima + r_func_max[ias, None] * angular_pts[ias, :]
            oas_pts = maxima + r_func_max[oas, None] * angular_pts[oas, :]
            ias_hull = ConvexHull(ias_pts)
            oas_hull = ConvexHull(oas_pts)
            ias_bnd = ias_hull.points[ias_hull.vertices]
            oas_bnd = oas_hull.points[oas_hull.vertices]

            # Compute the distance matrix
            dist_mat = cdist(ias_bnd, oas_bnd)
            i_smallest_dist = np.argmin(dist_mat, axis=1)
            smallest_dist = dist_mat[np.arange(len(ias_bnd)), i_smallest_dist]
            # print("Smallest distance between ias and oas ", np.sort(smallest_dist))
            # print("new pts ", new_pts)

            # fig = plt.figure()
            # ax = plt.axes(projection='3d')
            # p = ias_bnd
            # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
            # p = ias_bnd[smallest_dist < np.exp(np.mean(np.log(smallest_dist)) - 2.0 * np.std(np.log(smallest_dist))), :]
            # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="m")
            # p = oas_bnd
            # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
            # plt.show()

            # for each point in say ias take the closest two points in oas.
            new_ang_pts = np.zeros((0, 3), dtype=np.float64)  # usually 7 points per ias boundary are added.
            # Assumes the distances have a log-normal distribution, and taking the second quantile to get the
            #   points closest to the OAS from the IAS.
            extreme_ends = np.exp(np.mean(np.log(smallest_dist)) - 1.2 * np.std(np.log(smallest_dist)))
            # print("Distance tolerance ", extreme_ends)
            indices = np.where(smallest_dist < extreme_ends)[0]
            i_angular = 0
            for i_ias in indices:
                pt_ias = ias_bnd[i_ias, :]
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
                new_pts = new_pts - maxima
                # Delete points that are on the maxima.
                direction = np.delete(new_pts, np.all(np.abs(new_pts) < 1e-10, axis=1), axis=0)
                radiuses = np.linalg.norm(direction, axis=1)
                direction = direction / radiuses[:, None]
                # Delete directions that are the same
                direction, indices = np.unique(direction, axis=0, return_index=True)
                radiuses = radiuses[indices]
                new_ang_pts = np.vstack((new_ang_pts, direction))

                # Add the correct IAS parameters to `ias_parameters`
                for i in range(len(new_pts)):
                    # Construct lower and upper bound such that it is at the midpoint
                    l_bnd = radiuses[i] - 0.1
                    u_bnd = radiuses[i] + 0.1
                    ss = 0.05
                    ias_parameters.append(
                        [i_basin, i_angular, l_bnd, u_bnd, ss, -1, i_angular]
                    )
                    i_angular += 1

            all_angular_pts.append(new_ang_pts)
            indices_for_each_basin.append(sum(indices_for_each_basin) + len(new_ang_pts))
        indices_for_each_basin.append(len(all_angular_pts))
        ias_parameters = np.array(ias_parameters)

        print("Solve for the new refinement")
        # Solve for the IAS

        angular_pts = [[0.0, 0.0, 0.0]] * len(self.maximas)
        ias_lengths = [1] * len(self.maximas)
        for i, i_basin in enumerate(basin_ids):
            angular_pts[i_basin] = all_angular_pts[i]
            ias_lengths[i_basin] = len(all_angular_pts[i])
        r_func_new, _ = solve_intersection_of_ias_point(
            self.maximas, ias_parameters, angular_pts, dens_func, grad_func, self.beta_spheres,
            bnd_err=1e-5, ias_lengths=ias_lengths, ss_0=ss_0, max_ss=max_ss, tol=tol,
        )

        # For each basin_id, check if their density values are not less than the isosurface value.
        all_pts = []
        for i, i_basin in enumerate(basin_ids):
            new_pts = self.maximas[i_basin] + r_func_new[i_basin][:, None] * all_angular_pts[i]
            print(f"Number of new points to add: {len(new_pts)}")

            # Check if the new points are less than isosurface value and project them so that they do have .
            dens_vals = dens_func(new_pts)
            indices = np.where(dens_vals < self.iso_val)[0]
            if len(indices) != 0:
                # Construct bounded interval to solve for the root.
                solve_for_oas_points(
                    np.array([self.maximas[i_basin]]), [0], [indices], [all_angular_pts[i]],
                    dens_func, self.iso_val, iso_err, [r_func_new[i_basin]]
                )
                new_pts[indices] = self.maximas[i_basin] + r_func_new[i_basin][indices, None] * all_angular_pts[i][indices, :]

            all_pts.append(new_pts)
        return all_pts

    def plot_basins(self, i_basin, include_other_surfaces=False):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        p = self.maximas
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60, label="Maximas")
        p = self.get_ias_pts_of_basin(i_basin, include_other_surfaces=include_other_surfaces)
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
        p = self.get_oas_pts_of_basin(i_basin)
        ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
        plt.show()

    def interpolate_radial_func(self, method="smooth", ias=False, oas=False):
        # if method not in ["smooth", ]
        if ias and oas:
           raise ValueError(f"Both {ias} and {oas} cannot be true.")
        if ias:
            #TODO
            pass
        raise NotImplementedError(f"Not implemented yet.")
