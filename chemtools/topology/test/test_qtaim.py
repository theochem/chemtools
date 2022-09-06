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
from scipy.spatial import ConvexHull
from scipy.stats import special_ortho_group
from scipy.spatial.transform.rotation import Rotation

from chemtools.topology.qtaim import qtaim, _get_area_of_coplanar_polygon, qtaim_surface
from grid.cubic import Tensor1DGrids, UniformGrid
from grid.onedgrid import OneDGrid, GaussLaguerre
from grid.becke import BeckeWeights
from grid.molgrid import MolGrid

import pytest


def _get_cubic_grid(l_bnd, u_bnd, ss):
    oned = np.arange(l_bnd, u_bnd, ss)
    oned_grid = OneDGrid(oned, oned)
    return Tensor1DGrids(oned_grid, oned_grid, oned_grid)


def _get_molecular_grid(centers):
    # Construct two atomic grid whose Lebedev degrees increase.
    oned = np.arange(0.001, 2.0, 0.25)
    rgrid = OneDGrid(oned, oned)
    numbs = np.array([1] * centers.shape[0])
    return MolGrid.from_preset(atnums=numbs, atcoords=centers, rgrid=rgrid,
                               preset="coarse", aim_weights=BeckeWeights())


@pytest.mark.parametrize("shape", np.random.randint(5, 30, size=(10, 3)))
def test_with_simple_zero_flux_surface_of_two_exponentials_on_cubic_grid(shape):
    r"""Test two exponentials symmetrically spaced apart on a cubic grid.

    One is multipled by 0.98 to remove symmetry. Zero-flux surface occurs at (0, y, z).
    """
    centers = np.array([[-1, 0, 0], [1, 0, 0]])
    gaussian_func = lambda pts: np.exp(-np.linalg.norm(pts - centers[0], axis=1)) + \
                                0.98 * np.exp(-np.linalg.norm(pts - centers[1], axis=1))

    # Define Grid and evaluate the density on the grid
    origin = np.array([-1.5, -1.5, -1.5])
    shape = np.array(shape)
    axes = np.eye(3) * (1.5 + 1.5) / (shape - 1).T
    grid = UniformGrid(origin, axes, shape=shape)
    gaussians = gaussian_func(grid.points)
    result = qtaim(grid, gaussians, bounding_box=False)

    assert result["basin_cont"].shape[1] == 2
    for i, pt in enumerate(grid.points):
        basin_assigned = result["basin_cont"][i].argmax() + 1
        if basin_assigned == 1:
            # Assert it is part of basin one.
            assert pt[0] <= 0.0
        elif basin_assigned == 2:
            # assert it is part of basin two.
            assert pt[0] >= 0.0


@pytest.mark.parametrize("shape", np.random.randint(5, 30, size=(3, 3)))
def test_basin_are_correctly_assigned_against_accurate_scipy_ode_solver(shape):
    centers = np.array([[-1.5, 0, 0], [1.5, 0, 0]])
    gaussian_func = lambda pts: np.exp(-np.linalg.norm(pts - centers[0], axis=1)**2.0) + \
                                np.exp(-np.linalg.norm(pts - centers[1], axis=1)**2.0)

    gradient_func = lambda pts: (
         -2.0 * ((pts - centers[0]) *  np.exp(-np.linalg.norm(pts - centers[0], axis=1) ** 2.0).T
                + (pts - centers[1]) * np.exp(-np.linalg.norm(pts - centers[1], axis=1)**2.0).T)
    )

    # Define Grid and evaluate the density on the grid
    origin = np.array([-1.5, -1.5, -1.5])
    shape = np.array(shape)
    axes = np.eye(3) * (1.5 + 1.5) / (shape - 1).T
    grid = UniformGrid(origin, axes, shape=shape)
    gaussians = gaussian_func(grid.points)

    # Do qtaim on the grid.
    result = qtaim(grid, gaussians, bounding_box=False, grad_func=gradient_func)
    maximas = grid.points[result["maxima_indices"]]

    # take random sample of points
    numb_samples = 1000
    sample_indices = np.random.randint(0, grid.points.shape[0], size=numb_samples)
    for i_samp in sample_indices:
        pt_samp = grid.points[i_samp]

        basin_weights = result["basin_cont"][i_samp].toarray()

        # Only check with points that are not on the zero-flux surface.
        if np.all(np.abs(basin_weights - 0.5) > 0.01) and np.abs(pt_samp[0]) > axes[0, 0]:
            sol = solve_ivp(
                lambda t, x: gradient_func(np.array([x]))[0].T,
                y0=pt_samp,
                t_span=(0, 1000),
                method="DOP853",
                max_step=50
            )
            # print("solution ", sol,  " maximas ", )
            print("Pt Sample", pt_samp, "Basin of it ", result["basin_cont"][i_samp], basin_weights)

            # basin assigned by the algorithnm
            basin_assigned = result["basin_cont"][i_samp].toarray().argmax() + 1
            # basin assigned by the ode
            basin_assigned_ode = np.linalg.norm(sol["y"][:, -1] - maximas, axis=1).argmin() + 1
            print(basin_assigned, basin_assigned_ode)
            assert basin_assigned == basin_assigned_ode


@pytest.mark.parametrize("num_pts", np.random.randint(4, 8, size=(200,)))
def test_get_area_of_coplanar_points_against_scipy_convexhull(num_pts):
    r"""Test finding the area of coplanar points against SciPy convex hull algorithm."""
    # seems that this is only accurate up to seven points, couldn't get it working past 7
    # unless the convex, coplanar polygon was a "nice" polygon.
    origin, pt = np.random.random((2,3))
    vertices = np.zeros((num_pts, 3))
    vertices[0] = pt

    # Rotate the points from finding a rotation matrix that rotates based on total_deg
    total_deg = 360 / (num_pts + np.random.randint(2, 9))
    # rotate x,y,z by total_deg
    rot_mat = Rotation.from_euler('xyz', [total_deg, total_deg, total_deg], degrees=True)
    rot_mat = rot_mat.as_matrix()
    for i in range(1, num_pts):
        vertices[i] = origin + rot_mat.dot(vertices[i - 1] - origin)
    desired = _get_area_of_coplanar_polygon(vertices)
    convex = ConvexHull(vertices, qhull_options="QJ")
    print(desired, convex.area, convex.area / 2.0)
    assert np.abs(desired - convex.area / 2.0) < 1e-8


def test_get_area_of_coplanar_points_against_perfect_shapes():
    r"""Test get area of copolanar against squares and rectangles."""
    square = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]]) * 2
    desired = _get_area_of_coplanar_polygon(square)
    print(desired,  ConvexHull(square, qhull_options="QJ").area / 2.0)
    assert np.abs(desired - ConvexHull(square, qhull_options="QJ").area / 2.0) < 1e-8
    assert np.abs(desired - 2 * 2) < 1e-8

    # Rotate square
    rot_matrix = special_ortho_group.rvs(3)
    square = square.dot(rot_matrix)
    desired = _get_area_of_coplanar_polygon(square)
    print(desired,  ConvexHull(square, qhull_options="QJ").area / 2.0)
    assert np.abs(desired - ConvexHull(square, qhull_options="QJ").area / 2.0) < 1e-8
    assert np.abs(desired - 2 * 2) < 1e-8

    # Test on rectangle
    rectangle = np.array([[0, 0, 0], [1, 0, 0], [1, 5, 0], [0, 5, 0]])
    desired = _get_area_of_coplanar_polygon(rectangle)
    print(desired,  ConvexHull(rectangle, qhull_options="QJ").area / 2.0)
    assert np.abs(desired - ConvexHull(rectangle, qhull_options="QJ").area / 2.0) < 1e-8
    assert np.abs(desired - 5) < 1e-8

    # Rotate rectangle
    rectangle = rectangle.dot(rot_matrix)
    desired = _get_area_of_coplanar_polygon(rectangle)
    print(desired,  ConvexHull(rectangle, qhull_options="QJ").area / 2.0)
    assert np.abs(desired - ConvexHull(rectangle, qhull_options="QJ").area / 2.0) < 1e-8
    assert np.abs(desired - 5) < 1e-8


@pytest.mark.parametrize("shape", np.random.randint(12, 30, size=(10, 3)))
def test_qtaim_cubic_vs_qtaim_voronoi_algorithms(shape):
    r"""Test QTAIM algorithm using a cubic grid and voronoi style."""
    centers = np.array([[-1, 0, 0], [1, 0, 0]])
    # multiply by 0.98 to order the points uniquely when you sort in the qtaim algorithm.
    gaussian_func = lambda pts: np.exp(-np.linalg.norm(pts - centers[0], axis=1)) + \
                                0.99 * np.exp(-np.linalg.norm(pts - centers[1], axis=1))

    # Define Grid and evaluate the density on the grid
    origin = np.array([-1.5, -1.5, -1.5])
    shape = np.array(shape)
    axes = np.eye(3) * (1.5 + 1.5) / (shape - 1).T
    grid = UniformGrid(origin, axes, shape=shape)
    gaussians = gaussian_func(grid.points)

    result_voronoi = qtaim(grid.points, gaussians, num_centers=2)
    result_cubic = qtaim(grid, gaussians, num_centers=2)

    assert result_voronoi["basin_cont"].shape[1] == 2
    assert result_cubic["basin_cont"].shape[1] == 2

    indices = np.argsort(gaussians)[::-1]
    gaussians = gaussians[indices]
    for i, pt in enumerate(grid.points[indices, :]):
        # Points on the boundary of the cube may have different areas
        if np.all(np.abs(np.abs(pt) - 1.5) > 1e-4):
            basin_weights_cubic = result_cubic["basin_cont"][indices[i]].toarray()
            basin_weights_voronoi = result_voronoi["basin_cont"][indices[i]].toarray()
            print(i, indices[i], pt, basin_weights_cubic, basin_weights_voronoi, gaussians[i])
            assert np.all(np.abs(basin_weights_cubic - basin_weights_voronoi) < 1e-8)
        else:
            # Atleast check if their assigned basins are the same.
            basin_weights_cubic = result_cubic["basin_cont"][indices[i]].toarray().argmax()
            basin_weights_voronoi = result_voronoi["basin_cont"][indices[i]].toarray().argmax()
            assert basin_weights_cubic == basin_weights_voronoi


@pytest.mark.parametrize("use_gradient", [True, False])
def test_integral_of_gaussians_using_cubic_grid_up_to_three_decimal(use_gradient):
    r"""Test the integral of the basins of two Gaussians that are far apart."""
    centers = np.array([[-1, 0, 0], [1, 0, 0]])
    # multiply by 0.98 to order the points uniquely when you sort in the qtaim algorithm.
    alpha = 30
    gaussian_func = lambda pts: np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1)**2.0) + \
        0.99 * np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1)**2.0)


    gradient_func = lambda pts: (
         -2.0 * ((pts - centers[0]) *  np.exp(-np.linalg.norm(pts - centers[0], axis=1) ** 2.0).T
                + (pts - centers[1]) * np.exp(-np.linalg.norm(pts - centers[1], axis=1)**2.0).T)
    )

    # Define Grid and evaluate the density on the grid
    origin = np.array([-1.5, -1.5, -1.5])
    shape = np.array([50, 45, 40])
    axes = np.eye(3) * (1.5 + 1.5) / (shape - 1).T
    print(axes)
    grid = UniformGrid(origin, axes, shape=shape, weight="Rectangle")
    gaussians = gaussian_func(grid.points)

    if use_gradient:
        result_voronoi = qtaim(grid, gaussians, num_centers=2, grad_func=gradient_func)
    else:
        result_voronoi = qtaim(grid, gaussians, num_centers=2, grad_func=None)

    for i in range(2):
        integral = grid.integrate(result_voronoi["basin_cont"][:, i].toarray().ravel() * gaussians)
        print(integral)
        factor = 1.0 if i == 0 else 0.99
        print(np.sqrt(np.pi / alpha)**3.0 * factor)
        assert np.abs(integral - factor * np.sqrt(np.pi / alpha)**3.0) < 1e-6


def test_qtaim_line_search():
    r"""TODO."""
    centers = np.array([[-1, 0, 0], [1, 0, 0]])
    # multiply by 0.98 to order the points uniquely when you sort in the qtaim algorithm.
    alpha = 3
    gaussian_func = lambda pts: np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1)**2.0) + \
        np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1)**2.0)


    gradient_func = lambda pts: (
         -2.0 * alpha * (
            (pts - centers[0]) *  np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1) ** 2.0).T
            + (pts - centers[1]) * np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1)**2.0).T
        )
    )


    oned = np.arange(1e-4, 2, 0.1)
    rgrid = OneDGrid(oned, np.ones(len(oned)) * 0.5)
    iso_val = 1e-5
    result = qtaim_surface(rgrid, 10, centers, gaussian_func, gradient_func,
                           iso_val=iso_val,
                           bnd_err=1e-5,
                           iso_err=1e-6, dens_cutoff=1e-9, beta_sphere=[0.8, 0.8])
    import matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits import mplot3d
    matplotlib.use("Qt5Agg")
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    q = result.generate_pts_on_surface(0)
    p = result.get_ias_pts_of_basin(0)
    ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
    p = result.get_oas_pts_of_basin(0)
    ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
    plt.show()



class TestQTAIMSurfaceOnTwoBodyGaussian():
    def gaussian_func(self, pts, centers, alpha):
        return np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1)**2.0) + \
                    np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1)**2.0)

    def gradient_gaussian(self, pts, centers, alpha):
        return (
            -2.0 * alpha * (
            (pts - centers[0]) *  np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1) ** 2.0).T
            + (pts - centers[1]) * np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1)**2.0).T
            )
        )

    @pytest.mark.parametrize(
        "centers", [
            np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            np.array([[-0.5, 0., 0.], [0.0, 0.5, 0.]]),
            np.vstack((np.random.uniform(-1, 0, size=(3,)), np.random.uniform(0, 1, size=(3,))))
        ],
    )
    @pytest.mark.parametrize("iso_val", [1e-5, 1e-4, 1e-2])
    def test_outer_atomic_surface_has_correct_isosurface_values(self, centers, iso_val):
        r"""Test outer atomic surface has correct isosurface value."""
        rgrid = GaussLaguerre(15)
        alpha = 2
        gaussian_func = lambda pts: self.gaussian_func(pts, centers, alpha)
        gradient_func = lambda pts: self.gradient_gaussian(pts, centers, alpha)
        iso_err = 1e-6
        result = qtaim_surface(rgrid, 15, centers, gaussian_func, gradient_func,
                               iso_val=iso_val,
                               bnd_err=1e-5, iso_err=iso_err, dens_cutoff=1e-9,
                               optimize_centers=False)
        # Test that the outer surface gives the correct
        for i in range(0, 2):
            oas_0 = result.get_oas_pts_of_basin(i)
            np.set_printoptions(threshold=np.inf)
            print(gaussian_func(oas_0))
            assert np.all(np.abs(gaussian_func(oas_0) - iso_val) < iso_err)

    @pytest.mark.parametrize("beta_sphere", [None, [0.8, 0.8]])
    @pytest.mark.parametrize("bnd_err", [1e-5, 1e-3])
    @pytest.mark.parametrize("alpha,refine", [[5, True], [1, True], [0.5, False]])
    def test_inner_atomic_surface_is_correct_on_simple_example(
            self, beta_sphere, bnd_err, alpha, refine
    ):
        r"""Test inner atomic surface lies exactly on x-axis on this example."""
        centers = np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        gaussian_func = lambda pts: self.gaussian_func(pts, centers, alpha)
        gradient_func = lambda pts: self.gradient_gaussian(pts, centers, alpha)

        rgrid = GaussLaguerre(10)
        result = qtaim_surface(rgrid, 20, centers, gaussian_func, gradient_func,
                               iso_val=1e-4,
                               bnd_err=bnd_err, iso_err=1e-6, dens_cutoff=1e-9,
                               optimize_centers=False, refine=refine)
        for i in range(0, 2):
            ias_0 = result.get_ias_pts_of_basin(i)
            assert np.all(np.abs(ias_0[:, 0]) < bnd_err)

    @pytest.mark.parametrize(
        "centers, refine", [
            [np.array([[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]), True],
            [np.vstack(
                (np.random.uniform(-1, -0.1, size=(3,)), np.random.uniform(0.1, 1, size=(3,)))
            ), False],
            [np.vstack(
                (np.random.uniform(-1, -0.1, size=(3,)), np.random.uniform(0.1, 1, size=(3,)))
            ), True],
            [np.vstack(
                (np.random.uniform(-1, -0.1, size=(3,)), np.random.uniform(0.1, 1, size=(3,)))
            ), False]
        ],
    )
    def test_outer_atomic_surface_is_correctly_assigned_to_basin(self, centers, refine):
        alpha = 0.75
        gaussian_func = lambda pts: self.gaussian_func(pts, centers, alpha)
        gradient_func = lambda pts: self.gradient_gaussian(pts, centers, alpha)

        rgrid = OneDGrid(np.arange(0., 5, 0.5), np.arange(0., 5, 0.5))
        result = qtaim_surface(rgrid, 15, centers, gaussian_func, gradient_func,
                               iso_val=1e-4,
                               bnd_err=1e-6, iso_err=1e-6, dens_cutoff=1e-6,
                               refine=refine)

        # Test that points on oas all converge to the maxima and no other.
        for i in range(0, 2):
            oas_0 = result.get_oas_pts_of_basin(i)

            for pt in oas_0:
                sol = solve_ivp(
                    lambda t, x: gradient_func(np.array([x]))[0].T,
                    y0=pt,
                    t_span=(0, 10000),
                    method="Radau",  # DOP853
                    max_step=10,
                    atol=1e-9,
                    rtol=1e-5
                )

                print(pt, sol["y"][:, -1], centers)
                if not np.all(np.abs(sol["y"][:, -1] - result.maximas[i]) < 1e-2):
                    import matplotlib
                    import matplotlib.pyplot as plt
                    from mpl_toolkits import mplot3d
                    matplotlib.use("Qt5Agg")
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    q = result.generate_pts_on_surface(0)
                    p = result.get_ias_pts_of_basin(0)
                    ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
                    p = result.get_oas_pts_of_basin(0)
                    ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
                    ax.scatter(sol["y"][:, -1][0], sol["y"][:, -1][1], sol["y"][:, -1][2], color="g", s=60)
                    ax.scatter(pt[0], pt[1], pt[2], color="y", s=60)
                    plt.show()
                assert np.all(np.abs(sol["y"][:, -1] - result.maximas[i]) < 1e-2)

    def test_integration_of_basin(self):
        centers = np.array([[-1, 0, 0], [1, 0, 0]])
        alpha = 3
        gaussian_func = lambda pts: np.exp(
            -alpha * np.linalg.norm(pts - centers[0], axis=1) ** 2.0) + \
                                    np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1) ** 2.0)

        gradient_func = lambda pts: (
                -2.0 * alpha * (
                (pts - centers[0]) * np.exp(-alpha * np.linalg.norm(pts - centers[0], axis=1) ** 2.0).T
                + (pts - centers[1]) * np.exp(-alpha * np.linalg.norm(pts - centers[1], axis=1) ** 2.0).T
            )
        )

        oned = np.arange(1e-4, 2, 0.1)
        rgrid = OneDGrid(oned, np.ones(len(oned)) * 0.5)
        result = qtaim_surface(rgrid, 10, centers, gaussian_func, gradient_func,
                               iso_val=1e-5,
                               bnd_err=1e-5, iso_err=1e-6, dens_cutoff=1e-9,
                               beta_sphere=[0.8, 0.8])

        # Test integration
        desired = np.sqrt(np.pi / alpha) ** 3.0
        for i in range(2):
            atomgrid_basin_0 = result.get_atom_grid_over_basin(i)
            true = atomgrid_basin_0.integrate(gaussian_func(atomgrid_basin_0.points))
            assert np.abs(true - desired) < 1e-8

    # def test_ch4(self):
    #     alpha = 3
    #     from chemtools.wrappers import Molecule
    #     mol = Molecule.from_file(r"/home/pally/PythonProjects/chemtools/chemtools/data/ch4_uhf_ccpvdz.fchk")
    #     centers = mol.coordinates
    #     gaussian_func = lambda pts: mol.compute_density(pts)
    #     gradient_func = lambda pts: mol.compute_gradient(pts)
    #
    #     oned = np.arange(1e-4, 5, 0.5)
    #     rgrid = OneDGrid(oned, np.ones(len(oned)) * 0.5)
    #     result = qtaim_surface(rgrid, 20, centers, gaussian_func, gradient_func,
    #                            iso_val=1e-5, ss_watershed=1e-2,
    #                            bnd_err=1e-5, iso_err=1e-6, dens_cutoff=1e-9,
    #                            optimize_centers=False, refine=True)
    #
    #     import matplotlib
    #     import matplotlib.pyplot as plt
    #     from mpl_toolkits import mplot3d
    #     matplotlib.use("Qt5Agg")
    #     for i in range(0, centers.shape[0]):
    #         fig = plt.figure()
    #         ax = plt.axes(projection='3d')
    #         p = centers
    #         ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
    #         p = result.get_ias_pts_of_basin(i)
    #         ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
    #         p = result.get_oas_pts_of_basin(i)
    #         ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
    #         plt.show()
