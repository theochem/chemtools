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
import pytest
import pathlib
import numpy as np
from scipy.integrate import solve_ivp

from chemtools.wrappers import Molecule
from chemtools.topology.qtaim import qtaim_surface_vectorize

from grid.onedgrid import UniformInteger
from grid.rtransform import PowerRTransform


def _run_qtaim_algorithm(fchk, degs, iso_val=1e-10, iso_err=1e-5, bnd_err=1e-5, ss_0=0.01, max_ss=0.1, tol=1e-7):
    file_path = pathlib.Path(__file__).parent.resolve().__str__()[:-13]
    file_path += "data/examples/" + fchk

    mol = Molecule.from_file(file_path)
    centers = mol.coordinates
    gaussian_func = lambda pts: mol.compute_density(pts)
    gradient_func = lambda pts: mol.compute_gradient(pts)

    result = qtaim_surface_vectorize(
        degs, centers, gaussian_func, gradient_func,
        iso_val=iso_val, bnd_err=bnd_err, iso_err=iso_err, optimize_centers=True,
        ss_0=ss_0, max_ss=max_ss, tol=tol
    )
    return mol, result


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("h2o.fchk", [30, 10, 10]),
        ("nh3.fchk", [30, 10, 10, 10]),
        ("ch4.fchk", [25, 15, 15, 15, 15])
     ]
)
def test_atomic_density_sum_to_numb_electrons(fchk, degs):
    r"""The sum of the atomic charges should equal to the charge of the molecule."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)

    numb = 350
    oned = UniformInteger(numb)
    density_integral = 0.0
    for i in range(len(mol.coordinates)):
        b = max(np.max(qtaim.r_func[i][qtaim.ias[i]]), np.max(qtaim.r_func[i][qtaim.oas[i]])) + 1.0
        rgrid = PowerRTransform(1e-10, b).transform_1d_grid(oned)
        atomgrid_basin_0 = qtaim.get_atom_grid_over_basin(i, rgrid)

        dens = mol.compute_density(atomgrid_basin_0.points)
        print("Density Integral ", atomgrid_basin_0.integrate(dens))
        density_integral += atomgrid_basin_0.integrate(dens)

    print("Total Density Integral ", density_integral)
    assert np.abs(density_integral - np.sum(mol.numbers)) < 1e-2


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("h2o.fchk", [30, 10, 10]),
        ("nh3.fchk", [30, 10, 10, 10]),
        ("ch4.fchk", [25, 15, 15, 15, 15])
     ]
)
def test_laplacian_is_small(fchk, degs):
    r"""Laplacian over each basin should be close to zero."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)

    numb = 500
    oned = UniformInteger(numb)
    for i in range(len(mol.coordinates)):
        b = max(np.max(qtaim.r_func[i][qtaim.ias[i]]), np.max(qtaim.r_func[i][qtaim.oas[i]])) + 1.0
        rgrid = PowerRTransform(1e-10, b).transform_1d_grid(oned)
        atomgrid_basin_0 = qtaim.get_atom_grid_over_basin(i, rgrid)
        laplacian = 0.25 * mol.compute_laplacian(atomgrid_basin_0.points)
        integral = atomgrid_basin_0.integrate(laplacian)

        print("Laplacian Integral ", integral)
        assert np.abs(integral) < 1e-3, "Laplacian Integral should be close to zero."


@pytest.mark.parametrize(
    "fchk, degs, iso_val, iso_err",
     [
        ("h2o.fchk", [25, 15, 15], 0.001, 1e-5),
        ("h2o.fchk", [10, 25, 20], 1e-10, 1e-6),
        ("h2o.fchk", [50, 10, 15], 1e-12, 1e-7),
        ("h2o.fchk", [25, 15, 15], 0.01, 1e-7),
        ("nh3.fchk", [25, 15, 15, 15], 0.001, 1e-6),
        ("ch4.fchk", [25, 15, 15, 15, 15], 1e-10, 1e-5)
     ]
)
def test_oas_isosurface_value(fchk, degs, iso_val, iso_err):
    r"""Test the isosurface value of the OAS points are correct."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs, iso_val, iso_err)
    for i in range(len(mol.coordinates)):
        oas_pts = qtaim.get_oas_pts_of_basin(i)
        density = mol.compute_density(oas_pts)
        print(np.abs(density - iso_val))
        assert np.all(np.abs(density - iso_val) < iso_err)

        # test ias pts density value is greater than isosurface value
        ias_pts = qtaim.get_ias_pts_of_basin(i)
        if len(ias_pts) != 0:  # atom_kr would not  have any ias pts.
            density = mol.compute_density(ias_pts)
            assert np.all(np.abs(density - iso_val) > iso_err)


@pytest.mark.parametrize(
    "fchk, degs, bnd_err",
     [
        ("h2o.fchk", [15, 8, 8], 1e-5),
        ("h2o.fchk", [15, 8, 8], 1e-4),
        ("h2o.fchk", [15, 8, 8], 1e-3),
        ("nh3.fchk", [15, 8, 8, 8], 1e-5),
        ("ch4.fchk", [15, 8, 8, 8, 8], 1e-5)
     ]
)
def test_ias_basin_values(fchk, degs, bnd_err):
    r"""Test IAS basin value assignment is correctly assigned."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs, bnd_err=bnd_err, ss_0=0.01, max_ss=0.1, tol=1e-10)

    def norm_grad_func(x):
        grad = mol.compute_gradient(x)
        return grad / np.linalg.norm(grad, axis=1)[:, None]

    coords = qtaim.maximas
    print("Coordinates ", coords)
    for i in range(0, len(coords)):
        # test ias pts density value is greater than isosurface value
        ias_indices = qtaim.ias[i]
        print("atom i ", i)
        print(ias_indices)
        ias_ang = qtaim.generate_angular_pts_of_basin(i)[ias_indices, :]

        # None of the basin values should be the current maxima
        basin_vals_ias = np.array(qtaim.basins_ias[i])
        assert np.all(basin_vals_ias != i)

        for j in range(len(ias_ang)):
            basin_pt = basin_vals_ias[j]
            ias_pt = coords[i] + ias_ang[j] * qtaim.r_func[i][ias_indices[j]]
            # Should converge to the other basin
            ias_pt_basin = coords[i] + ias_ang[j] * (qtaim.r_func[i][ias_indices[j]] + bnd_err * 10)
            # Should converge to the current maxima
            ias_pt_inner = coords[i] + ias_ang[j] * (qtaim.r_func[i][ias_indices[j]] - bnd_err * 10)
            print(ias_pt, basin_pt)

            # Should ODE on both ias_pt_basin and ias_pt_inner and make sure it converges
            #  to the current maxima
            sol = solve_ivp(
                lambda t, x: norm_grad_func(np.array([x]))[0].T,
                y0=ias_pt_basin,
                t_span=(0, 8),
                method="RK45",
                first_step=bnd_err,
                max_step=0.23,
                atol=1e-7,
                rtol=1e-4,
            )["y"][:, -1]
            print(sol)
            assert np.all(np.abs(sol - coords[basin_pt]) < 1e-1)

            sol = solve_ivp(
                lambda t, x: norm_grad_func(np.array([x]))[0].T,
                y0=ias_pt_inner,
                t_span=(0, 8),
                method="RK45",
                first_step=bnd_err,
                max_step=0.23,
                atol=1e-7,
                rtol=1e-4,
            )["y"][:, -1]
            print(sol)
            assert np.all(np.abs(sol - coords[i]) < 1e-1)

            print("")
