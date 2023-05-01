import pytest
import pathlib
import numpy as np
from scipy.integrate import solve_ivp

from chemtools.wrappers import Molecule
from chemtools.topology.qtaim_gpu import qtaim_surface_vectorize

from grid.onedgrid import GaussChebyshev
from grid.rtransform import BeckeRTransform


def _run_qtaim_algorithm(fchk, degs):
    file_path = pathlib.Path(__file__).parent.resolve().__str__()[:-13]
    file_path += "data/examples/" + fchk

    mol = Molecule.from_file(file_path)
    centers = mol.coordinates
    gaussian_func = lambda pts: mol.compute_density(pts)
    gradient_func = lambda pts: mol.compute_gradient(pts)

    result = qtaim_surface_vectorize(
        degs, centers, gaussian_func, gradient_func,
        iso_val=1e-8, bnd_err=1e-5, iso_err=1e-6, optimize_centers=True
    )
    return mol, result


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("atom_kr.fchk", [20]),
        ("h2o.fchk", [30, 10, 10]),
        ("nh3.fchk", [30, 15, 15, 15]),
        ("ch4.fchk", [30, 15, 15, 15, 15])
     ]
)
def test_atomic_density_sum_to_numb_electrons(fchk, degs):
    r"""The sum of the atomic charges should equal to the charge of the molecule."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)

    numb = 500
    oned = GaussChebyshev(numb)
    rgrid = BeckeRTransform(1e-8, 4).transform_1d_grid(oned)
    density_integral = 0.0
    for i in range(len(mol.coordinates)):
        atomgrid_basin_0 = qtaim.get_atom_grid_over_basin(i, rgrid)
        dens = mol.compute_density(atomgrid_basin_0.points)
        print("Density Integral ", atomgrid_basin_0.integrate(dens))
        density_integral += atomgrid_basin_0.integrate(dens)
        print()
    print("Total Density Integral ", density_integral)
    assert np.abs(density_integral - np.sum(mol.numbers)) < 1e-2


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("atom_kr.fchk", [20]),
        ("h2o.fchk", [30, 10, 10]),
        ("nh3.fchk", [30, 15, 15, 15]),
        ("ch4.fchk", [30, 15, 15, 15, 15])
     ]
)
def test_laplacian_is_small(fchk, degs):
    r"""Laplacian over each basin should be close to zero."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)

    numb = 500
    oned = GaussChebyshev(numb)
    rgrid = BeckeRTransform(1e-8, 4).transform_1d_grid(oned)
    for i in range(len(mol.coordinates)):
        atomgrid_basin_0 = qtaim.get_atom_grid_over_basin(i, rgrid)
        laplacian = 0.25 * mol.compute_laplacian(atomgrid_basin_0.points)
        integral = atomgrid_basin_0.integrate(laplacian)

        print("Laplacian Integral ", integral)
        assert np.abs(integral) < 1e-3, "Laplacian Integral should be close to zero."


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("h2o.fchk", [15, 8, 8]),
        ("nh3.fchk", [15, 8, 8, 8]),
        ("atom_kr.fchk", [10]),
        ("ch4.fchk", [15, 8, 8, 8, 8])
     ]
)
def test_oas_isosurface_value(fchk, degs):
    r"""Test the isosurface value of the OAS points are correct."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)
    iso_val = 1e-8
    for i in range(len(mol.coordinates)):
        oas_pts = qtaim.get_oas_pts_of_basin(i)
        density = mol.compute_density(oas_pts)
        assert np.all(np.abs(density - iso_val) < 1e-6)

        # test ias pts density value is greater than isosurface value
        ias_pts = qtaim.get_ias_pts_of_basin(i)
        if len(ias_pts) != 0:  # atom_kr would not  have any ias pts.
            density = mol.compute_density(ias_pts)
            assert np.all(density > 1e-8)


@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("h2o.fchk", [15, 8, 8]),
        # ("nh3.fchk", [15, 8, 8, 8]),
        # ("ch4.fchk", [15, 8, 8, 8, 8])
     ]
)
def test_ias_basin_values(fchk, degs):
    r"""Test IAS basin value assignment is correctly assigned."""
    mol, qtaim = _run_qtaim_algorithm(fchk, degs)

    def norm_grad_func(x):
        grad = mol.compute_gradient(x)
        return grad / np.linalg.norm(grad, axis=1)[:, None]

    coords = mol.coordinates
    print("Coordinates ", coords)
    for i in range(len(coords)):
        # test ias pts density value is greater than isosurface value
        ias_indices = qtaim.ias[i]
        print(ias_indices)
        ias_ang = qtaim.generate_angular_pts_of_basin(i)[ias_indices, :]

        # None of the basin values should be the current maxima
        basin_vals_ias = np.array(qtaim.basins_ias[i])
        assert np.all(basin_vals_ias != i)

        for j in range(len(ias_ang)):
            basin_pt = basin_vals_ias[j]
            ias_pt = coords[i] + ias_ang[j] * qtaim.r_func[i][ias_indices[j]]
            # Should converge to the other basin
            ias_pt_basin = coords[i] + ias_ang[j] * (qtaim.r_func[i][ias_indices[j]] + 1e-4)
            # Should converge to the current maxima
            ias_pt_inner = coords[i] + ias_ang[j] * (qtaim.r_func[i][ias_indices[j]] - 1e-4)
            print(ias_pt, basin_pt)

            # Should ODE on both ias_pt_basin and ias_pt_inner and make sure it converges
            #  to the current maxima
            sol = solve_ivp(
                lambda t, x: norm_grad_func(np.array([x]))[0].T,
                y0=ias_pt_basin,
                t_span=(0, 100),
                method="DOP853",
                max_step=np.inf
            )["y"][:, -1]
            print(sol)
            assert np.all(np.abs(sol - coords[basin_pt]) < 1e-1)

            sol = solve_ivp(
                lambda t, x: norm_grad_func(np.array([x]))[0].T,
                y0=ias_pt_inner,
                t_span=(0, 100),
                method="DOP853",
                max_step=np.inf
            )["y"][:, -1]
            print(sol)
            assert np.all(np.abs(sol - coords[i]) < 1e-1)

            print("")



@pytest.mark.parametrize(
    "fchk, degs",
     [
        ("h2o.fchk", [15, 8, 8]),
        ("nh3.fchk", [15, 8, 8, 8]),
        ("atom_kr.fchk", [10]),
        ("ch4.fchk", [15, 8, 8, 8, 8])
     ]
)
def test_outer_atomic_surface_is_correctly_assigned(fchk, degs):
    pass
