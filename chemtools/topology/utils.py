import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import root_scalar

from chemtools.topology.ode import steepest_ascent_rk45, gradient_path

r"""
Utility functions that is common between the QTAIM algorithms.
"""

__all__ = [
    "solve_for_oas_points",
    "construct_radial_grids",
    "find_optimize_centers",
    "determine_beta_spheres",
    "solve_for_isosurface_pt"
]


def solve_for_oas_points(
        maximas, oas, radial_grid, angular_pts, dens_func, iso_val, iso_err, r_func
):
    r"""
    For each index in outer-atomic surface (OAS) solves for the isovalue point along a ray.

    This is stored inside `r_func`.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        The maximas of the density
    oas: list[list]
        List of indices that correspond to angular points whose ray intersect the isosurface
        of the electron density.
    radial_grid: list[ndarray]
        List of radial grids (arrays on zero to infinity) correspond to each maxima.
    angular_pts: ndarray(N,)
        The angular points on the sphere.
    dens_func: callable()
        The electron density function.
    iso_val: float
        The isosurface value that is to be solved
    iso_err: float
        The isosurface error
    r_func: list[ndarray()]
        This holds the radial coordinate on the ray that intersects the OAS.

    """
    for i_maxima in range(len(maximas)):
        maxima = maximas[i_maxima]
        radial = radial_grid[i_maxima]
        ang_pts = angular_pts[i_maxima]
        for i_oas in oas[i_maxima]:
            # Construct upper and lower bound of the isosurface equation
            ang_pt = ang_pts[i_oas]
            iso_eq = np.abs(dens_func(maxima + ang_pt * radial[:, None]) - iso_val)
            i_iso = np.argsort(iso_eq)[0]
            l_bnd = radial[i_iso - 1] if i_iso >= 0 else radial[i_iso] / 2.0
            u_bnd = radial[i_iso + 1] if i_iso + 1 < len(radial) else radial[i_iso] * 2.0
            # Solve for the isosurface point
            oas_pt = solve_for_isosurface_pt(
                l_bnd, u_bnd, maxima, ang_pt, dens_func, iso_val, iso_err
            )
            # print("Check isosurface pt", oas_pt, dens_func(np.array([oas_pt])))
            # Record them
            r_func[i_maxima][i_oas] = np.linalg.norm(oas_pt - maxima)


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


def construct_radial_grids(maximas, min=0.2, pad=5.0, ss0=0.23):
    r"""
    Construct radial grids

    Parameters
    ----------
    maximas: ndarray(M, 3)
        Coordinates of the maximas
    min: float
        The minimum radial value on [0, \infty) to construct a radial grid.
    pad: float
        Extra padding to add to make sure the radial grid covers the intersection
        with the inter-atomic and outer-atomic surfaces.
    ss0: float
        The step-size of the uniform radial grid.

    Returns
    -------
    list[ndarray]:
        List of radial grid of length number of maximas. The radial grid are uniform
        grids that start at `min` and end on the maximum distance to the fifth atom plus
        an extra padding with stepsize `ss0`.

    """
    # Construct a radial grid for each atom by taking distance to the closest five atoms.
    #  Added an extra padding in the case of carbon in CH4
    #  TODO: the upper-bound should depend on distance to isosurface value and distance
    #         between atoms
    dist_maxs = cdist(maximas, maximas)
    distance_maximas = np.sort(dist_maxs, axis=1)[:, min(5, maximas.shape[0] - 1)]
    print(cdist(maximas, maximas))
    radial_grid = [
        np.arange(min, x + pad, ss0) for x in distance_maximas
    ]
    return radial_grid


def determine_beta_spheres(beta_spheres, maximas, radial_grid, angular_pts, dens_func, grad_func):
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
    for i_maxima, maxima in enumerate(maximas):
        if beta_spheres[i_maxima] == initial_beta_sph:
            optimal_rad = -np.inf
            for rad_pt in radial_grid[i_maxima]:
                if rad_pt > initial_beta_sph:
                    # Determine the points on the sphere with this radius
                    pts = maxima + rad_pt * angular_pts
                    print(pts)
                    basins = steepest_ascent_rk45(
                        pts, dens_func, grad_func, beta_spheres, maximas
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
