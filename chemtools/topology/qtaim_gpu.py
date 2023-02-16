from collections import OrderedDict
import itertools
import numpy as np
from scipy.spatial.distance import cdist

from grid.lebedev import AngularGrid

from chemtools.topology.surface import SurfaceQTAIM
from chemtools.topology.utils import (
    determine_beta_spheres,
    find_optimize_centers,
    solve_for_oas_points
)
from chemtools.topology.ode import steepest_ascent_rk45


__all__ = ["qtaim_surface_vectorize"]


def _classify_rays_as_ias_or_oas(
        maximas, all_points, all_basins, index_to_atom,
        numb_rays_to_atom, numb_rad_to_radial_shell,
):
    r"""
    Classify all rays in a molecule as either crossing the outer or inner atomic surface.

    Also provides the interval limits [r_0, r_1] of each ray that crosses the IAS.

    Parameters
    ----------
    maximas: ndarray(M, 3)
        Optimized centers of the electron density.
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
    ias, oas, ias_bnds: list[list[int]], list[list[int]], list[OrderedDict]
        A list of size `M`, that holds a list of indices of which angular pt/ray corresponds to
        either intersecting the ias or oas. The final element is a list of size `M`, of
        ordered dictionary whose keys are the indices of ias and items are the lower
        and upper-bound of the radius where the intersection occurs somewhere inbetween.

    """
    numb_maximas = len(maximas)
    oas = [[] for _ in range(numb_maximas)]  # outer atomic surface
    ias = [[] for _ in range(numb_maximas)]  # inner atomic surface.

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
    print("Index to atom ", index_to_atom, all_points.shape)
    ias_bnds = [OrderedDict() for _ in range(0, numb_maximas)]  # Keys are Points index
    np.set_printoptions(threshold=np.inf)
    for i_maxima in range(0, numb_maximas):
        print("ATom i ", i_maxima)
        print("Starting and Final index", index_to_atom[i_maxima], index_to_atom[i_maxima + 1])
        basins_a = all_basins[index_to_atom[i_maxima]:index_to_atom[i_maxima + 1]]  # Basin of atom
        points_a = all_points[index_to_atom[i_maxima]:index_to_atom[i_maxima + 1]]  # Points of atom
        print("Basins of atom ", basins_a)
        numb_rad_pts = numb_rad_to_radial_shell[i_maxima]

        i_ray = 0
        print(index_to_atom[i_maxima], numb_rad_pts)
        print("Number of angular points in this atom", numb_rays_to_atom[i_maxima])
        for i_ang in range(numb_rays_to_atom[i_maxima]):
            print("Angular pt j", i_ang)
            print("Number of radial points in this angular pt ", numb_rad_pts[i_ang])

            # Get the basin of the ray
            basins_ray = basins_a[i_ray:i_ray + numb_rad_pts[i_ang]]
            print("Basin of the ray ", basins_ray)

            # Classify basins as either OAS and IAS, if IAS, then count the number of
            #     intersections of the IAS. In addition, get the l_bnd, u_bnd of each intersection.
            # Groups basins i.e. [1, 1, 1, 0, 0, 0, 1, 1, 2, 2] -> [1, 0, 1, 2]
            group_by = [(k, list(g)) for k, g in itertools.groupby(basins_ray)]
            unique_basins = np.array([x[0] for x in group_by])
            print(basins_ray == i_maxima)
            print(unique_basins)

            # All pts in the ray got assigned to the same basin
            if len(unique_basins) == 1:
                # This implies it is an OAS point, else then it is an IAS with a bad ray.
                if unique_basins[0] == i_maxima:
                    print("OAS Point")
                    oas[i_maxima].append(i_ang)
                else:
                    # This is IAS with a bad ray, would have to re-determine the l_bnd
                    raise RuntimeError("Fix later, bad ray")
            else:
                # The point is an IAS, determine the number of intersections.
                conv_to_atom = unique_basins == i_maxima
                numb_intersections = np.sum(conv_to_atom)
                if numb_intersections >= 1:
                    print("IAS With one Intersection.")
                    # Determine lower and upper-bound Point on ray.
                    index_u_bnd = len(group_by[0][1])
                    # Determine radius from the upper and lower bound.
                    r_ubnd = np.linalg.norm(points_a[i_ray + index_u_bnd] - maximas[i_maxima])
                    r_lbnd = np.linalg.norm(
                        points_a[i_ray + index_u_bnd - 1] - maximas[i_maxima])
                    # Update containers
                    ias_bnds[i_maxima][i_ang] = [r_lbnd, r_ubnd]
                    ias[i_maxima].append(i_ang)
                    print("Radius Lower and Upper bound ", r_lbnd, r_ubnd)
                # else:
                #     r"""
                #     There are rays for example CH4, where the ray goes from basin 1 to 0 to 1
                #     again, it doesn't make much sense why this is the case, because the ray
                #     is very unlikely to do this and should have gone to the other hydrogen. But
                #     the density values is incredibly small here.
                #
                #     """
                #     print("IAS With Multiple Intersections")
                #     print(dens_func(points_atoms[i_ray:i_ray + numb_rad_pts[i_ang]]))
                # from chemtools.topology.qtaim import gradient_path_all_pts
                # basins = gradient_path_all_pts(
                #     points_atoms[i_ray:i_ray + numb_rad_pts[i_ang]], grad_func, beta_spheres, i_maxima, maximas,
                #     t_span=(0, 100), max_step=np.inf, method="LSODA",
                #     first_step=1e-7
                # )

            # import matplotlib
            # import matplotlib.pyplot as plt
            # from mpl_toolkits import mplot3d
            # matplotlib.use("Qt5Agg")
            # fig = plt.figure()
            # ax = plt.axes(projection='3d')
            # p = points_atoms[i_ray: i_ray + numb_rad_pts[i_ang]]
            # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
            # ax.scatter(maximas[:, 0], maximas[:, 1], maximas[:, 2], color="r", s=60)
            # plt.show()

            i_ray += numb_rad_pts[i_ang]
    return ias, oas, ias_bnds


def construct_all_points_of_rays_of_atoms(
    maximas, angular_pts, radial_grid, dens_func, iso_val
):
    r"""
    Construct all points of all rays across all molecules.
    """
    #  Need a way to track which points correspond to which maxima,
    #  Need a way to track which sets of points correspond to a ray
    #  index_to_atom = [0, i_1, i_1 + i_2, ..., \sum_j^M i_j] first index always zero and last
    #      always the number of points.
    numb_maximas = len(maximas)
    index_to_atom = [0] * (numb_maximas + 1)  # First index is always zero
    NUMB_RAYS_TO_ATOM = [len(ang_grid) for ang_grid in angular_pts]
    numb_rad_to_radial_shell = []  # List of List: Number of radius points per ray
    points = []
    for i in range(0, numb_maximas):
        # Construct all points on the atomic grid around atom i
        radial_shells = np.einsum("i,jk->jik", radial_grid[i], angular_pts[i])
        print("Number of radial points", len(radial_grid[i]))
        rs = maximas[i, None, None, :] + radial_shells   # has shape (K, N, 3)
        rs = np.reshape(rs, (rs.shape[0] * rs.shape[1], 3)) # has shape (KN, 3)
        print("Total number of points ", rs.shape)

        # Record information what indices it corresponds to
        numb_rad_to_radial_shell.append([len(radial_grid[i])] * NUMB_RAYS_TO_ATOM[i])

        # First remove the density values that are less than isosurface values.
        density_vals = dens_func(rs)
        indices = np.where(density_vals < iso_val)[0]
        if len(indices) != 0:
            rs = np.delete(rs, indices, axis=0)
            # Convert from index I to (i) where i is the angular index and j is the radial.
            for k in indices:
                numb_rad_to_radial_shell[i][k // len(radial_grid[i])] -= 1

        index_to_atom[i + 1] = index_to_atom[i] + rs.shape[0]  # Add what index it is
        points.append(rs)
    points = np.vstack(points)  # has shape (Product_{i=1}^M K_i N_i, 3)
    print("Total number of points ", points.shape)
    return points, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell


def _solve_intersection_of_ias(
    maximas, ias_indices, angular_pts, dens_func, grad_func, beta_spheres, bnd_err, ss_0, max_ss, tol
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
    basin_ias = [[] for _ in range(len(maximas))]  # basin ids for inner atomic surface.
    while len(ias_indices) != 0:
        # Construct New Points
        points = []
        numb_pts_per_ray = []
        for (i_maxima, i_ang, l_bnd, u_bnd, ss) in ias_indices:
            print((i_maxima, i_ang, l_bnd, u_bnd, ss))
            ray = (
                maximas[int(i_maxima)] +
                np.arange(l_bnd, u_bnd + ss, ss)[:, None] * angular_pts[int(i_maxima)][int(i_ang), :]
            )
            points.append(ray)
            numb_pts_per_ray.append(len(ray))
        points = np.vstack(points)

        # Solve for basins
        basins = steepest_ascent_rk45(
            points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0
        )
        print("Basins", basins)

        # Refine the rays further
        index_basins = 0  # Index to iterate through basins
        converge_indices = []
        for i, (i_maxima, i_ang, l_bnd, u_bnd, ss) in enumerate(ias_indices):
            basins_ray = basins[index_basins:index_basins + numb_pts_per_ray[i]]
            print((i_maxima, i_ang, l_bnd, u_bnd, ss))
            print("Basins ", i, basins_ray)

            # Basins that switch index
            i_switch = np.argmax(basins_ray != i_maxima)
            print("Index of switch ", i_switch)
            if i_switch == 0:
                raise ValueError(f"Fix this.")

            # If ss was less than bnd_err, then we converge and should stop.
            if ss <= bnd_err:
                # Take midpoint to be the radius of intersection
                radius_mid_pt = (2.0 * l_bnd + ss * i_switch) / 2.0
                r_func[int(i_maxima)][int(i_ang)] = radius_mid_pt
                basin_ias[int(i_maxima)].append(basins_ray[i_switch])
                converge_indices.append(i)  # Put in list to remove indices.
            else:
                # Update ias_indices for the next iteration for example l_bnd, u_bnd, step-size
                new_l_bnd = l_bnd + ss * (i_switch - 1)
                new_u_bnd = l_bnd + ss * (i_switch)
                new_ss = max(ss / 10.0, bnd_err)
                print("New step-size ", new_ss)
                ias_indices[i] = [i_maxima, i_ang, new_l_bnd, new_u_bnd, new_ss]

            # Update index for the next ray
            index_basins += numb_pts_per_ray[i]
            print("COnvergence indices", converge_indices)
        # Remove converged indices
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
    iso_err=1e-6,
    beta_spheres=None,
    beta_sphere_deg=21,
    ss_0=0.23,
    max_ss=0.5,
    tol=1e-7,
    optimize_centers=True
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
        Tolerance for the adaptive step-size.
    optimize_centers: bool
        If true, then the steepest-ascent is performed on the centers to find the local maximas.

    Returns
    -------
    SurfaceQTAIM:
        Object that holds all information regarding the surface of each atom.

    Notes
    -----
    The algorithm is as follows:
    1. Optimize the centers provided to obtain the local maximas.
    2. Determine the beta-spheres over all atoms.
    3. Using an angular grid and radial grid, construct all rays propogating across all atoms.
    4. Solve for each basin value for each point.
    5. Analyze the basin values and classifty each ray as either an outer-atomic or inner-atomic
        surface point.
    6. For the inner-atomic rays, find the point of intersection to the surface boundary.

    """
    if len(angular) != len(centers):
        raise ValueError(f"Length of angular {len(angular)} should be the same as the"
                         f"number of centers {len(centers)}.")
    if beta_spheres is not None and len(centers) != len(beta_spheres):
        raise ValueError(
            f"Beta sphere length {len(beta_spheres)} should match the"
            f" number of centers {len(centers)}"
        )

    # Using centers, update to the maximas
    maximas = centers
    if optimize_centers:
        # Using ODE solver to refine the maximas further.
        maximas = find_optimize_centers(maximas, grad_func)

    # Construct a radial grid for each atom by taking distance to the closest five atoms.
    #  Added an extra padding in the case of carbon in CH4
    #  TODO: the upper-bound should depend on distance to isosurface value and distance
    #         between atoms
    dist_maxs = cdist(maximas, maximas)
    distance_maximas = np.sort(dist_maxs, axis=1)[:, min(5, maximas.shape[0] - 1)]
    ss0 = 0.2
    radial_grid = [
        np.arange(0.2, x + 5.0, ss0) for x in distance_maximas
    ]
    # input("Hello")

    numb_maximas = len(maximas)
    angular_pts = []
    for ang in angular:
        if isinstance(ang, int):
            angular_pts.append(AngularGrid(degree=ang).points)
        else:
            angular_pts.append(ang.points)

    # Determine beta-spheres from a smaller angular grid
    #  Degree can't be too small or else the beta-radius is too large and IAS point got classified
    #  as OAS point. TODO: Do the angle/spherical trick then do the beta-sphere
    ang_grid = AngularGrid(degree=beta_sphere_deg)
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
    # Reduce the number of radial points that are greater than the beta-sphere.
    for i in range(0, numb_maximas):
        radial_grid[i] = radial_grid[i][radial_grid[i] >= beta_spheres[i]]

    # First step is to construct a grid that encloses all radial shells across all atoms
    points, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell = \
        construct_all_points_of_rays_of_atoms(
            maximas, angular_pts, radial_grid, dens_func, iso_val
    )
    print("Index to atom ", index_to_atom)

    # Then assign basins values for all the points.
    import time
    start = time.time()
    basins = steepest_ascent_rk45(
        points, dens_func, grad_func, beta_spheres, maximas, tol=tol, max_ss=max_ss, ss_0=ss_0
    )
    final = time.time()
    print("Basins", basins)
    print("Length of basins ", len(basins))
    print("Difference ", final - start)

    # Using basin values, classify each ray as either IAS or OAS and if IAS, find the interval
    #   along the ray that intersects the IAS.
    ias, oas, ias_bnds = _classify_rays_as_ias_or_oas(
        maximas, points, basins, index_to_atom, NUMB_RAYS_TO_ATOM, numb_rad_to_radial_shell
    )

    # The IAS is just refining the ray, till you find the exact intersection with the surface.
    # Checks docs of `_solve_intersection_of_ias` for what ias_indices is.
    ias_indices = np.array(list(
        itertools.chain.from_iterable(
            [[(i, y, ias_bnds[i][y][0], ias_bnds[i][y][1], max(ss0 / 10.0, bnd_err, ))
              for y in ias[i]] for i in range(numb_maximas)]
        )
    ))
    print("ias indices", ias_indices)
    r_func, basin_ias = _solve_intersection_of_ias(
        maximas, ias_indices, angular_pts, dens_func, grad_func, beta_spheres, bnd_err,
        tol=tol, max_ss=max_ss, ss_0=ss_0
    )

    # Solve OAS Points and updates r_func
    solve_for_oas_points(maximas, oas, radial_grid, angular_pts, dens_func, iso_val, iso_err, r_func)

    return SurfaceQTAIM(r_func, angular, maximas, oas, ias, basin_ias)
