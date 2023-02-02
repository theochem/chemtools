from collections import OrderedDict
import itertools
import numpy as np

from scipy.optimize import root_scalar
from scipy.spatial.distance import cdist
from scipy.sparse import lil_matrix

from grid.atomgrid import AtomGrid
from grid.cubic import UniformGrid, _HyperRectangleGrid
from grid.lebedev import AngularGrid
from grid.utils import convert_cart_to_sph

from chemtools.topology.qtaim import determine_beta_spheres, _optimize_centers, solve_for_isosurface_pt, SurfaceQTAIM

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d

from chemtools.topology.ode import steepest_ascent_rk45


class EntireGrid():
    r"""
    Class that holds the entire angular grid over all points

    """
    def __init__(self):
        pass


def qtaim_surface_vectorize(angular, centers, dens_func, grad_func, iso_val=0.001,
                            dens_cutoff=1e-5, bnd_err=1e-4, iso_err=1e-6,
                            beta_spheres=None, optimize_centers=True, refine=False):
    r"""

    Notes
    -----
    1. Determine the beta-spheres.
    2.


    - Highest angular degree is 131 for Lebedev-Laikov grid with 5810 points. For
      symmetric spherical t-design degree is 325 with 52978 points.  The number of radial points
      is N with number of maximas is M, then the number of points is 5810*M*N, 52978*M*N.
      Let N=100, then we have 58100*M and 529780*M, which you'll need 80 atoms to reach
      1 GB of memory.

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
    #  Degree can't be too small or else the beta-radius is too large and IAS point got classified
    #  as OAS point. TODO: Do the angle/spherical trick then do the beta-sphere
    ang_grid = AngularGrid(degree=10)
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


    # TODO: Better to turn these into classes
    r_func = [np.zeros((numb_ang_pts,), dtype=np.float64) for _ in range(numb_maximas)]
    oas = [[] for _ in range(numb_maximas)]  # outer atomic surface
    ias = [[] for _ in range(numb_maximas)]  # inner atomic surface.
    basin_ias = [[] for _ in range(numb_maximas)]  # basin ids for inner atomic surface.
    refined_ang = [] if refine else None
    maxima_to_do = range(0, numb_maximas) if type(refine) == type(True) else [refine]  # refining

    # First step is to construct a grid that encloses all radial shells across all atoms
    # Probably best to make a datastructure class for this, because different atoms
    #  can have different radial grids, angular grids and converge differently.
    #  Need a way to track which points correspond to which maxima,
    #  Need a way to track which sets of points correspond to a ray
    #  index_to_atom = [0, i_1, i_1 + i_2, ..., \sum_j^M i_j] first index always zero and last
    #      always the number of points.
    index_to_atom = [0] * (numb_maximas + 1)  # First index is always zero
    NUMB_RAYS_TO_ATOM = [len(angular_pts)] * numb_maximas
    numb_rad_to_radial_shell = []  # List of List: Number of radius points per ray
    points = []
    for i in range(0, numb_maximas):
        # Construct all points on the atomic grid around atom i
        radial_shells = np.einsum("i,jk->jik", radial_grid[i], angular_pts)
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
            print(density_vals)
            print(indices)
            # Convert from index I to (i) where i is the angular index and j is the radial.
            for k in indices:
                numb_rad_to_radial_shell[i][k // len(radial_grid[i])] -= 1
            input("Density vals small")

        index_to_atom[i + 1] = index_to_atom[i] + rs.shape[0]  # Add what index it is
        points.append(rs)
    points = np.vstack(points)  # has shape (Product_{i=1}^M K_i N_i, 3)
    print(points.shape)
    print("Index to atom ", index_to_atom)

    # Then solve for the ODE so that they all converge.
    import time
    start = time.time()
    basins = steepest_ascent_rk45(
        points, dens_func, grad_func, beta_spheres, maximas, tol=1e-7, max_ss=0.5, ss_0=0.23
    )
    final = time.time()
    print("Basins", basins)
    print("Length of basins ", len(basins))
    print("Difference ", final - start)
    input("Start Classification")

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
    #     The later method would be to use two nested for-loops, one goes through each maxima
    #     Then goes through each ray, then you would refine to find the IAS intersection.
    #     In the scenario it would cross twice, one can do this scenario twice playing sepcial
    #     attention.
    print("Index to atom ", index_to_atom, points.shape)
    ias_bnds = [OrderedDict() for _ in range(0, numb_maximas)]  # Keys are Points index
    np.set_printoptions(threshold=np.inf)
    for i_maxima in range(0, numb_maximas):
        print("ATom i ", i_maxima)
        print("Starting and Final index", index_to_atom[i_maxima], index_to_atom[i_maxima + 1])
        basins_atoms = basins[index_to_atom[i_maxima]:index_to_atom[i_maxima + 1]]  # Basin of atom
        points_atoms = points[index_to_atom[i_maxima]:index_to_atom[i_maxima + 1]]  # Points of atom
        print("Basins of atom ", basins_atoms)
        numb_rad_pts = numb_rad_to_radial_shell[i_maxima]

        i_ray = 0
        print(index_to_atom[i_maxima], numb_rad_pts)
        print("Number of angular points in this atom", NUMB_RAYS_TO_ATOM[i_maxima])
        for i_ang in range(NUMB_RAYS_TO_ATOM[i_maxima]):
            print("Angular pt j", i_ang)
            print("Number of radial points in this angular pt ", numb_rad_pts[i_ang])

            # Get the basin of the ray
            basins_ray = basins_atoms[i_ray:i_ray + numb_rad_pts[i_ang]]
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
                    r_ubnd = np.linalg.norm(points_atoms[i_ray + index_u_bnd] - maximas[i_maxima])
                    r_lbnd = np.linalg.norm(points_atoms[i_ray + index_u_bnd - 1] - maximas[i_maxima])
                    # Update containers
                    ias_bnds[i_maxima][i_ang] = [r_lbnd, r_ubnd]
                    ias[i_maxima].append(i_ang)
                    print("Radius Lower and Upper bound ", r_lbnd, r_ubnd)
                # else:
                    r"""
                    There are rays for example CH4, where the ray goes from basin 1 to 0 to 1 
                    again, it doesn't make much sense why this is the case, because the ray
                    is very unlikely to do this and should have gone to the other hydrogen. But
                    the density values is incredibly small here.
                    
                    """
                    #     print("IAS With Multiple Intersections")
                    #     print(dens_func(points_atoms[i_ray:i_ray + numb_rad_pts[i_ang]]))
                    #     from chemtools.topology.qtaim import gradient_path_all_pts
                    #     basins = gradient_path_all_pts(
                    #         points_atoms[i_ray:i_ray + numb_rad_pts[i_ang]], grad_func, beta_spheres, i_maxima, maximas,
                    #         t_span=(0, 100), max_step=np.inf, method="LSODA",
                    #         first_step=1e-7
                    #     )

                    # import matplotlib
                    # import matplotlib.pyplot as plt
                    # from mpl_toolkits import mplot3d
                    # matplotlib.use("Qt5Agg")
                    # fig = plt.figure()
                    # ax = plt.axes(projection='3d')
                    # p = points_atoms[i_ray: i_ray + numb_rad_pts[i_ang]]
                    # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
                    # ax.scatter(maximas[:, 0], maximas[:, 1], maximas[:, 2], color="r", s=60)
                    # ax.set_zlabel("Z")
                    # ax.set_ylabel("Y")
                    # plt.show()


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

    # The IAS and OAS were determined, the OAS would be solving for the isosurface value.
    #     The IAS is then just refining the ray, till you find the exact intersection with the
    #     surface.
    # Construct all the points to be solved together.  Honestly, can do these all at once since
    #  the gpu code is fairly fast enough.
    # All of these lists are loop-specific.
    ias_indices = np.array(list(
        itertools.chain.from_iterable(
            [[(i, y, ias_bnds[i][y][0], ias_bnds[i][y][1], max(ss0 / 10.0, bnd_err, ))
              for y in ias[i]] for i in range(numb_maximas)]
        )
    )) # Concatenate all indices together, first index is which index of maxima, then second
    #  index is index of angular point, third index is the lower bound radius and fourth
    #  index is the upper-bound radius, fifth index is the step-size.
    print("ias indices", ias_indices)
    input('Start Solving IAS')
    while len(ias_indices) != 0:
        # Construct New Points
        points = []
        numb_pts_per_ray = []
        for (i_maxima, i_ang, l_bnd, u_bnd, ss) in ias_indices:
            print((i_maxima, i_ang, l_bnd, u_bnd, ss))
            ray = maximas[int(i_maxima)] + np.arange(l_bnd, u_bnd + ss, ss)[:, None] * angular_pts[int(i_ang), :]
            points.append(ray)
            numb_pts_per_ray.append(len(ray))
        points = np.vstack(points)

        # Solve for basins
        start = time.time()
        basins = steepest_ascent_rk45(
            points, dens_func, grad_func, beta_spheres, maximas, tol=1e-7, max_ss=0.5, ss_0=0.23
        )
        final = time.time()
        print("Basins", basins)
        print("Difference ", final - start)

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
                basin_ias[int(i_maxima)].append(basins[i_switch])
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
        input("Next ray")
        # Remove converged indices
        ias_indices = np.delete(ias_indices, converge_indices, axis=0)


    # Solve for multiple intersections
    # TODO


    # Solve OAS Points
    for i_maxima in range(numb_maximas):
        maxima = maximas[i_maxima]
        radial = radial_grid[i_maxima]
        for i_oas in oas[i_maxima]:
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

    return SurfaceQTAIM(r_func, angular, maximas, oas, ias, basin_ias, refined_ang)
