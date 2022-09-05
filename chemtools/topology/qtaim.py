
import numpy as np

from scipy.integrate import solve_ivp
from scipy.interpolate import LSQSphereBivariateSpline, SmoothSphereBivariateSpline
from scipy.spatial import ConvexHull, Delaunay, Voronoi
from scipy.spatial.distance import cdist

from scipy.sparse import lil_matrix

from grid.atomgrid import AtomGrid
from grid.cubic import UniformGrid, _HyperRectangleGrid
from grid.lebedev import AngularGrid
from grid.utils import convert_cart_to_sph

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits import mplot3d

r"""
QTAIM Algorithm 1 (General Grid) with Yu-Trinkle
-------------------------------------------
1) Sort based on density values
2) Calculate the Voronoi Diagram 
3) Find and store all neighbors of each point using the Voronoi Diagram, may miss some points on boundary.
4) Go through each point
    i) See if the neighbors of that point are assigned.
    ii) If no neighbours are assigned and it has neighbors then it is a maxima point
    iii) If no neighbors are assigned and no neighbor information is found, then assing
        its weight based on the average of its closest three points.
    iv) If all neighbors that were assigned, are assigned to a single basin, assign this to that basin.
    v) If some neighbors are assigned to different basins, then this point is a watershed point. 
        - Solve using scipy ode
        - Yu-trinkle algorithm.



QTAIM ALGORITHM 2
-----------------
Assumes if you're within the beta-sphere, then we're done.

"""


class _BasinContainer(object):
    __slots__ = ["basin", "numb_basins_found", "num_centers"]

    def __init__(self, num_pts, num_centers=None):
        self.numb_basins_found = 0
        self.num_centers = 1 if num_centers is None else num_centers
        self.basin = lil_matrix((num_pts, self.num_centers), dtype=np.float64)

    def __getitem__(self, index):
        # Get the basin values for the `index`th point based on the maximum.
        # If the basin value returned is -1.0, it means it wasn't found yet.
        arr = self.basin.getrow(index).toarray()[0]  # [0] converts ndim to 1.
        if any(x != 0.0 for x in arr):
            # Plus one because basin values defined/starts at one.
            return arr.argmax() + 1
        return -1.0

    def __setitem__(self, index, basin):
        # This is the case of using Henklemenb/Bader's method on watershed points.
        if isinstance(basin, (int, float, np.float64, np.int64)):
            if basin > 0.0:
                # Assign to the `index` point to basin number with one.
                self.basin[index, int(basin) - 1] = 1.0
            else:
                raise ValueError(f"Basin value {basin} to the point {index} should be "
                                 f"greater than zero.")
        # This is the case of when you use Yu-trinkle algorithm on watershed points.
        elif isinstance(basin, (list, np.ndarray)):
            self.basin[index, :] = basin
        else:
            raise TypeError(f"Basin {type(basin)} should be a number of a float/list/array.")

    def get_basins_from_indices(self, indices):
        # Get the basins from the indices of the points removing all zero elements.
        # FIXME : Add error exception if hte indices from watershed backtracing are outside the
        # grid.
        # This removes -1 because the __get_item__ returns -1.0 if it doesn't have any neighbors.
        return {self.__getitem__(i) for i in indices} - {-1.0}

    def get_basin_weights_of_points(self, indices):
        # Given a set of point indices, get the basin weights for all weights
        return self.basin[indices].toarray()

    def increase_numb_atoms(self):
        self.numb_basins_found += 1
        # If the numb basins found is greater than num_centers (specified by user)
        # then it resizes the sparse array.
        if self.numb_basins_found > self.num_centers:
            self.num_centers += 1
            shape = self.basin.shape
            # Resize is better than reshape as it changes it in-place.
            self.basin.resize((shape[0], self.num_centers))


def _get_area_of_coplanar_polygon(points):
    # see math stackexchange: how-do-you-calculate-the-area-of-a-2d-polygon-in-3d
    # points (M, 3) array, assumes the points all lie on a plane, i.e. coplanar
    # this assumes that the point are ordered adjacently.
    num_verts = points.shape[0]
    center = np.sum(points, axis=0) / num_verts  # get the center of points
    area = 0
    for i in range(num_verts):
        v_i_plues_one = points[0] if i == num_verts - 1 else points[i + 1]
        area += np.linalg.norm(np.cross(
            (points[i] - center), (v_i_plues_one - center)
        )) / 2.0
    return area


def _get_area_of_voronoi_ridge(i_pt, i_nbh, index_to_voronoi_ridge, voronoi):
    # index_to_voronoi_ridge list of lists
    #  find the row index r_{ij} in voronoi.ridge_points that contains (i, j)
    ith_voronoi_ridges = index_to_voronoi_ridge[i_pt]
    i_nbh_voronoi_ridges = index_to_voronoi_ridge[i_nbh]
    ridge_index = (set(ith_voronoi_ridges) & set(i_nbh_voronoi_ridges)) # Take intersection
    assert len(ridge_index) == 1
    ridge_index = ridge_index.pop()
    # Get the voronoi vertices  via : voronoi.vertices[delaunay.ridge_vertices[r_{ij}]].
    #   voronoi.ridge_vertices[r_{ij}] makes sure it doesn't have -1 in it.
    ridge_vertices = voronoi.vertices[voronoi.ridge_vertices[ridge_index]]
    print("ridge vertices", ridge_vertices, "voronoi ridge vertices", voronoi.ridge_vertices[ridge_index])
    assert -1 not in voronoi.ridge_vertices[ridge_index]
    # Get the area defined by the polygon of size 4, this assumes the polygon is coplanar,
    # i.e. lies on a plane. For a calculation of the formula see:
    # stack exchange article: how-do-you-calculate-the-area-of-a-2d-polygon-in-3d
    # return _get_area_of_coplanar_polygon(ridge_vertices)
    if len(ridge_vertices) <= 3:
        return _get_area_of_coplanar_polygon(ridge_vertices)
    return ConvexHull(ridge_vertices, qhull_options="QJ").area / 2.0


def _assign_weight_yu_trinkle_voronoi(index, basin_cont, density_vals, voronoi,
                              neighbors_index, index_to_voronoi_ridge):
    total_fraction = 0.0
    weights = np.zeros((basin_cont.num_centers,))
    weights_nbhs = basin_cont.get_basin_weights_of_points(neighbors_index)
    print("basin weights of neighbors", weights_nbhs)
    # Go through each neighbour X` of the `index`th point X.
    for k, nbh_index in enumerate(neighbors_index):
        density_diff = density_vals[nbh_index] - density_vals[index]
        print("Nbh Index ", nbh_index, "Nbh Point ", voronoi.points[nbh_index], "Density Diff ", density_diff, "Densities ",  density_vals[nbh_index], density_vals[index])
        # Only look at neighbours X` whose density values is greater
        if 0.0 < density_diff:
            # Calc flux-probability: J_{X, X`} = a_{X, X`} (p_{X`} - p_{X}) / l_{X, X`}
            length = np.linalg.norm(voronoi.points[index] - voronoi.points[nbh_index])
            area = _get_area_of_voronoi_ridge(index, nbh_index, index_to_voronoi_ridge, voronoi)
            flux_probability = area * density_diff / length
            print("Flux Probability ", flux_probability, "length ", length)
            total_fraction += flux_probability
            print("Weight of neighbor", weights_nbhs[k], "area", area)
            # Calculate w^A(X) = \sum_{X`} J_{X, X`} w^A(X`)
            weights += flux_probability * weights_nbhs[k]

    weights /= total_fraction
    assert total_fraction != 0.0, "The neighbors most likely have the same density values, finer grid " \
                                  "which can avoid points with identical neighbors is recommended."
    print(total_fraction, weights)
    return weights


def _assign_weight_yu_trinkle_cubic(index, basin_cont, density_vals, closest_nbh_indices,
                                    uniform_grid, areas, grad_func=None):
    total_fraction = 0.0
    weights = np.zeros((basin_cont.num_centers,))
    weights_nbhs = basin_cont.get_basin_weights_of_points(closest_nbh_indices)

    print("basin weights of neighbors", weights_nbhs)
    # Go through each neighbour X` of the `index`th point X.
    for k, nbh_index in enumerate(closest_nbh_indices):
        density_diff = density_vals[nbh_index] - density_vals[index]
        print("Nbh Index ", nbh_index, "Nbh Point ", uniform_grid.points[nbh_index], "Density diff",
              density_diff)
        # Only look at neighbours X` whose density values is greater
        if 0.0 < density_diff:
            # Calc flux-probability: J_{X, X`} = a_{X, X`} (p_{X`} - p_{X}) / l_{X, X`}
            length = np.linalg.norm(uniform_grid.points[index] - uniform_grid.points[nbh_index])
            # area = _get_area_of_voronoi_ridge(index, nbh_index, index_to_voronoi_ridge, voronoi)
            area = areas[k]
            print("area", area, "length ", length)
            if grad_func is None:
                flux_probability = area * density_diff / length
            else:
                # calculate normal to the Voronoi Facet/Boundary
                normal = uniform_grid.points[nbh_index] - uniform_grid.points[index]
                normal /= length
                midpoint = (uniform_grid.points[index] + uniform_grid.points[nbh_index]) / 2.0
                flux_probability = area * normal.dot(grad_func(np.array([midpoint]))[0])
            print("Flux Probability ", flux_probability)
            total_fraction += flux_probability
            print("Weight of neighbor", weights_nbhs[k])
            # Calculate w^A(X) = \sum_{X`} J_{X, X`} w^A(X`)
            weights += flux_probability * weights_nbhs[k]

    weights /= total_fraction
    return weights


def _get_neighbor_and_ridge_information_from_voronoi(voronoi):
    r"""
     Voronoi data structure doesn't give a good data-structure to figure out which neighbors are
       between a point in `points`.  So here, we convert that based on the attribute
       voronoi.ridge_points (Indices of the points between which each Voronoi ridge lies.).
      Neighbours[i_pt]=[i_1, ..., i_n] gets the neighbors indices of the i_pt.
      index_to_voronoi_ridge[i_pt] = [r_1, .., r_N] gives the voronoi ridges index r_k.
    """
    neighbors_indices = [[] for _ in range(0, voronoi.points.shape[0])]
    index_to_voronoi_ridge = [[] for _ in range(0, voronoi.points.shape[0])]
    for i_ridge, (x, y) in enumerate(voronoi.ridge_points):
        neighbors_indices[x] += [y]
        neighbors_indices[y] += [x]
        index_to_voronoi_ridge[x] += [i_ridge]
        index_to_voronoi_ridge[y] += [i_ridge]
    return neighbors_indices, index_to_voronoi_ridge


def _assign_weight_average_neighbors(index, basin_cont, voronoi, original_num_pts):
    # QHull voronoi algorithm has trouble with corner points and finding neighbors
    # stackoverflow: questions/25754145/scipy-voronoi-3d-not-all-ridge-points-are-shown .
    # Solution:find the three closest points to this point and take the average of the
    # weights to define the weight of this point.  Three was chosen because a corner
    # point in a rectangular grid has three neighbors.
    print("QHull/Voronoi couldn't find neighbors, manually find average.")
    distance = cdist(voronoi.points[index:index + 1], voronoi.points)[0]
    min_index = distance.argsort()
    min_index = np.delete(min_index, min_index >= original_num_pts)
    min_index = min_index[1:4]  # ignore first point cause it is itself
    basin_wghts = basin_cont.get_basin_weights_of_points(min_index)
    print("Average ", basin_wghts)
    weights = np.average(basin_wghts, axis=0)  # Take average
    print("Average Weights", weights)
    return weights


def voronoi_volumes(voronoi):
    # Given Voronoi, this calculates the volume of each Voronoi region.
    vol = np.zeros(voronoi.npoints)
    for i, reg_num in enumerate(voronoi.point_region):
        indices = voronoi.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(voronoi.vertices[indices]).volume
    return vol


def close_neighbors_step():
    # Given a point coordinate in 3D grid (i, j, k). Adding these give the coordinates of its
    # close neighbors.
    return np.array([
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
        [0, 1, 0],
        [1, 0, 0]
    ])

def close_diagonal_neighbors_step():
    # Given a point coordinate in 3D grid (i, j, k). Adding these give the coordinates of its
    # close diagonal neighbors.
    return np.array([
        [-1, -1, 0],
        [-1, 0, -1],
        [-1, 0, 1],
        [-1, 1, 0],
        [0, -1, -1],
        [0, -1, 1],
        [0, 1, -1],
        [0, 1, 1],
        [1, -1, 0],
        [1, 0, -1],
        [1, 0, 1],
        [1, 1, 0],
    ])

def diagonal_neighbors_step():
    # Given a point coordinate in 3D grid (i, j, k). Adding these give the coordinates of its
    # diagonal neighbors.
    return np.array([
        [-1, -1, -1],
        [-1, -1, 1],
        [-1, 1, -1],
        [-1, 1, 1],
        [1, -1, -1],
        [1, -1, 1],
        [1, 1, -1],
        [1, 1, 1]
    ])


def _points_on_bounding_box(points, step_size=0.25, extension=0.01):
    r"""Get the points on the surface of a bounding box of a specified grid."""
    # Place bounding box over the points
    l_bnd = np.min(points, axis=0)
    u_bnd = np.max(points, axis=0)

    # Compute the required number of points along x, y, and z axis
    shape = np.ceil((u_bnd - l_bnd + 2.0 * extension) / step_size)
    axes = np.eye(3) * step_size

    # construct x-y plane bottom and top
    coords = np.array(np.meshgrid(np.arange(shape[0] - 1), np.arange(shape[1] - 1), [0.]))
    coords = np.swapaxes(coords, 1, 2)
    coords = coords.reshape(3, -1)
    new_pts_bottom = coords.T.dot(axes) + np.array([l_bnd[0], l_bnd[1], l_bnd[2] - extension])
    points = np.vstack((points, new_pts_bottom))
    new_pts_up = coords.T.dot(axes) + np.array([l_bnd[0], l_bnd[1], u_bnd[2] + extension])
    points = np.vstack((points, new_pts_up))


    # construct y-z plane left and right
    coords = np.array(np.meshgrid([0.], np.arange(shape[1] - 1), np.arange(shape[2] - 1)))
    coords = np.swapaxes(coords, 1, 2)
    coords = coords.reshape(3, -1)
    new_pts_left = coords.T.dot(axes) + np.array([l_bnd[0] - extension, l_bnd[1], l_bnd[2]])
    points = np.vstack((points, new_pts_left))
    new_pts_right = coords.T.dot(axes) + np.array([u_bnd[0] + extension, l_bnd[1], l_bnd[2]])
    points = np.vstack((points, new_pts_right))



    # construct x-z plane towards and back
    coords = np.array(np.meshgrid(np.arange(shape[0] - 1), [0.], np.arange(shape[2] - 1)))
    coords = np.swapaxes(coords, 1, 2)
    coords = coords.reshape(3, -1)
    new_pts_down = coords.T.dot(axes) + np.array([l_bnd[0], l_bnd[1] - extension, l_bnd[2]])
    points = np.vstack((points, new_pts_down))
    new_pts_up = coords.T.dot(axes) + np.array([l_bnd[0], u_bnd[1] + extension, l_bnd[2]])
    points = np.vstack((points, new_pts_up))

    unique_pts, indices = np.unique(points, return_index=True, axis=0)
    assert unique_pts.shape == points.shape, "Bounding box is not unique."
    return points


def _get_neighbor_indices_for_cubic_grid(index, type, uniform_grid, return_area=False):
    coord = uniform_grid.index_to_coordinates(index)
    print(coord)
    if type == "closest-neighbors":
        nbh_coords = close_neighbors_step() + coord
    elif type == "all-neighbors":
        nbh_coords = np.vstack(
            (close_neighbors_step(), diagonal_neighbors_step(), close_diagonal_neighbors_step())
        ) + coord
    else:
        raise ValueError(f"Could not recognize {type}.")

    # -1 or num_pts means neighbors doesn't exist.
    nbh_coords = np.delete(nbh_coords, np.any(nbh_coords == -1, axis=1), axis=0)
    nbh_coords = np.delete(nbh_coords, np.any(nbh_coords == uniform_grid.shape, axis=1), axis=0)
    closest_nbh_indices = [uniform_grid.coordinates_to_index(x) for x in nbh_coords]
    if return_area:
        if type == "closest-neighbors":
            ss = np.array([np.linalg.norm(axis) for axis in uniform_grid.axes])
            ss = (1 - np.abs(close_neighbors_step())) * ss
            ss[ss == 0] = 1
            return closest_nbh_indices, np.prod(ss, axis=1)
        else:
            raise ValueError(f"`return_area` is true only when type == 'closest-neighbors'.")
    return closest_nbh_indices


def qtaim(grid_pts, density_vals, grad_func=None, num_centers=None, method="yu-trinkle",
          remove_duplicates=True, bounding_box=True):
    r"""
    Define what a watershed point is.

    Qtaim fails when density_Vals is very small, like  a peak Gaussian with alpha=500
    
    """
    # Assert method values
    if not isinstance(grid_pts, (_HyperRectangleGrid, np.ndarray)):
        raise TypeError(f"Points should either be a numpy array or a UniformGrid object.")
    is_cubic_grid = isinstance(grid_pts, _HyperRectangleGrid)
    points = grid_pts if not is_cubic_grid else grid_pts.points

    # if remove_duplicates:
    #     points, indices = np.unique(points, return_index=True, axis=0)
    #     density_vals = density_vals[indices]

    original_num_pts = points.shape[0]
    if not is_cubic_grid:
        if bounding_box:
            points = _points_on_bounding_box (points, extension=0.25, step_size=0.1)
        voronoi = Voronoi(points=points, qhull_options="Qbb Qc Qz")
        assert np.all(
            np.abs(voronoi.points - points) < 1e-8), "Voronoi points should be the same as points."
        # neighbors_indices: mapping of point index to the neighbors using Voronoi diagram
        # neighbors_indices: mapping of point index to the Voronoi ridges it's part of.
        neighbors_indices, index_to_voronoi_ridge = _get_neighbor_and_ridge_information_from_voronoi(
            voronoi
        )

    # num_centers speed calculations up dramatically due to sparsity structure
    num_centers = 1 if num_centers is None else num_centers
    maxima_indices = []  # Indices of the centers.
    watershed_indices = []  # indices of points that are on the boundary of the basin.
    basin_cont = _BasinContainer(original_num_pts, num_centers)
    # TODO: When density values are sorted, maybe group them based on similar values and pick the ones closest to the points.
    sorted_density_indices = np.argsort(density_vals)[::-1]

    # Go through each point with the highest density values to the smallest.
    for i, index in enumerate(sorted_density_indices):
        print("Index ", index, " Point ", points[index], " Density value ", density_vals[index])
        if not is_cubic_grid:
            # Get the closest neighbor indices and remove those that are part of the bounding box.
            closest_nbh_indices = neighbors_indices[index]
            closest_nbh_indices = [i for i in closest_nbh_indices if i < original_num_pts]
            print("Voronoi ", voronoi.regions[voronoi.point_region[index]])
        else:
            closest_nbh_indices, areas = _get_neighbor_indices_for_cubic_grid(
                index, "closest-neighbors", grid_pts, return_area=True
            )
        # Get the basin-values of the closest points.
        basin_vals = basin_cont.get_basins_from_indices(closest_nbh_indices)
        print("Neighbours Indices ", closest_nbh_indices)
        print("Basin of neighbours ", basin_vals)
        # Closest neighbours were all not assigned yet a basin, means it is a maximum.
        if len(basin_vals) == 0:
            if len(closest_nbh_indices) > 0:  # If the  neighbors were found
                found_maxima = not is_cubic_grid
                if is_cubic_grid and not found_maxima:
                    # Check all neighbors rather than close neighbors. Voronoi already checks all.
                    all_nbh_indices = _get_neighbor_indices_for_cubic_grid(
                        index, "all-neighbors", grid_pts
                    )
                    print("All neighbors", all_nbh_indices)
                    basin_vals = basin_cont.get_basins_from_indices(all_nbh_indices)
                    found_maxima = len(basin_vals) == 0

                if found_maxima:
                    print("Maximum found")
                    maxima_indices.append(index)
                    basin_cont.increase_numb_atoms()
                    basin_cont[index] = basin_cont.numb_basins_found
                    print("Number of basins founds ", basin_cont.numb_basins_found,
                          " Number Centers Specified At Beginning ", basin_cont.num_centers)
                elif len(basin_vals) == 1:
                    # (Cubic grid only)Diagonal element probably assigned, so assign it to that.
                    print("Point assigned to the basin of the neighbors")
                    basin_cont[index] = basin_vals.pop()
                else:
                    # TODO this case do watershed
                    # Most likely occurs due to exact density vals and sorting is unordered.
                    raise NotImplementedError("TODO")
            elif not is_cubic_grid:
                # Assign weight based on average of its neighbors due to problems with QHull.
                weights = _assign_weight_average_neighbors(index, basin_cont, voronoi, original_num_pts)
                # assert the weights are not all zero
                if np.all(weights == 0.0):
                    print("weights ", weights)
                    raise RuntimeError("Weights are all zero")
                basin_cont[index] = weights

        # All neighbours were assigned to a single basin, assign this point to that basin.
        elif len(basin_vals) == 1:
            print("Point assigned to the basin of the neighbors")
            basin_cont[index] = basin_vals.pop()
        else:
            # It is assigned to multiple basins.
            print("Found watershed point")
            watershed_indices.append(index)
            # Consider the case it is a critical point, how do you check for this?
            # check for the gradient.  Consider that one could assign a special criteria for this.
            if not is_cubic_grid:
                neighbor_index = neighbors_indices[index]
                neighbor_index = [i for i in neighbor_index if i < original_num_pts]
                weights = _assign_weight_yu_trinkle_voronoi(
                    index, basin_cont, density_vals, voronoi, neighbor_index,
                    index_to_voronoi_ridge
                )
            else:
                weights = _assign_weight_yu_trinkle_cubic(
                    index, basin_cont, density_vals, closest_nbh_indices, grid_pts, areas,
                    grad_func
                )

            print(weights)
            print("Sum of weights ", np.sum(weights))
            basin_cont[index] = weights
            if np.abs(np.sum(weights) - 1.0) > 1e-10:
                raise RuntimeError(
                    f"The Weights {weights} did not sum {np.sum(weights)} up to one."
                )
        # if index == 729:
        #     input("sd")
        print("")

    # TODO: Update watershed indices
    # Calculate Voronoi volumes
    volume = np.prod(np.array([np.linalg.norm(axis) for axis in grid_pts.axes])) if is_cubic_grid \
        else voronoi_volumes(voronoi)[:original_num_pts]
    return {"basin_cont" : basin_cont.basin.tocsc(), "maxima_indices" : maxima_indices,
            "watershed_indices" : watershed_indices, "voronoi_volumes" : volume}


def gradient(pt, grad_func):
    grad= grad_func(np.array([pt]))[0].T
    return grad / np.linalg.norm(grad)


def _project_pt_on_line(pt, ray_pt1, ray_pt2):
    r"""Project point onto a line defined by maxima, point on sphere and radial points."""
    ap = pt - ray_pt1
    ab = ray_pt2 - ray_pt1
    return ray_pt1 + np.dot(ap, ab) / np.dot(ab, ab) * ab


class SurfaceQTAIM():
    def __init__(self, r_func, rad_grids, angular_degs, maximas, oas, ias, basins_ias,
                 refined_ang=None):
        self._r_func = r_func
        self._maximas = maximas
        self._rad_grids = rad_grids
        self._angular_degs = angular_degs
        self._oas = oas
        self._ias = ias
        self._basins_ias = basins_ias
        self._refined_ang = refined_ang

    @property
    def r_func(self):
        # List[M, np.ndarray(N_i,)] ::math:`r_j(\theta_i, \phi_i)` for each M basins and N_i angular points.
        return self._r_func

    @property
    def oas(self):
        # List[List[int]] : First list is over basins, second over indices of points of outeratomic surface.
        return self._oas

    @property
    def ias(self):
        # List[List[int]] : First list is over basins, second over indices of points of interatomic surface.
        return self._ias

    @property
    def maximas(self):
        # ndarray(M, 3) : The maxima of each basin.
        return self._maximas

    @property
    def angular_degs(self):
        # int or List[int] : List of angular degrees over each basin.
        return self._angular_degs

    @property
    def rad_grids(self):
        # List[OneDGrids] : List of M OneDGrids for integrating over radial component in [0, \inty).
        return self._rad_grids

    @property
    def basins_ias(self):
        return self._basins_ias

    @property
    def refined_ang(self):
        # List[M, np.ndarray(N_i, 2)] : Additional Angular points to append to the angular grid
        return self._refined_ang

    def generate_angular_grid_of_basin(self, i_basin):
        # Note this doesn't include the extra angular points generated by refinement.
        deg = self.angular_degs
        deg = deg[i_basin] if isinstance(deg, list) else deg
        return AngularGrid(degree=deg)

    def generate_angular_pts_of_basin(self, i_basin):
        angular_grid = self.generate_angular_grid_of_basin(i_basin)
        points = angular_grid.points
        if self.refined_ang is not None:
            points = np.vstack((points, self.refined_ang[i_basin]))
        return points

    def get_atom_grid_over_basin(self, i_basin):
        # integrate over a basin.
        deg = self.angular_degs
        deg = deg[i_basin] if isinstance(deg, list) else deg
        atom_grid = AtomGrid(self.rad_grids[i_basin], degrees=[deg], center=self.maximas[i_basin])

        # Go through each spherical point and get the r(theta_i, phi_i) limit
        for i_sph in range(atom_grid.get_shell_grid(0).size):
            r_limit = self.r_func[i_basin, i_sph]
            # Go through each radial point and see if it is larger than the limit.
            for i_rad in range(atom_grid.rgrid.size):
                if atom_grid.rgrid.points[i_rad] > r_limit:
                    # Find where (r_{ij}, theta_i, phi_i) is and set the weights to zero.
                    i_start, i_final = atom_grid.indices[i_rad], atom_grid.indices[i_rad + 1]
                    atom_grid.weights[i_start + i_sph] = 0.0
        # atom_grid.weights[inequality] = 0.0
        return atom_grid

    def generate_pts_on_surface(self, i_basin):
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        return self.maximas[i_basin] + self.r_func[i_basin][:, None] * sph_pts

    def get_ias_of_basin(self, i_basin):
        ias = self.ias[i_basin]
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        return self.maximas[i_basin] + self.r_func[i_basin][ias, None] * sph_pts[ias]

    def get_oas_of_basin(self, i_basin):
        oas = self.oas[i_basin]
        sph_pts = self.generate_angular_pts_of_basin(i_basin)
        return self.maximas[i_basin] + self.r_func[i_basin][oas, None] * sph_pts[oas]

    def interpolate_radial_func(self, method="smooth", ias=False, oas=False):
        # if method not in ["smooth", ]
        if ias and oas:
           raise ValueError(f"Both {ias} and {oas} cannot be true.")
        if ias:
            pass



def get_closest_points_between_ias_and_oas(
    ias: list, oas: int, angular_pts: np.ndarray, r_func_max: np.ndarray, maxima: np.ndarray
):
    r"""

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
    new_ang_pts = np.zeros((len(ias_bnd) * 7, 3))  # 7 points per ias boundary are added.
    for i_ias, pt_ias in enumerate(ias_bnd):
        two_indices = dist_mat[i_ias].argsort()[:2]
        pt1, pt2 = oas_bnd[two_indices[0]], oas_bnd[two_indices[1]]

        print(i_ias, pt_ias, pt1, pt2)

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
        direction = new_pts - maxima
        t = np.linalg.norm(direction, axis=1)
        direction = direction / t[:, None]
        print(direction)

        new_ang_pts[i_ias * 7:(i_ias + 1) * 7] = direction

        # import matplotlib
        # import matplotlib.pyplot as plt
        # from mpl_toolkits import mplot3d
        # matplotlib.use("Qt5Agg")
        # fig = plt.figure()
        # ax = plt.axes(projection='3d')
        # p = np.array([pt1, pt2, pt_ias])
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
        # p = ias_bnd
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
        # p = oas_bnd
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
        # p = maxima + t[:, None] * direction
        # ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="y", s=60)
        # plt.show()
    return new_ang_pts


def gradient_path(pt, grad_func, t_span=(0, 1000), method="LSODA", max_step=100, t_inc=200,
                  max_tries=10, first_step=1e-3, beta_sphere_maxima=-np.inf, maxima=None):
    # TODO: If the density value is low, gradient low and trying ODE did not move much, then
    #  an option is to turn max_step tp np.inf and change t_span to 10,000.
    is_converged = False
    y0 = pt.copy()
    numb_times = 0
    while not is_converged and numb_times < max_tries:
        sol = solve_ivp(
            lambda t, x: grad_func(np.array([x]))[0].T,
            y0=y0,
            t_span=t_span,
            method=method,
            max_step=max_step,
            first_step=first_step,
        )
        # print(sol)
        assert sol["success"], "ODE was not successful."
        convergence = np.abs(sol["y"][:, -2] - sol["y"][:, -1])
        if np.all(convergence < 0.01):
            return sol["y"][:, -1]
        elif beta_sphere_maxima != -np.inf:
            # If the point converged to the beta sphere of the maxima, then we stop.
            radial = np.linalg.norm(sol["y"][:, -1] - maxima)
            if radial < beta_sphere_maxima:
                print("Close to beta sphere.")
                return sol["y"][:, -1]
        else:
            # This isn't good for finding isosurfaces, because it would keep on going longer than expected.
            # Also I can do the beta sphere trick here for convegence rather than going all the towards.
            print(sol["y"][:, -1], t_span)
            t_span = (t_span[1], t_span[1] + t_inc)
            y0 = sol["y"][:, -1]
            numb_times += 1
    if numb_times == max_tries:
        raise RuntimeError(f"No convergence in gradient path pt {pt},"
                           f" solution {sol['y'][:, -1]}, t_span {t_span}")


def solve_for_isosurface_pt(index_iso, rad_pts, maxima, cart_sphere_pt, density_func,
                            iso_val, iso_err):
    # Given a series of points based on a maxima defined by angles `cart_sphere_pt` with
    #  radial pts `rad_pts`.   The `index_iso` tells us where on these points to construct another
    #  refined grid from finding l_bnd and u_bnd.  This assumes the lower bound and upper bound
    #  contains the isosurface point.  This point is solved using a root-finding algorithm to
    #  solve for the root of the density function.
    print(rad_pts, index_iso)
    l_bnd = rad_pts[index_iso - 1] if index_iso >= 0 else rad_pts[index_iso] / 2.0
    u_bnd = rad_pts[index_iso + 1] if index_iso + 1 < len(rad_pts) else rad_pts[index_iso] * 2.0
    dens_l_bnd = density_func(np.array([maxima + l_bnd * cart_sphere_pt]))
    dens_u_bnd = density_func(np.array([maxima + u_bnd * cart_sphere_pt]))
    if iso_val < dens_u_bnd or dens_l_bnd < iso_val:
        raise ValueError(f"Radial grid {l_bnd, u_bnd} did not bound {dens_l_bnd, dens_u_bnd} "
                         f"the isosurface value {iso_val}. Use larger radial grid.")

    # Use Root-finding algorithm to find the isosurface point.
    root_func = lambda t: density_func(np.array([maxima + t * cart_sphere_pt]))[0] - iso_val
    from scipy.optimize import root_scalar
    sol = root_scalar(root_func, method="toms748", bracket=(l_bnd, u_bnd), xtol=iso_err)
    print(sol)
    assert sol.converged, f"Root function did not converge {sol}."
    bnd_pt = maxima + sol.root * cart_sphere_pt
    print(bnd_pt, density_func(np.array([bnd_pt])))
    return bnd_pt


def solve_for_basin_bnd_pt(
    dens_cutoff, maxima, radial, cart_sphere_pt, density_func, grad_func, other_maximas, bnd_err,
    ss_watershed, iso_val, beta_sphere_maxima, beta_sphere_others
):
    # Construct the ray and compute its density values based on a maxima defined by angles
    #  `cart_sphere_pt` with radial pts `rad_pts`.    It goes through each point on the ray
    #   if the ray density value is greater than dens_cutoff, then it is likely this ray
    #   tends towards infinity and has no basin boundary.  If density value is larger, then
    #   it solves for the gradient path via solving gradient ode.  If this ode solution,
    #   is close to other basins, then we found when it switched basins.  Then we take
    #  the two points where it switches basin and compute the distance, if this distance
    #  is less than `bnd_err`, then we take the midpoint to be the boundary point on the ray
    #  that intersects the ias.  If not, then we construct a new ray with different l_bnd
    #  and u_bnd and reduce the step-size further and repeat this process.
    rad_pts = radial.copy()
    ss_ray = ss_watershed * 100   # Stay with a coarse ray then refine further.
    index_iso = None  # Needed to refine if the ray tends towards infinity.
    bnd_pt = None  # Boundary or Isosurface Point
    is_ray_to_inf = False  # Does this ray instead go towards infinity

    # TODO: Consider increase upper bound if it fails.
    found_watershed_on_ray = False
    basin_id = None
    while not found_watershed_on_ray:
        ray = maxima + rad_pts[:, None] * cart_sphere_pt
        ray_density = density_func(ray)
        print("Start of Ray ", ray[0], " Cartesian pt of Sphere ", cart_sphere_pt, "Final Ray Pt: ",
              ray[-1])

        # Compute the boundary point
        for i_ray, pt in enumerate(ray):
            if ray_density[i_ray] > dens_cutoff:
                # Multiply the integration span by the radial component so it scales.
                # Make the first step proportional to the step-size of the ray
                y_val = gradient_path(pt, grad_func,
                                      t_span=(0, max(1000 * rad_pts[i_ray], 75)),
                                      max_step=50,
                                      beta_sphere_maxima=beta_sphere_maxima,
                                      maxima=maxima)#, first_step=ss_ray / 10)

                print("Pt ", pt, "Y ", y_val, "Maxima ", maxima)
                # If the next point is closer to other maximas or in beta spheres
                dist_maxima = np.linalg.norm(y_val - other_maximas, axis=1)
                if np.any(dist_maxima < 1e-1) or np.any(dist_maxima < beta_sphere_others):
                    print("Close to other maxima")
                    # If dist between the basin switch is less than boundary err
                    dist = np.linalg.norm(ray[i_ray] - ray[i_ray - 1])
                    if dist < bnd_err:
                        # Take the midpoint to be the boundary point.
                        found_watershed_on_ray = True
                        bnd_pt = (ray[i_ray] + ray[i_ray - 1]) / 2.0
                        basin_id = np.where(dist_maxima < 1e-1)[0][0] + 1  # basins starts at 1
                        print("Found the Watershed point ", bnd_pt, basin_id)
                    else:
                        # Refine Grid Further
                        l_bnd = np.linalg.norm(ray[i_ray - 1] - maxima) if i_ray != 0 else 1e-3
                        u_bnd = np.linalg.norm(ray[i_ray] - maxima)
                        ss_ray /= 10.0
                        rad_pts = np.arange(l_bnd, u_bnd + ss_ray, ss_ray)
                        print("Refine the ray further with l_bnd, u_bnd, ss: ", l_bnd, u_bnd,
                              ss_ray)
                    break  # break out of this ray and either quit or do refined ray.
            else:
                # The density values became less than the isosurface cut-off
                print("Stop: Density value is less than isosurface cut-off.")
                is_ray_to_inf = True
                index_iso = i_ray
                break  # break out of ray loop
            is_ray_to_inf = True if i_ray == len(ray) - 1 else is_ray_to_inf

        if is_ray_to_inf:
            index_iso = np.argsort(np.abs(ray_density - iso_val))[0]
            print("Ray to infinity with index ", index_iso)
            break
    return bnd_pt, is_ray_to_inf, index_iso, found_watershed_on_ray, basin_id


def qtaim_line_search(rgrids, angular, centers, density_func, grad_func, dens_cutoff=1e-5, iso_val=0.001,
                      ss_watershed=1e-3, bnd_err=1e-4, iso_err=1e-6, beta_sphere=None,
                      optimize_centers=True, refine=False):
    r"""
    One approach is to initially only test the endpoint of the ray.  If this ray ends up at
    the same basin, then it is a ray to infinity. If not, then this ray crosses the IAS.  This
    assumes the ray is far away.  The downside of this approach is the how slow the convergence
    is for points that are far away.  The other downside is for very complicated basins, this
    approach might not work.

    - refine is not useful for integration purposes other than using the interpolate function.

    Parameters
    ----------
    angular : int or np.ndarray[M, 2]
        Either integer specifying the degree to construct angular/Lebedev grid around each maxima
        or a array of spherical coordinates :math:`(\theta, \phi)`.

    beta_sphere : list[float]
        List of size `M` of radius of the sphere centered at each maxima. It avoids backtracing
        of points within the circle. If None is provided, then it doesn't use beta sphere
        for that maxima.

    refine : (bool, int)
        Either False or an integer indicating which center is being runned on.

    Returns
    -------

    Notes
    -----
    - Optimizing centers would affect the beta-sphere radius TODO.
    - It is possible for a Ray to intersect the zero-flux surface but this algorihtm will
        classify it as a ray to infinity because the points on the other side of the basin have
        density values so small that the ode doesn't converge to the maxima of the other basin.
        In this scenario it might be worthwhile to have a denser radial grid with less points
        away from infinity or have a smaller density cut-off.  Alternative for the developer,
        is to implement highly accurate ode solver at the expense of computation time.
    - Radeau helps for stiff problems far from the Nucleus, whereas for watershed points
        DOP843 is better

    """
    # TODO assert about refine statement.
    # TODO assert density cutofff is not greater than iso_val
    # Assert length of basins matchces the centers.
    # Using ODE solver find the actual maximas
    maximas = centers
    if optimize_centers:
        print(maximas)
        maximas = np.array(
            [gradient_path(x, grad_func, t_span=(0, 10), method="Radau",
                           first_step=1e-9, max_step=1e-2) for x in centers],
            dtype=np.float64
        )
        print("New maximas: \n ", maximas)
        # Check duplicates
        distance = cdist(maximas, maximas)
        distance[np.diag_indices(len(maximas))] = 1.0  # Set diagonal elements to one
        if np.any(distance < 1e-8):
            raise RuntimeError(f"Optimized maximas contains duplicates: \n {maximas}.")

    numb_maximas = len(maximas)
    angular_pts = AngularGrid(degree=angular).points if isinstance(angular, int) else angular
    r, thetas, phis = convert_cart_to_sph(angular_pts).T
    numb_ang_pts = len(thetas)

    # TODO: have all zeros is inefficient due to refinement
    r_func = [np.zeros((numb_ang_pts,), dtype=np.float64) for _ in range(numb_maximas)]
    oas = [[] for _ in range(numb_maximas)]  # outer atomic surface
    ias = [[] for _ in range(numb_maximas)]  # inner atomic surface.
    basin_ias = [[] for _ in range(numb_maximas)]  # basin ids for inner atomic surface.
    refined_ang = [] if refine else None
    maxima_to_do = range(0, numb_maximas) if type(refine) == type(True) else [refine]  # for refinement
    for i_maxima, maxima in enumerate(maximas):
        # Maximas aren't usually large, so doing this is okay. Quick fix to use refinement without
        #  re-writing this function into seperate functions.
        if i_maxima in maxima_to_do:
            print("Start: Maxima ", maxima)
            other_maximas = np.delete(maximas, i_maxima, axis=0)
            other_beta_sph = -np.inf  # Infinity so that the if-statement holds true
            beta_sph_max = -np.inf
            if beta_sphere is not None:
                beta_sph_max = beta_sphere[i_maxima]
                other_beta_sph = [beta_sphere[i] for i in range(0, numb_maximas) if i != i_maxima]

            for i_ang in range(0, numb_ang_pts):  # Go through each point of the sphere
                print("I_ang ", i_ang)
                cart_sphere_pt, theta, phi = angular_pts[i_ang], thetas[i_ang], phis[i_ang]

                # Do backtracing on the ray
                radial = rgrids.points
                radial = radial if beta_sphere is None else radial[radial > beta_sphere[i_maxima]]

                bnd_pt, is_ray_to_inf, i_iso, found_watershed_on_ray, basin_id = solve_for_basin_bnd_pt(
                    dens_cutoff,maxima, radial, cart_sphere_pt, density_func, grad_func,
                    other_maximas, bnd_err, ss_watershed, iso_val, beta_sph_max, other_beta_sph
                )
                # If the ray tends towards infinity instead, solve for the isosurface value.
                if is_ray_to_inf:
                    bnd_pt = solve_for_isosurface_pt(
                        i_iso, radial, maxima, cart_sphere_pt, density_func, iso_val,
                        iso_err
                    )

                r_func[i_maxima][i_ang] = np.linalg.norm(bnd_pt - maxima)
                if is_ray_to_inf:
                    oas[i_maxima].append(i_ang)
                elif found_watershed_on_ray:
                    ias[i_maxima].append(i_ang)
                    basin_ias[i_maxima].append(basin_id)

                print("")

            if type(refine) == type(True) and refine:  # refine can be integer, so this ignores it.
                # Take convex hull between ias and oas and construct additional points in that region.
                #  `new_pts` is concatenated to angular grids and is in cartesian coordinates.
                print("IAS ", ias[i_maxima])
                print("OAS", oas[i_maxima])
                new_pts = get_closest_points_between_ias_and_oas(
                    ias[i_maxima], oas[i_maxima], angular_pts, r_func[i_maxima], maxima
                )
                print("new pts ", new_pts, np.linalg.norm(new_pts, axis=1))
                # Re-do this qtaim algortihm only on this center
                # TODO Fix this with other maximas.
                refined_qtaim = qtaim_line_search(rgrids, new_pts, maximas, density_func,
                                                  grad_func, dens_cutoff, iso_val, ss_watershed,
                                                  bnd_err, iso_err, beta_sphere=beta_sphere,
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
    result = SurfaceQTAIM(r_func, [rgrids], angular, maximas, oas, ias, basin_ias, refined_ang)
    # if type(refine) == type(True) and refine:
    #     import matplotlib
    #     import matplotlib.pyplot as plt
    #     from mpl_toolkits import mplot3d
    #     matplotlib.use("Qt5Agg")
    #     fig = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     i = 0
    #     p = centers
    #     ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="g", s=60)
    #     p = result.get_ias_of_basin(i)
    #     ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="k")
    #     p = result.get_oas_of_basin(i)
    #     ax.scatter(p[:, 0], p[:, 1], p[:, 2], color="r")
    #     plt.show()
    return SurfaceQTAIM(r_func, [rgrids], angular, maximas, oas, ias, basin_ias, refined_ang)
