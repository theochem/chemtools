import numpy as np
from scipy.sparse import lil_matrix
from scipy.spatial import ConvexHull, Voronoi
from scipy.spatial.distance import cdist
from grid.cubic import UniformGrid, _HyperRectangleGrid


__all__ = [""]

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
    ridge_index = (set(ith_voronoi_ridges) & set(i_nbh_voronoi_ridges))  # Take intersection
    assert len(ridge_index) == 1
    ridge_index = ridge_index.pop()
    # Get the voronoi vertices  via : voronoi.vertices[delaunay.ridge_vertices[r_{ij}]].
    #   voronoi.ridge_vertices[r_{ij}] makes sure it doesn't have -1 in it.
    ridge_vertices = voronoi.vertices[voronoi.ridge_vertices[ridge_index]]
    print("ridge vertices", ridge_vertices, "voronoi ridge vertices",
          voronoi.ridge_vertices[ridge_index])
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
        print("Nbh Index ", nbh_index, "Nbh Point ", voronoi.points[nbh_index],
              "Density Diff ", density_diff, "Densities ", density_vals[nbh_index],
              density_vals[index])
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
        if -1 in indices:  # some regions can be opened
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


def qtaim(grid_pts, density_vals, grad_func=None, num_centers=None,
          remove_duplicates=True, bounding_box=True):
    r"""
    Find basins from Yu-Trinkle algorithm on a cubic grid or a general grid using Voronoi diagrams.

    If a general grid is used, then problems may arise due to instabilities of constructing
    a Voronoi diagrams.  Providing a cubic grid is more robust.

    Parameters
    ----------
    grid_pts: theochem/grid._HyperRectangleGrid or np.ndarray
        If it is the latter, then it is a cubic grid and whose Voronoi diagrams is known.
        If it is an array of points then the Voronoi diagram is implemented.
    density_vals: ndarrray
        Density values of each of the grid points.
    grad_func: Callable(ndarray(N,3) -> ndarray(N,3))
        The gradient of the density values.  If provided, then it calculates the weights of
        the watershed points more accuaretely.
    num_centers: int, optional
        The number of centers/local maximas.
    remove_duplicates: bool
        If true, then it removes duplicates from `grid_pts`. This is due to construction
        of Voronoi diagrams.
    bounding_box: bool
        If true, then constructs a bounding box around the atom.

    Returns
    -------
    dict:
        Dictionary with the following keys:

        "basin_cont": CSC
            Sparse array (CSC format) with columns correspond to each basin.
        "maxima_indices": ndarray
            Array that holds which points correspond to the local maxima based on the grid.
        "watershed_indices": ndarray
            Array that holds indices of which points correspond to watershed points.
        "voronoi_volumes": ndarray or float
            Corresponds to the volume of the Voronoi diagram.

    Notes
    -----
    1) Sort based on density values
    2) Calculate the Voronoi Diagram
    3) Find and store all neighbors of each point using the Voronoi Diagram, may miss
        some points on boundary.
    4) Go through each point
        i) See if the neighbors of that point are assigned.
        ii) If no neighbours are assigned and it has neighbors then it is a maxima point
        iii) If no neighbors are assigned and no neighbor information is found, then assing
            its weight based on the average of its closest three points.
        iv) If all neighbors that were assigned, are assigned to a single basin, assign this
            to that basin.
        v) If some neighbors are assigned to different basins, then this point is a watershed point.
            and it solves using Yu-Trinkle method for the fractional points.

    """
    # Assert method values
    if not isinstance(grid_pts, (_HyperRectangleGrid, np.ndarray)):
        raise TypeError(f"Points should either be a numpy array or a UniformGrid object.")
    is_cubic_grid = isinstance(grid_pts, _HyperRectangleGrid)
    points = grid_pts if not is_cubic_grid else grid_pts.points

    if remove_duplicates:
        points, indices = np.unique(points, return_index=True, axis=0)
        density_vals = density_vals[indices]

    original_num_pts = points.shape[0]
    if not is_cubic_grid:
        if bounding_box:
            points = _points_on_bounding_box(points, extension=0.25, step_size=0.1)
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
    # TODO: When density values are sorted, maybe group them based on similar values and pick
    #  the ones closest to the points.
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
                weights = _assign_weight_average_neighbors(
                    index, basin_cont, voronoi, original_num_pts
                )
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
        print("")

    # TODO: Update watershed indices
    # Calculate Voronoi volumes
    volume = np.prod(np.array([np.linalg.norm(axis) for axis in grid_pts.axes])) if is_cubic_grid \
        else voronoi_volumes(voronoi)[:original_num_pts]
    return {"basin_cont": basin_cont.basin.tocsc(), "maxima_indices": maxima_indices,
            "watershed_indices": watershed_indices, "voronoi_volumes": volume}