
r"""
Data structure that holds the outer-atomic (OAS) and intra-atomic (IAS) surfaces.

Can be used for
- analyzing the IAS and OAS.
- integration over basins.
"""
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

from grid.angular import AngularGrid
from grid.atomgrid import AtomGrid
from grid.basegrid import Grid

__all__ = ["SurfaceQTAIM"]


class SurfaceQTAIM():
    def __init__(self, r_func, angular_degs, maximas, oas, ias, basins_ias, iso_val, beta_spheres,
                 refined_ang=None):
        self._r_func = r_func
        self._maximas = maximas
        self._angular_degs = angular_degs
        self._oas = oas
        self._ias = ias
        self._basins_ias = basins_ias
        self._refined_ang = refined_ang
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
        # angular points.
        return self._r_func

    @property
    def oas(self):
        # List[List[int]] : First list is over basins, second over indices of points of outeratomic
        # surface.
        return self._oas

    @property
    def ias(self):
        # List[List[int]] : First list is over basins, second over indices of points of interatomic
        # surface.
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

    def save(self, filename):
        save_dict = {
            "ias" : np.array(self.ias, dtype=np.object),
            "oas" : np.array(self.oas, dtype=np.object),
            "basin_ias": np.array(self.basins_ias, dtype=np.object),
            "maximas": np.array(self.maximas),
            "angular_degs" : np.array(self.angular_degs),
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
        if self.refined_ang is not None:
            points = np.vstack((points, self.refined_ang[i_basin]))
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
        points = self.maximas[i_basin] + self.r_func[i_basin][:, None] * sph_pts
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
                if i != i_basin:
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

    def interpolate_radial_func(self, method="smooth", ias=False, oas=False):
        # if method not in ["smooth", ]
        if ias and oas:
           raise ValueError(f"Both {ias} and {oas} cannot be true.")
        if ias:
            #TODO
            pass
        raise NotImplementedError(f"Not implemented yet.")
