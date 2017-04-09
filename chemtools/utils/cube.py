# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
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
"""The Cube Module."""

import numpy as np
from horton import log, IOData

__all__ = ['CubeGen']


class CubeGen(object):
    """Class for generating a cubic grid and writing cube files."""

    def __init__(self, numbers, pseudo_numbers, coordinates, origin, axes, shape):
        """Initialize ``CubeGen`` class based on the origin, axes and shape of the cube.

        Parameters
        ----------
        numbers : np.ndarray, shape=(M,)
            Atomic number of `M` atoms in the molecule.
        pseudo_numbers : np.ndarray, shape=(M,)
            Pseudo-number of `M` atoms in the molecule.
        coordinates : np.ndarray, shape=(M, 3)
            Cartesian coordinates of `M` atoms in the molecule.
        origin : np.ndarray, shape=(3,)
            Cartesian coordinates of the cubic grid origin.
        axes : np.ndarray, shape=(3, 3)
            The three vectors, stored as rows of axes array,
            defining the Cartesian coordinate system used to build the
            cubic grid.
        shape : np.ndarray, shape=(3,)
            Number of grid points along `x`, `y`, and `z` axis.
        """
        self._numbers = numbers
        self._pseudo_numbers = pseudo_numbers
        self._coordinates = coordinates
        if origin.shape[0] != 3:
            raise ValueError('Argument origin should be an np.ndarray with shape=(3,)')
        self._origin = origin
        if axes.shape[0] != 3 or axes.shape[1] != 3:
            raise ValueError('Argument axes should be an np.ndarray with shape=(3, 3)')
        self._axes = axes
        if shape.shape[0] != 3:
            raise ValueError('Argument shape should be an np.ndarray with shape=(3,)')
        self._shape = shape
        #
        # Make cubic grid
        #
        # Number of points along x, y and z axis
        npoints_x, npoints_y, npoints_z = self._shape
        # Total number of grid points
        self._npoints = npoints_x * npoints_y * npoints_z
        # Make an array to store coordinates of grid points
        self._points = np.zeros((self._npoints, 3))
        # Compute coordinates of grid points relative to the origin
        self._points += self._origin
        count = 0
        for nx in range(npoints_x):
            for ny in range(npoints_y):
                for nz in range(npoints_z):
                    coordinate = np.dot(np.array([nx, ny, nz]), self._axes)
                    self._points[count, :] += coordinate
                    count += 1

        # log information
        self._log_init()

    @classmethod
    def from_molecule(cls, numbers, pseudo_numbers, coordinates, spacing=0.2,
                      threshold=5.0, rotate=True):
        """
        Initialize ``CubeGen`` class based on the Cartesian coordinates of the molecule.

        Parameters
        ----------
        numbers : np.ndarray, shape=(M,)
            Atomic number of `M` atoms in the molecule.
        pseudo_numbers : np.ndarray, shape=(M,)
            Pseudo-number of `M` atoms in the molecule.
        coordinates : np.ndarray, shape=(M, 3)
            Cartesian coordinates of `M` atoms in the molecule.
        spacing : float, default=0.2
            Increment between grid points along `x`, `y` and `z` direction.
        threshold : float, default=5.0
            The extension of the cube on each side of the molecule.
        rotate : bool, default=True
            When True, the molecule is rotated so the axes of the cube file are
            aligned with the principle axes of rotation of the molecule.
        """
        # calculate center of mass of the nuclear charges:
        totz = np.sum(pseudo_numbers)
        com = np.dot(pseudo_numbers, coordinates) / totz

        if rotate:
            # calculate moment of inertia tensor:
            itensor = np.zeros([3, 3])
            for i in range(pseudo_numbers.shape[0]):
                xyz = coordinates[i] - com
                r = np.linalg.norm(xyz)**2.0
                tempitens = np.diag([r, r, r])
                tempitens -= np.outer(xyz.T, xyz)
                itensor += pseudo_numbers[i] * tempitens

            _, v = np.linalg.eigh(itensor)
            new_coordinates = np.dot((coordinates - com), v)
            axes = spacing * v

        else:
            # Just use the original coordinates
            new_coordinates = coordinates
            # Compute the unit vectors of the cubic grid's coordinate system
            axes = np.diag([spacing, spacing, spacing])

        # maximum and minimum value of x, y and z coordinates
        max_coordinate = np.amax(new_coordinates, axis=0)
        min_coordinate = np.amin(new_coordinates, axis=0)
        # Compute the required number of points along x, y, and z axis
        shape = (max_coordinate - min_coordinate + 2.0 * threshold) / spacing
        shape = np.ceil(shape)
        shape = np.array(shape, int)
        # Compute origin
        origin = com - np.dot((0.5 * shape), axes)

        return cls(numbers, pseudo_numbers, coordinates, origin, axes, shape)

    @classmethod
    def from_cube(cls, filename):
        r"""
        Initialize ``CubeGen`` class based on the grid specifications of a cube file.

        Parameters
        ----------
        filename : str
            Cube file name with \*.cube extension.
        """
        if not filename.endswith('.cube'):
            raise ValueError('Argument filename should be a cube file with *.cube extension!')

        # Extract the specifications of the cubic grid from cube file's header
        numbers, pseudo_numbers, coordinates, origin, axes, shape = cls._read_cube_header(filename)

        return cls(numbers, pseudo_numbers, coordinates, origin, axes, shape)

    @classmethod
    def from_file(cls, filename, spacing=0.2, threshold=5.0, rotate=True):
        """
        Initialize ``CubeGen`` class based on the grid specifications of a file.

        Parameters
        ----------
        filename : str
            file name with, readable with HORTON's IOData.
        spacing : float, default=0.2
            Increment between grid points along `x`, `y` and `z` direction.
        threshold : float, default=5.0
            The extension of the cube on each side of the molecule.
        rotate : bool, default=True
            When True, the molecule is rotated so the axes of the cube file are
            aligned with the principle axes of rotation of the molecule.
        """
        # Load file
        mol = IOData.from_file(filename)
        return cls.from_molecule(mol.numbers, mol.pseudo_numbers, mol.coordinates, spacing,
                                 threshold, rotate)

    @property
    def numbers(self):
        """Atomic number of the atoms in the molecule."""
        return self._numbers

    @property
    def pseudo_numbers(self):
        """Pseudo-number of the atoms in the molecule."""
        return self._pseudo_numbers

    @property
    def coordinates(self):
        """Cartesian coordinates of the atoms in the molecule."""
        return self._coordinates

    @property
    def origin(self):
        """Cartesian coordinate of the cubic grid origin."""
        return self._origin

    @property
    def axes(self):
        """Array with axes of the cube.

        The three vectors, stored as rows of axes array, defining the Cartesian
        coordinate system used to build the cubic grid.
        """
        return self._axes

    @property
    def shape(self):
        """Number of grid points along `x`, `y`, and `z` axis."""
        return self._shape

    @property
    def npoints(self):
        """Total number of grid points."""
        return self._npoints

    @property
    def points(self):
        """Cartesian coordinates of the cubic grid points."""
        return self._points

    def _log_init(self):
        """Log an overview of the cube's properties."""
        if log.do_medium:
            log('Initialized cube: %s' % self.__class__)
            log.deflist([('Origin ', self._origin),
                         ('Axes[0]', self._axes[0]),
                         ('Axes[1]', self._axes[1]),
                         ('Axes[2]', self._axes[2]),
                         ('Shape  ', self._shape)])
            log.blank()

    def dump_cube(self, filename, data):
        r"""Write the data evaluated on grid points into a cube file.

        Parameters
        ----------
        filename : str
            Cube file name with \*.cube extension.
        data : np.ndarray, shape=(npoints,)
            An array containing the evaluated scalar property on the grid points.
        """
        if not filename.endswith('.cube'):
            raise ValueError('Argument filename should be a cube file with `*.cube` extension!')
        if data.size != self._npoints:
            raise ValueError('Argument data should have the same size as the grid. ' +
                             '{0}!={1}'.format(data.size, self._npoints))

        # Write data into the cube file
        with open(filename, 'w') as f:
            # writing the cube header:
            f.write('Cubefile created with HORTON CHEMTOOLS\n')
            f.write('OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n')
            natom = len(self._numbers)
            x, y, z = self._origin
            f.write('{0:5d} {1:11.6f} {2:11.6f} {3:11.6f}\n'.format(natom, x, y, z))
            rvecs = self._axes
            for i, (x, y, z) in zip(self._shape, rvecs):
                f.write('{0:5d} {1:11.6f} {2:11.6f} {3:11.6f}\n'.format(i, x, y, z))
            for i, q, (x, y, z) in zip(self._numbers, self._pseudo_numbers, self._coordinates):
                f.write('{0:5d} {1:11.6f} {2:11.6f} {3:11.6f} {4:11.6f}\n'.format(i, q, x, y, z))
            # writing the cube data:
            num_chunks = 6
            for i in range(0, data.size, num_chunks):
                row_data = data.flat[i:i+num_chunks]
                f.write((row_data.size*' {:12.5E}').format(*row_data))
                f.write('\n')


    def weights(self, method='R'):
        """
        Return integration weights at every point on the cubic grid.

        Parameters
        ----------
        method : str, default='R0'
            The method for computing the integration weights at every point on the grid. Options:

                - 'R' method perfors rectangle/trapezoidal rule, without assuming that the function
                  is close to zero at the edges of the grid.
                - 'R0' method performing rectangle/trapezoidal rule, assuming that the function is
                  very close to zero at the edges of the grid.
        """
        if method == 'R':
            volume = np.linalg.norm(self._shape[0] * self._axes[0])
            volume *= np.linalg.norm(self._shape[1] * self._axes[1])
            volume *= np.linalg.norm(self._shape[2] * self._axes[2])
            numpnt = 1.0 * self._npoints
            weights = np.full(self._npoints, volume / numpnt)

        elif method == 'R0':
            volume = np.linalg.norm((self._shape[0] + 1.0) * self._axes[0])
            volume *= np.linalg.norm((self._shape[1] + 1.0) * self._axes[1])
            volume *= np.linalg.norm((self._shape[2] + 1.0) * self._axes[2])

            numpnt = (self._shape[0] + 1.0) * (self._shape[1] + 1.0) * (self._shape[2] + 1.0)
            weights = np.full(self._npoints, volume / numpnt)

        else:
            raise ValueError('Argument method {0} is not known.'.format(method))
        return weights

    def integrate(self, data, method='R0'):
        """
        Integrate the data on a cubic grid.

        Parameters
        ----------
        data : np.ndarray, shape=(npoints, m)
            Data at every point on the grid given as an array. The size of axis=0 of this array
            should equal the number of grid points.

        method : str, default='R0'
            The method for computing the integration weights at every point on the grid. Options:

                - 'R' method perfors rectangle/trapezoidal rule, without assuming that the function
                  is close to zero at the edges of the grid.
                - 'R0' method performing rectangle/trapezoidal rule, assuming that the function is
                  very close to zero at the edges of the grid.
        """
        if data.shape[0] != self._npoints:
            raise ValueError('Argument data should have the same size as the grid for axis=0. ' +
                             '{0}!={1}'.format(data.shape[0], self._npoints))
        value = np.tensordot(self.weights(method=method), data, axes=(0, 0))
        return value

    @staticmethod
    def _read_cube_header(filename):
        """
        Return specifications of the cubic grid from the given cube file.

        Parameters
        ----------
        filename : str
            Cube file name with *.cube extension.
        """
        with open(filename) as f:
            # skip the title
            f.readline()
            # skip the second line
            f.readline()

            def read_grid_line(line):
                """Read a number and (x, y, z) coordinate from the cube file line."""
                words = line.split()
                return (
                    int(words[0]),
                    np.array([float(words[1]), float(words[2]), float(words[3])], float)
                    # all coordinates in a cube file are in atomic units
                )

            # number of atoms and origin of the grid
            natom, origin = read_grid_line(f.readline())
            # numer of grid points in A direction and step vector A, and so on
            shape0, axis0 = read_grid_line(f.readline())
            shape1, axis1 = read_grid_line(f.readline())
            shape2, axis2 = read_grid_line(f.readline())
            shape = np.array([shape0, shape1, shape2], int)
            axes = np.array([axis0, axis1, axis2])

            def read_coordinate_line(line):
                """Read atomic number and (x, y, z) coordinate from the cube file line."""
                words = line.split()
                return (
                    int(words[0]), float(words[1]),
                    np.array([float(words[2]), float(words[3]), float(words[4])], float)
                    # all coordinates in a cube file are in atomic units
                )

            numbers = np.zeros(natom, int)
            pseudo_numbers = np.zeros(natom, float)
            coordinates = np.zeros((natom, 3), float)
            for i in xrange(natom):
                numbers[i], pseudo_numbers[i], coordinates[i] = read_coordinate_line(f.readline())
                # If the pseudo_number field is zero, we assume that no effective core
                # potentials were used.
                if pseudo_numbers[i] == 0.0:
                    pseudo_numbers[i] = numbers[i]

        return numbers, pseudo_numbers, coordinates, origin, axes, shape
