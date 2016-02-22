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
'''The Utility Module.'''


import numpy as np


__all__ = ['doc_inherit', 'CubeGen']


def doc_inherit(base_class):
    """Docstring inheriting method descriptor

       doc_inherit decorator

       Usage:

       .. code-block:: python

            class Foo(object):
                def foo(self):
                    "Frobber"
                    pass

            class Bar(Foo):
                @doc_inherit(Foo)
                def foo(self):
                    pass

       Now, ``Bar.foo.__doc__ == Bar().foo.__doc__ == Foo.foo.__doc__ ==
       "Frobber"``
    """
    def decorator(method):
        overridden = getattr(base_class, method.__name__, None)
        if overridden is None:
            raise NameError('Can\'t find method \'%s\' in base class.')
        method.__doc__ = overridden.__doc__
        return method
    return decorator



class CubeGen(object):
    '''
    Class for generating a cubic grid and writing cube files.
    '''
    def __init__(self, numbers, pseudo_numbers, coordinates, origin, axes, shape):
        '''
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
        '''
        self._numbers = numbers
        self._pseudo_numbers = pseudo_numbers
        self._coordinates = coordinates
        self._origin = origin
        self._axes = axes
        self._shape = shape
        #
        # Make cubic grid
        #
        # Number of points along x, y and z axis
        npoints_x, npoints_y, npoints_z = self._shape
        # Total number of grid points
        self._npoints = npoints_x * npoints_y * npoints_z
        # Make an array to store coordinates of grid points
        self._gridpoints = np.zeros((self._npoints, 3))
        # Compute coordinates of grid points relative to the origin
        self._gridpoints += self._origin
        count = 0
        for nx in range(npoints_x):
            for ny in range(npoints_y):
                for nz in range(npoints_z):
                    coordinate = np.dot(np.array([nx, ny, nz]), self._axes)
                    self._gridpoints[count, :] += coordinate
                    count += 1

    @classmethod
    def from_molecule(cls, numbers, pseudo_numbers, coordinates, spacing=0.2, threshold=5.0):
        '''
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
        '''
        # maximum and minimum value of x, y and z coordinates
        max_coordinate = np.amax(coordinates, axis=0)
        min_coordinate = np.amin(coordinates, axis=0)
        # Compute the required number of points along x, y, and z axis
        shape = (max_coordinate - min_coordinate + 2.0 * threshold) / spacing
        shape = np.ceil(shape)
        shape = np.array(shape, int)
        # Compute coordinates of the cubic grid origin
        middle = (max_coordinate + min_coordinate) / 2.0
        origin = middle - 0.5 * shape * spacing
        # Compute the unit vectors of the cubic grid's coordinate system
        axes = np.diag([spacing, spacing, spacing])

        return cls(numbers, pseudo_numbers, coordinates, origin, axes, shape)

    @classmethod
    def from_cube(cls, filename):
        '''
        Initialize ``CubeGen`` class based on the grid specifications of a cube file.

        Parameters
        ----------
        filename : str
            Cube file name with *.cube extension.
        '''
        if not filename.endswith('.cube'):
            raise ValueError('Arguemnt filename should be a cube file with *.cube extension!')

        # Extract the specifications of the cubic grid from cube file's header
        numbers, pseudo_numbers, coordinates, origin, axes, shape = cls._read_cube_header(filename)

        return cls(numbers, pseudo_numbers, coordinates, origin, axes, shape)

    @property
    def numbers(self):
        '''
        Atomic number of the atoms in the molecule.
        '''
        return self._numbers

    @property
    def pseudo_numbers(self):
        '''
        Pseudo-number of the atoms in the molecule.
        '''
        return self._pseudo_numbers

    @property
    def coordinates(self):
        '''
        Cartesian coordinates of the atoms in the molecule.
        '''
        return self._coordinates

    @property
    def origin(self):
        '''
        Cartesian coordinate of the cubic grid origin.
        '''
        return self._origin

    @property
    def axes(self):
        '''
        The three vectors, stored as rows of axes array, defining the Cartesian
        coordinate system used to build the cubic grid.
        '''
        return self._axes

    @property
    def shape(self):
        '''
        Number of grid points along `x`, `y`, and `z` axis.
        '''
        return self._shape

    @property
    def npoints(self):
        '''
        Total number of grid points.
        '''
        return self._npoints

    @property
    def gridpoints(self):
        '''
        Cartesian coordinates of the cubic grid points.
        '''
        return self._gridpoints


    def dump_cube(self, filename, data):
        '''
        Write the data evaluated on grid points into a *.cube file.

        Parameters
        ----------
        filename : str
            Cube file name with *.cube extension.
        data : np.ndarray, shape=(npoints,)
            An array containing the evaluated scalar property on the grid points.
        '''
        if not filename.endswith('.cube'):
            raise ValueError('Arguemnt filename should be a cube file with `*.cube` extension!')
        if data.size != self._npoints:
            raise ValueError('Argument data should have the same size as the grid. {0}!={1}'.format(data.size, self._npoints))

        # Write data into the cube file
        with open(filename, 'w') as f:
            title = 'Cubefile created with HORTON CHEMTOOLS'
            # writing the cube header:
            print >> f, title
            print >> f, 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
            natom = len(self._numbers)
            x, y, z = self._origin
            print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (natom, x, y, z)
            rvecs = self._axes
            for i in xrange(3):
                x, y, z = rvecs[i]
                print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (self._shape[i], x, y, z)
            for i in xrange(natom):
                q = self._pseudo_numbers[i]
                x, y, z = self._coordinates[i]
                print >> f, '%5i % 11.6f % 11.6f % 11.6f % 11.6f' % (self._numbers[i], q, x, y, z)
            # writing the cube data:
            counter = 0
            for value in data.flat:
                f.write(' % 12.5E' % value)
                if counter % 6 == 5:
                    f.write('\n')
                counter += 1

    @staticmethod
    def _read_cube_header(filename):
        '''
        Return specifications of the cubic grid from the given cube file.

        Parameters
        ----------
        filename : str
            Cube file name with *.cube extension.
        '''
        with open(filename) as f:
            # skip the title
            f.readline()
            # skip the second line
            f.readline()

            def read_grid_line(line):
                '''Read a number and (x, y, z) coordinate from the cube file line.'''
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
                '''Read atomic number and (x, y, z) coordinate from the cube file line.'''
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
