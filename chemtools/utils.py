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
from horton.io.iodata import IOData
import numpy as np

__all__ = ['doc_inherit', 'CubeGen']

'''The Utility Module.'''

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
    '''An easy-to-use method to generate cubic grids and cube files)
    '''
    def __init__(self, molecule, origin, axes, shape):
        '''
        Parameters
        ----------
        molecule :
        instance of ``horton.IOData``.
        origin :
        origin of the cube grid, coordinates of gridpoint[0].
        axes:
        vectors along which the cube is constructed.
        shape:
        number of points to sample along each direction.
        '''
        # add test to see if isinstance of horton.IOData?
        self._mol = molecule
        self._origin = origin
        self._axes = axes
        self._shape = shape

        self._npoints = np.prod(self._shape)
        self._gridpoints = np.full((self._npoints,3), self._origin)
        cnt = 0
        for i in range(0,self._shape[0]):
            for j in range(0,self._shape[1]):
                for k in range(0,self._shape[2]):
                    v = np.array([1.0*i,1.0*j,1.0*k])
                    self._gridpoints[cnt, :] +=  np.dot(v,self._axes)
                    cnt += 1

    @classmethod
    def from_molecule(cls, molecule, spacing=0.2, treshold=5.0):
        '''
        Parameters
        ----------
        molecule :
        instance of ``horton.IOData``.
        spacing :
        size of the step taken in each direction.
        treshold:
        extra treshhold in each direction to ensure the plotted function is within the cubefile.
        '''

        shape = (np.amax(molecule.coordinates, axis=0) - np.amin(molecule.coordinates, axis=0) + 2.0*treshold)/spacing

        shape = np.ceil(shape)

        origin = - 0.5*shape*spacing

        shape = np.array(shape, int)

        axes=np.zeros((3, 3),float)
        np.fill_diagonal(axes,spacing)

        return cls(molecule,origin,axes,shape)

    @classmethod
    def from_file(cls, filename):

        if filename.endswith('.cube'):
            molecule = IOData.from_file(filename)
            with open(filename) as f:
                origin, axes, shape = _read_cube_header(f)


        else:
            raise ValueError('Unknown file format for reading: %s' % filename)

        return cls(molecule, origin, axes, shape)

    @property
    def mol(self):
        '''
        instance of ``horton.IOData``.
        '''
        return self._mol

    @property
    def origin(self):
        '''
        origin of the cube grid, coordinates of gridpoint[0].
        '''
        return self._origin

    @property
    def axes(self):
        '''
        vectors along which the cube is constructed.
        '''
        return self._axes

    @property
    def shape(self):
        '''
        shape of the cubic grid.
        '''
        return self._shape

    @property
    def npoints(self):
        '''
        number of moints in the cubic grid.
        '''
        return self._npoints

    @property
    def gridpoints(self):
        '''
        xyz coordinates of the gridpoints. The shape of this numpy array is [npoints, 3].
        '''
        return self._gridpoints


    def dump_cube(self, filename, data):
        '''Write a IOData to a .cube file.

           **Arguments:**

           filename
                The name of the file to be written. This usually the extension
                ".cube".

           data
            data for the cubefile.
        '''
        with open(filename, 'w') as f:
            title = getattr(self._mol, 'title', 'Created with HORTON CHEMTOOLS')
            # writing the cube header:
            print >> f, title
            print >> f, 'OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z'
            natom = len(self._mol.numbers)
            x, y, z = self._origin
            print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (natom, x, y, z)
            rvecs = self._axes
            for i in xrange(3):
                x, y, z = rvecs[i]
                print >> f, '%5i % 11.6f % 11.6f % 11.6f' % (self._shape[i], x, y, z)
            for i in xrange(natom):
                q = self._mol.pseudo_numbers[i]
                x, y, z = self._mol.coordinates[i]
                print >> f, '%5i % 11.6f % 11.6f % 11.6f % 11.6f' % (self._mol.numbers[i], q, x, y, z)
            # writing the cube data:
            counter = 0
            for value in data.flat:
                f.write(' % 12.5E' % value)
                if counter%6 == 5:
                    f.write('\n')
                counter += 1


def _read_cube_header(f):
    # skip the title
    f.readline()
    # skip the second line
    f.readline()

    def read_grid_line(line):
        """Read a grid line from the cube file"""
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

    return origin, axes, shape
