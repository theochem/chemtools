# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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
"""Grid Wrapper Module."""


from horton import BeckeMolGrid
from chemtools.wrappers.molecule import Molecule


__all__ = ['BeckeGrid']


class BeckeGrid(object):
    """Grid class for wrapping grid module from HORTON package.

    Usage
    -----
    Initialization:
    >>> grid_model = BeckeGrid.from_molecule(mol, 'fine')
    >>> my_grid = grid_model.grid
    Or
    >>> grid_model = BeckeGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers, 'fine')
    >>> my_grid = grid_model.grid

    Change grid:
    >>> grid_model.grid_type = 'coarse'
    >>> new_grid_1 = grid_model.grid
    Or
    >>> grid_model.grid_type = 'linear:1e-5:2e1:120:110'
    >>> new_grid_2 = grid_model.grid

    Attributes
    ----------
    grid_type : str
        Information about current grid properties.
    rname : str
        Type of customized radial grid.
    rpoint : int
        The number of points for the angular Lebedev-Laikov grid.
    rrad : int
        The number of radial grid points.
    rrange : tuple(float, float)
        The first and the last radial grid point in angstroms.
    grid : BeckeMolGrid, read-only
        Generated Becke integration grid for given parameters.

    """

    def __init__(self, coordinates, numbers, pseudo_numbers, specification='medium', k=3, random_rotate=False,
                 mode='discard'):
        """Initialize Grid object.

        Parameters
        ----------
        coordinates : np.ndarray(N, 3)
            Coordinates of given molecule.
        numbers : np.ndarray(N,)
            Atomic numbers of given molecule.
        pseudo_numbers : np.ndarray(N,)
            Pseudo-potential core charges.
        specification : str, optional
            A specification of Grid property.
        k : int, optional
            The order of the switching function in Becke's weighting scheme.
        random_rotate : bool, optional
            Flag to control random rotation of spherical grids.
        mode : str, optional
            Select one of the following options regarding atomic subgrids.
            Avail choices: ['discard', 'keep', 'only'].

        """
        self._coordinates = coordinates
        self._numbers = numbers
        self._pseudo_n = pseudo_numbers
        self._grid_type = None  # used to store grid info for default type
        self._k = k
        self._random_rotate = random_rotate
        if mode in ['discard', 'keep', 'only']:
            self._mode = mode
        else:
            raise ValueError('Argument mode={0} is not valid!'.format(mode))
        # grid type specification
        self._custom_type = [None] * 5
        self.grid_type = specification

    @property
    def grid_type(self):
        """Type of current grid.

        Returns
        -------
        str
            Specific preset types or a string of properties.
        """
        if self._grid_type:
            return self._grid_type
        else:
            return ':'.join(map(str, self._custom_type))

    @grid_type.setter
    def grid_type(self, value):
        """Set the Property of Grid.

        Parameters
        ----------
        value : str
            One of preset types from ('coarse', 'medium', 'fine', 'veryfine',
            'ultrafine', 'insane') or a string of five properties split by ':'

        Raises
        ------
        ValueError
            Invalid input type name
        """
        valid_types = ('coarse', 'medium', 'fine', 'veryfine', 'ultrafine', 'insane')
        if value in valid_types:
            # set custom to None
            self._custom_type = [None] * 5
            # set new grid type
            self._grid_type = value
        else:
            # split input str
            ind_set = value.split(':')
            if len(ind_set) != 5:
                raise ValueError('Input grid_type: {} is not valid'.format(value))
            ind_set[1:3] = map(float, ind_set[1:3])
            ind_set[3:] = map(int, ind_set[3:])
            self.rname = ind_set[0]
            self.rrange = ind_set[1:3]
            self.rrad = ind_set[3]
            self.rpoint = ind_set[4]
            # map str to float and int

    @property
    def rname(self):
        """Type of the radial grid.

        Returns
        -------
        str
            the name of the grid type.
        """
        return self._custom_type[0]

    @rname.setter
    def rname(self, value):
        """Set the type of the radial grid from the preset types.

        Parameters
        ----------
        value : str
            The value of radial grid name,
            valid choices are ['linear', 'exp', 'power'].

        Raises
        ------
        ValueError
            If the given value is not one of the above
        """
        valid_inputs = ('linear', 'exp', 'power')
        if value not in valid_inputs:
            raise ValueError('Input value: {} is not valid'.format(value))
        self._custom_type[0] = value
        self._reset_grid_type()

    @property
    def rrange(self):
        """Specify the first and the last radial grid point in angstroms.

        Returns
        -------
        tuple(float, float)
            the start and stop grid point
        """
        return self._custom_type[1:3]

    @rrange.setter
    def rrange(self, value):
        """Specify the first and the last radial grid point in angstroms.

        Parameters
        ----------
        value : TYPE
            Description
        value : tuple(float, float)

        Raises
        ------
        TypeError
            A tuple of two positive numbers in the increasing sequence.
        """
        start, stop = value[:]
        if isinstance(start, (int, float)) and isinstance(stop, (int, float)):
            self._custom_type[1:3] = start, stop
            self._reset_grid_type()
        else:
            raise TypeError("Given values are not int")

    @property
    def rrad(self):
        """Return the number of radial grid points.

        Returns
        -------
        int
            Number of radial grid points
        """
        return self._custom_type[3]

    @rrad.setter
    def rrad(self, value):
        """Return the number of radial grid points.

        Parameters
        ----------
        value : int
            Number of radial points

        Raises
        ------
        TypeError
            The input need to be an integer
        """
        if not isinstance(value, int):
            raise TypeError('Given value is not an int')
        self._custom_type[3] = value
        self._reset_grid_type()

    @property
    def rpoint(self):
        """Return the number of points for the angular Lebedev-Laikov grid.

        Returns
        -------
        int
            The number of points on the angular grid
        """
        return self._custom_type[4]

    @rpoint.setter
    def rpoint(self, value):
        """Return the number of points for the angular Lebedev-Laikov grid.

        Parameters
        ----------
        value : int
            The value need to be in the following list:
            (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266,
             302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030,
             2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)

        Raises
        ------
        TypeError
            The value needs to be an integer
        ValueError
            The value needs to be one from the above tuple.
        """
        if not isinstance(value, int):
            raise TypeError('Given value is not an int')
        valid_ill = (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266,
                     302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030,
                     2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)
        if value not in valid_ill:
            raise ValueError("""Given value is not a valid value,
                check 'http://theochem.github.io/horton/2.1.0/te\
                ch_ref_grids.html#ref-grids' for reference""")
        self._custom_type[4] = value
        self._reset_grid_type()

    @property
    def grid(self):
        """Return a grid object based on set property.

        Returns
        -------
        BeckeMolGrid
            The BeckeMolGrid object generated by all the given properties.
        """
        if None in self._custom_type and self._grid_type is None:
            raise ValueError("Don't have enough info to generate a grid")
        return BeckeMolGrid(
            self._coordinates,
            self._numbers,
            self._pseudo_n,
            agspec=self.grid_type,
            k=self._k,
            random_rotate=self._random_rotate,
            mode=self._mode)

    @classmethod
    def from_molecule(cls,
                      mol,
                      grid_type='medium',
                      k=3,
                      random_rotate=False,
                      mode='discard'):
        """Construct a grid for given molecule.

        Parameters
        ----------
        mol : instance of Molecule
            Instance of Molecule class.
        grid_type : str, optional
            Information about current grid properties.
        k : int, optional
            The order of the switching function in Becke's weighting scheme.
        random_rotate : bool, optional
            Flag to control random rotation of spherical grids.
        mode : str, optional
            Select one of the following options regarding atomic subgrids.

        Returns
        -------
        BeckeGrid
            The Grid object for constructing integral grid

        Raises
        ------
        ValueError
            The given molecule is not a valid instance
        """
        if not isinstance(mol, Molecule):
            raise TypeError('Given molecule object is not valid')
        return cls(
            mol.coordinates,
            mol.numbers,
            mol.pseudo_numbers,
            grid_type=grid_type,
            k=k,
            random_rotate=random_rotate,
            mode=mode)

    def _reset_grid_type(self):
        """Set self._grid_type back to None."""
        self._grid_type = None
