# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2018 The ChemTools Development Team
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
"""Wrapper Module."""

import numpy as np
from horton import BeckeMolGrid

__all__ = ["Grid"]


class Grid(object):
    def __init__(self,
                 coordinates,
                 numbers,
                 pseudo_numbers,
                 grid_type='medium',
                 k=3,
                 random_rotate=True,
                 mode='discard'):
        self._coordinates = coordinates
        self._numbers = numbers
        self._pseudo_n = pseudo_numbers
        self._grid_type = None
        self.grid_type = grid_type
        self._k = k
        self._random_rotate = random_rotate
        self._mode = mode
        # grid type specification
        self._custom_type = [None] * 5

    @property
    def grid_type(self):
        if self._grid_type:
            return self._grid_type
        else:
            return ':'.join(map(str, self._custom_type))

    @grid_type.setter
    def grid_type(self, value):
        PRE_SET = ('coarse', 'medium', 'fine', 'veryfine', 'ultrafine',
                   'insane')
        if value in PRE_SET:
            # set custom to None
            self._custom_type = [None] * 5
            # set new grid type
            self._grid_type = value
        else:
            # split input str
            ind_set = value.split(':')
            if len(ind_set) != 5:
                raise ValueError('Input type: {} is not valid'.format(value))
            ind_set[1:3] = map(float, ind_set[1:3])
            ind_set[3:] = map(int, ind_set[3:])
            self.rname = ind_set[0]
            self.rrange = ind_set[1:3]
            self.rrad = ind_set[3]
            self.rpoint = ind_set[4]
            # map str to float and int

    @property
    def rname(self):
        return self._custom_type[0]

    @rname.setter
    def rname(self, value):
        valid_inputs = ('linear', 'exp', 'power')
        if value not in valid_inputs:
            raise ValueError('Input value: {} is not valid'.format(value))
        self._custom_type[0] = value
        self._reset_grid_type()

    @property
    def rrange(self):
        return self._custom_type[1:3]

    @rrange.setter
    def rrange(self, value):
        start, stop = value[:]
        if isinstance(start, (int, float)) and isinstance(stop, (int, float)):
            self._custom_type[1:3] = start, stop
            self._reset_grid_type()
        else:
            raise TypeError("Given values are not int")

    @property
    def rrad(self):
        return self._custom_type[3]

    @rrad.setter
    def rrad(self, value):
        if not isinstance(value, int):
            raise TypeError('Given value is not an int')
        self._custom_type[3] = value
        self._reset_grid_type()

    @property
    def rpoint(self):
        return self._custom_type[4]

    @rpoint.setter
    def rpoint(self, value):
        if not isinstance(value, int):
            raise TypeError('Given value is not an int')
        valid_ill = (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266,
                     302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030,
                     2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)
        if value not in valid_ill:
            raise ValueError("""Given value is not a legit input,
                check 'http://theochem.github.io/horton/2.1.0/te\
                ch_ref_grids.html#ref-grids' for reference""")
        self._custom_type[4] = value
        self._reset_grid_type()

    @property
    def grid(self):
        """return a grid object based on set property"""
        return BeckeMolGrid(
            self._coordinates,
            self._numbers,
            self._pseudo_n,
            agspec=self.grid_type,
            k=self._k,
            random_rotate=self._random_rotate,
            mode=self._mode)

    def _reset_grid_type(self):
        self._grid_type = None
