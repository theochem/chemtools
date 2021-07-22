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
from grid.becke import BeckeWeights
from grid.molgrid import MolGrid
from grid.onedgrid import HortonLinear
from chemtools.wrappers.molecule import Molecule

__all__ = ['MolecularGrid']

#TODO: How to implement compute spherical average functions?
class MolecularGrid:
    """Becke-Lebedev molecular grid for numerical integrations."""

    def __init__(self, coordinates, numbers, pseudo_numbers,
                 specs='medium', k=3, rotate=False):

        self._coordinates = coordinates
        self._numbers = numbers
        self._pseudo_numbers = pseudo_numbers
        self._k = k
        self._rotate = rotate
        self._specs = specs

        # TODO: making points of angular an arg in __init__ generates an error
        onedg = HortonLinear(100)
        becke = BeckeWeights(order=self._k)
        points_of_angular = 110

        self._grid = MolGrid.horton_molgrid(self.coordinates, self.numbers,
                                            onedg, points_of_angular, becke)

    @classmethod
    def from_molecule(cls, molecule, specs='medium', k=3, rotate=False):
        if not isinstance(molecule, Molecule):
            raise TypeError('Argument molecule should be an instance of Molecule class.')
        coords, nums, pnums = molecule.coordinates, molecule.numbers, molecule.pseudo_numbers
        return cls(coords, nums, pnums, specs, k, rotate)

    @classmethod
    def from_file(cls, fname, *, specs='medium', k=3, rotate=False):
        mol = Molecule.from_file(fname)
        return cls.from_molecule(mol, specs, k, rotate)

    def __getattr__(self, item):
        return getattr(self._grid, item)

    @property
    def center(self):
        """Cartesian coordinates of atomic centers."""
        return self._coordinates

    @property
    def coordinates(self):
        """Cartesian coordinates of atomic centers."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic number of atomic centers."""
        return self._numbers

    @property
    def pseudo_numbers(self):
        """Pseudo atomic number of atomic centers."""
        return self._pseudo_numbers

    @property
    def npoints(self):
        """Number of grid points."""
        return self._grid.points.shape[0]

    @property
    def points(self):
        """Cartesian coordinates of grid points."""
        return self._grid.points

    @property
    def weights(self):
        """Integration weight of grid points."""
        return self._grid.weights

    def integrate(self, *value_arrays):
        r""" Integrate over the whole grid for given multiple value arrays.

        Product of all value_arrays will be computed element-wise then
        integrated on the grid with its weights.
        .. math::
            Integral = \int w(x) \prod_i f_i(x) dx

        Parameters
        ----------
        *value_arrays : np.ndarray(N, )
            One or multiple value array to integrate.
        """
        return self._grid.integrate(*value_arrays)

    def compute_spherical_average(self, value):
        raise NotImplementedError("Not yet implemented in grid")
