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
# pragma pylint: disable=invalid-name
"""Kinetic Energy Density Module."""


import numpy as np

from chemtools.utils.utils import doc_inherit
from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensGradTool, DensGradLapTool, DensGradLapKedTool


__all__ = ["KED"]


class KED(object):
    """Kinetic Energy Density Class."""

    def __init__(self, dens, grad, lap=None, kin=None):
        """Initialize class.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of grid points.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of grid points.
        lap : np.ndarray
            Laplacian of electron density evaluated on a set of grid points.
        kin : np.ndarray
            Positive-definite kinetic energy density evaluated on a set of grid points.

        """
        # initialize dens-based tools class
        if lap is None and kin is None:
            self._denstools = DensGradTool(dens, grad)
        elif kin is None:
            self._denstools = DensGradLapTool(dens, grad, lap)
        elif lap is None:
            self._denstools = DensGradTool(dens, grad)
            self._ked_pd = kin
        else:
            self._denstools = DensGradLapKedTool(dens, grad, lap, kin)

    @classmethod
    def from_molecule(cls, molecule, points, spin="ab", index=None):
        r"""Initialize class using instance of `Molecule` and grid points.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class.
        points : np.ndarray
            Cartesian coordinates, a 2D array with 3 columns, to calculate local properties.
        spin : str
            The type of occupied spin orbitals.
        index : sequence
            Sequence of integers representing the index of spin orbitals.

        """
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D array with 3 columns.")
        # compute density, gradient, & kinetic energy density on grid
        return cls(*molecule.compute_megga(points, spin=spin, index=index))

    @classmethod
    def from_file(cls, filename, points, spin="ab", index=None):
        """Initialize class from file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        points : np.ndarray
            Cartesian coordinates, a 2D array with 3 columns, to calculate local properties.
        spin : str
            The type of occupied spin orbitals.
        index : sequence
            Sequence of integers representing the index of spin orbitals.
        """
        molecule = Molecule.from_file(filename)
        return cls.from_molecule(molecule, points, spin, index)

    @property
    def points(self):
        """Coordinates of grid points."""
        return self._points

    @property
    @doc_inherit(DensGradTool, 'density')
    def density(self):
        return self._denstools.density

    @property
    @doc_inherit(DensGradLapKedTool, 'ked_positive_definite')
    def ked_positive_definite(self):
        if hasattr(self._denstools, 'ked_positive_definite'):
            return self._denstools.ked_positive_definite
        elif hasattr(self, '_ked_pd'):
            return self._ked_pd
        else:
            raise ValueError('Argument kin should be given when initializing the class.')

    @property
    @doc_inherit(DensGradTool, 'ked_thomas_fermi')
    def ked_thomas_fermi(self):
        return self._denstools.ked_thomas_fermi

    @property
    @doc_inherit(DensGradTool, 'ked_weizsacker')
    def ked_weizsacker(self):
        return self._denstools.ked_weizsacker
