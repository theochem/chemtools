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
# pragma pylint: disable=wrong-import-position
"""Module for Electron Localization Function (ELF) analysis of Quantum Chemistry Output Files."""


from chemtools.utils.cube import CubeGen
from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensityLocalTool
from chemtools.outputs import print_vmd_script_isosurface


__all__ = ['ELF']


class ELF(object):
    """Electron Localization Function (ELF) Class."""

    def __init__(self, dens, grad, kin, grid):
        """
        Parameters
        ----------
        dens : array_like
            Density evaluated on `cube.points`.
        grad
        kin
        grid
        """
        if dens.shape != (grid.npoints,):
            raise ValueError("Arguments dens should have the same size as grid.npoints!")
        if grad.shape != (grid.npoints, 3):
            raise ValueError("Arguments grad should have the same size as grid.npoints!")
        if kin.shape != (grid.shape,):
            raise ValueError("Arguments kin should have the same size as grid.npoints!")
        self._vals = DensityLocalTool(dens, grad, None, kin).electron_localization_function
        self._grid = grid
        self._basins = None
        self._topology = None

    @classmethod
    def from_file(cls, filename, spin='ab', index=None, grid=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        spin
        index
        grid : instance of `CubeGen`, default=None
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        molecule = Molecule.from_file(filename)
        return cls.from_molecule(molecule, spin=spin, index=index, grid=grid)

    @classmethod
    def from_molecule(cls, molecule, spin='ab', index=None, grid=None):
        """Initialize class from ``Molecule`` object.

        Parameters
        ----------
        molecule
        spin
        index
        grid
        """
        # Generate or check cubic grid
        if grid is None:
            cube = CubeGen.from_molecule(molecule.numbers, molecule.pseudo_numbers,
                                         molecule.coordinates, spacing=0.1, threshold=5.0)
        elif not hasattr(grid, 'points'):
            raise ValueError('Argument grid should have "points" attribute!')

        # Compute density, gradient & kinetic energy density on cubic grid
        dens = molecule.compute_density(cube.points, spin=spin, index=index)
        grad = molecule.compute_gradient(cube.points, spin=spin, index=index)
        kin = molecule.compute_kinetic_energy_density(cube.points, spin=spin, index=index)

        return cls(dens, grad, kin, cube)

    @property
    def values(self):
        """Electron Localization Function (ELF) evaluated on grid points."""
        return self._vals

    @property
    def basins(self):
        """Electron Localization Function (ELF) basin values of grid points."""
        if not self._basins:
            self._basins = self._compute_basins()
        return self._basins

    @property
    def topology(self):
        """Electron Localization Function (ELF) topology."""
        if not self._topology:
            self._topology = self._compute_topology()
        return self._topology

    def _compute_basins(self):
        raise NotImplementedError

    def _compute_topology(self):
        raise NotImplementedError

    def dump_isosurface_files(self, filename, isosurf=0.8):
        """
        Parameters
        ----------
        filename
        isosurf
        """
        if not isinstance(self._grid, CubeGen):
            raise ValueError("")
        # dump ELF cube files
        self._grid.dump_cube(filename + '-elf.cube', self.values)
        # write VMD script for visualization
        print_vmd_script_isosurface(filename + '.vmd', filename + '-elf.cube', isosurf=isosurf)
