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


from numpy.ma import masked_less

from chemtools.utils.cube import CubeGen
from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensGradBasedTool
from chemtools.outputs import print_vmd_script_isosurface


__all__ = ['ELF']


class ELF(object):
    r"""Electron Localization Function (ELF).

    Introduced by Becke and Edgecombe:

    .. math::
       ELF (\mathbf{r}) =
            \frac{1}{\left( 1 + \left(\frac{D_{\sigma}(\mathbf{r})}
            {D_{\sigma}^0 (\mathbf{r})} \right)^2\right)}

    with XXX, XXX, and positive definite kinetic energy density defined as, respectively,

    .. math::
        D_{\sigma} (\mathbf{r}) &= \tau_{\sigma} (\mathbf{r}) -
           \frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}}

       D_{\sigma}^0 (\mathbf{r}) &=
          \frac{3}{5} (6 \pi^2)^{2/3} \rho_{\sigma}^{5/3} (\mathbf{r})

       \tau_{\sigma} (\mathbf{r}) =
             \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2
    """

    def __init__(self, dens, grad, kin, grid=None):
        """Initialize ELF class.

        Parameters
        ----------
        dens : array_like
            Density evaluated on `cube.points`.
        grad
        kin
        grid
        """
        # if dens.shape != (grid.npoints,):
        #     raise ValueError("Arguments dens should have the same size as grid.npoints!")
        # if grad.shape != (grid.npoints, 3):
        #     raise ValueError("Arguments grad should have the same size as grid.npoints!")
        # if kin.shape != (grid.shape,):
        #     raise ValueError("Arguments kin should have the same size as grid.npoints!")
        self._grid = grid
        self._denstool = DensGradBasedTool(dens, grad)
        # compute elf value
        self._vals = kin - self._denstool.kinetic_energy_density_weizsacker
        self._vals /= masked_less(self._denstool.kinetic_energy_density_thomas_fermi, 1.0e-30)
        self._vals = 1.0 / (1.0 + self._vals**2.0)
        # assign basins and topology attributes
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
        grid : instance of `CubeGen`, optional
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
        # generate cubic grid or check grid
        if grid is None:
            grid = CubeGen.from_molecule(molecule.numbers, molecule.pseudo_numbers,
                                         molecule.coordinates, spacing=0.1, threshold=5.0)
        elif not hasattr(grid, 'points'):
            raise ValueError('Argument grid should have "points" attribute!')

        # compute density, gradient & kinetic energy density on grid
        dens = molecule.compute_density(grid.points, spin=spin, index=index)
        grad = molecule.compute_gradient(grid.points, spin=spin, index=index)
        kin = molecule.compute_kinetic_energy_density(grid.points, spin=spin, index=index)

        return cls(dens, grad, kin, grid)

    @property
    def density(self):
        r"""Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on grid points."""
        return self._denstool.density

    @property
    def values(self):
        r"""Electron localization function :math:`ELF(\mathbf{r})` evaluated on grid points."""
        return self._vals

    @property
    def basins(self):
        """Electron localization function basin values of grid points."""
        if not self._basins:
            self._basins = self._compute_basins()
        return self._basins

    @property
    def topology(self):
        """Electron localization function topology."""
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
