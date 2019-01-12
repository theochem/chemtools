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
# pragma pylint: disable=invalid-name
"""Kinetic Energy Density Module."""


import numpy as np

from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensityLocalTool


class KineticEnergyDensity(object):
    """Kinetic Energy Density Class."""

    def __init__(self, molecule, points, spin="ab", index=None):
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
        self._points = points
        # compute density, gradient, & kinetic energy density on grid
        dens = molecule.compute_density(self._points, spin, index)
        grad = molecule.compute_gradient(self._points, spin, index)
        ke = molecule.compute_kinetic_energy_density(self._points, spin, index)
        # initialize density-based tools class
        self._denstools = DensityLocalTool(dens, grad, None, ke)

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
        return cls(molecule, points, spin, index)

    @property
    def points(self):
        return self._points

    @property
    def density(self):
        r"""Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on the grid."""
        return self._denstools.density

    @property
    def positive_definite(self):
        """Positive definite kinetic energy density.

        See :py:meth:`DensityLocalTool.kinetic_energy_density
        <chemtools.denstools.densitybased.DensityLocalTool.kinetic_energy_density>`
        attribute.
        """
        return self._denstools.kinetic_energy_density

    @property
    def thomas_fermi(self):
        """Thomas-Fermi kinetic energy density.

        See :py:meth:`DensityLocalTool.kinetic_energy_density_thomas_fermi
        <chemtools.denstools.densitybased.DensityLocalTool.kinetic_energy_density_thomas_fermi>`
        attribute.
        """
        return self._denstools.kinetic_energy_density_thomas_fermi

    @property
    def weizsacker(self):
        """Weizsacker kinetic energy density.

        See :py:meth:`DensityLocalTool.kinetic_energy_density_weizsacker
        <chemtools.denstools.densitybased.DensityLocalTool.kinetic_energy_density_weizsacker>`
        attribute.
        """
        return self._denstools.kinetic_energy_density_weizsacker
