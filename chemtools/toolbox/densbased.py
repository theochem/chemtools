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
"""Density-Based Local Tools."""


from chemtools.wrappers2.molecule import Molecule
from chemtools.denstools.densbased import DensGradLapKedTool


class DensityLocalTool(DensGradLapKedTool):
    """Density Local Tool Class."""

    def __init__(self, dens, grad, lap, ked):
        r"""Initialize class from arrays.

        Parameters
        ----------
        dens : np.ndarray
            Electron density evaluated on a set of points, :math:`\rho(\mathbf{r})`.
        grad : np.ndarray
            Gradient vector of electron density evaluated on a set of points,
            :math:`\nabla \rho(\mathbf{r})`.
        lap : np.ndarray
            Laplacian of electron density evaluated on a set of points,
            :math:`\nabla^2 \rho(\mathbf{r})`.
        ked : np.ndarray
            Positive-definite or Lagrangian kinetic energy density evaluated on a set of
            points; :math:`\tau_\text{PD} (\mathbf{r})` or :math:`G(\mathbf{r})`.

        """
        super(DensityLocalTool, self).__init__(dens, grad, lap, ked)

    @classmethod
    def from_molecule(cls, molecule, points, spin='ab', index=None):
        r"""Initialize class using instance of `Molecule` and points.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class.
        points : np.ndarray
            The (npoints, 3) array of cartesian coordinates of points.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : sequence, optional
            Sequence of integers representing the index of spin orbitals.

        """
        dens = molecule.compute_density(points, spin, index)
        grad = molecule.compute_gradient(points, spin, index)
        lap = molecule.compute_laplacian(points, spin, index)
        ked = molecule.compute_ked(points, spin, index)
        return cls(dens, grad, lap, ked)

    @classmethod
    def from_file(cls, fname, points, spin='ab', index=None):
        r"""Initialize class from file.

        Parameters
        ----------
        fname : str
            Path to molecule's file.
        points : np.ndarray
            The (npoints, 3) array of cartesian coordinates of points.
        spin : str, optional
            Type of occupied spin orbitals; options are 'a', 'b' & 'ab'.
        index : sequence, optional
            Sequence of integers representing the index of spin orbitals.

        """
        molecule = Molecule.from_file(fname)
        return cls.from_molecule(molecule, points, spin, index)
