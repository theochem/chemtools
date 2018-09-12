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
"""Module for Density-Based and Orbital-Based Analysis of Quantum Chemistry Output Files.

This modules contains wrappers which take outputs of quantum chemistry software and
compute various descriptive tools based on the density and orbital information.
"""

from chemtools.wrappers.molecule import Molecule
from chemtools.orbtools.orbitalbased import OrbitalLocalTool

__all__ = ['OrbitalAnalysis']


class OrbitalAnalysis(OrbitalLocalTool):
    """Class for orbital-based analysis."""

    def __init__(self, molecule, points):
        super(OrbitalAnalysis, self).__init__(molecule, points)

    @classmethod
    def from_file(cls, filename, points):
        """
        Initialize class from file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        molecule = Molecule.from_file(filename)
        return cls.from_molecule(molecule, points)

    @classmethod
    def from_molecule(cls, molecule, points):
        """
        Initialize class from `Molecule` object.

        Parameters
        ----------
        molecule : ``Molecule``
            Instance of ``Molecule``.
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        return cls(molecule, points)
