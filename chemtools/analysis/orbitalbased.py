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

from horton import IOData
from chemtools.toolbox.orbitalbased import OrbitalLocalTool

__all__ = ['OrbitalAnalysis']


class OrbitalAnalysis(OrbitalLocalTool):
    """Class for orbital-based analysis."""

    def __init__(self, points, obasis, exp_alpha, exp_beta=None):
        super(OrbitalAnalysis, self).__init__(points, obasis, exp_alpha, exp_beta=exp_beta)

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
        # case of one file not given as a list
        if isinstance(filename, (str, unicode)):
            mol = IOData.from_file(filename)
            if hasattr(mol, 'exp_beta'):
                return cls(points, mol.obasis, mol.exp_alpha, mol.exp_beta)
            else:
                return cls(points, mol.obasis, mol.exp_alpha)
        # case of list of file(s)
        for _ in filename:
            raise ValueError('Multiple files are not supported')

    @classmethod
    def from_iodata(cls, iodata, points):
        """
        Initialize class from `IOData` object.

        Parameters
        ----------
        iodata : `IOData`
            Instance of `IOData`.
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        # check if iodata has exp_beta
        if hasattr(iodata, 'exp_beta'):
            return cls(points, iodata.obasis, iodata.exp_alpha, iodata.exp_beta)
        else:
            return cls(points, iodata.obasis, iodata.exp_alpha)
