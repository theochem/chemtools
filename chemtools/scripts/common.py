#!/usr/bin/env python
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
"""Common utility for scripts."""


import numpy as np

from chemtools.wrappers.molecule import Molecule
from chemtools.utils.cube import UniformGrid


help_cube = """
cubic grid used for evaluation and visualization.
This can be either a cube file with .cube extension, or a user-defined
cubic grid specified by comma-separated spacing and extension values,
e.g., 0.2,5.0 specifies 0.2 a.u. distance between grid points, and 5.0 a.u.
extension on each side of molecule. [default=%(default)s]
"""


def load_molecule_and_grid(fname, cube):
    """Return instances of molecule and uniform cubic grid.

    Parameters
    ----------
    fname : str
        Path to wave-function file.
    cube : str
       Uniform cubic grid specifications.

    """
    # load molecule
    mol = Molecule.from_file(fname)

    if cube.endswith(".cube"):
        # load & check cube file
        cube = UniformGrid.from_cube(cube)
        if np.allclose(mol.numbers, cube.numbers):
            raise ValueError("Atomic number in {0} & {1} should be the same!".format(fname, cube))
        if np.allclose(mol.coordinates, cube.coordinates):
            raise ValueError(
                "Atomic coordinates in {0} & {1} should be the same!".format(cube.fname, cube.cube)
            )
    elif len(cube.split(",")) == 2:
        # make a cubic grid
        spacing, extension = [float(item) for item in cube.split(",")]
        cube = UniformGrid.from_molecule(mol, spacing=spacing, extension=extension, rotate=True)

    else:
        raise ValueError("Argument cube={0} is not recognized!".format(cube))

    return mol, cube
