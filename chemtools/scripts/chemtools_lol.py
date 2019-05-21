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
# pragma pylint: disable=invalid-name
"""Localized Orbital Locator (LOL) Script."""


from chemtools import Molecule, UniformGrid, LOL


__all__ = ['parse_args_lol', 'main_lol']

# description message
lol_dest = """
Visualize Localized Orbital Locator (LOL) using VMD package.

The generated files include:
  output.vmd             The VMD script.
  output-lol.cube        The LOL cube file.

If VMD is setup on your system, you can visualize LOL with the command below:
    $ vmd -e output.vmd
For instruction on how to open the script from the VMD interactive environment,
please refer to ChemTools website.

Note: The output.vmd script requires output-lol.cube to plot LOL
      in VMD software (they files should be all in the same directory).
"""


def parse_args_lol(subparser):
    """Parse command-line arguments for computing LOL."""

    # required arguments
    subparser.add_argument(
        'fname',
        help='wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')

    subparser.add_argument(
        'output', help='name of generated cube file and vmd script.')

    # optional arguments
    subparser.add_argument(
        '--trans',
        default='inverse_rational',
        choices=['inverse_rational', 'inverse_hyperbolic'],
        type=str,
        help='type of transformation applied to LOL ratio. [default=%(default)s]')

    subparser.add_argument(
        '--trans_k',
        default=1,
        type=int,
        help='parameter k of transformation. [default=%(default)s]')

    subparser.add_argument(
        '--trans_a',
        default=1,
        type=int,
        help='parameter a of transformation. [default=%(default)s]')

    subparser.add_argument(
        '--cube',
        default='0.1,2.0',
        type=str,
        metavar='N',
        help='specify the cubic grid used for visualizing LOL. '
        'This can be either a cube file with .cube extension, or a '
        'user-defined cubic grid specified by spacing and extension parameters '
        'separated by a comma. For example, 0.2,5.0 which specifies 0.2 a.u. '
        'distance between grid points, and 5.0 a.u. extension of cubic grid '
        'on each side of the molecule. This cube is used for evaluating '
        'LOL and visualizing it using VMD program. [default=%(default)s]')

    subparser.add_argument(
        '--isosurface',
        default=0.5,
        type=float,
        help='iso-surface value of LOL to visualize. [default=%(default)s]')

    subparser.add_argument(
        '--denscut',
        default=0.0005,
        type=float,
        help='the LOL value of points with density < denscut is set to zero. '
             '[default=%(default)s]')


def main_lol(args):
    """Build LOL model and dump VMD script and cube files for visualizing LOL."""
    # load molecule
    mol = Molecule.from_file(args.fname)

    # make cubic grid
    if args.cube.endswith('.cube'):
        # load cube file
        cube = UniformGrid.from_cube(args.cube)
    elif len(args.cube.split(',')) == 2:
        # make a cubic grid
        spacing, extension = [float(item) for item in args.cube.split(',')]
        cube = UniformGrid.from_molecule(mol, spacing=spacing, extension=extension, rotate=True)
    else:
        raise ValueError('Argument cube={0} is not recognized!'.format(args.cube))

    # build LOL model
    lol = LOL.from_molecule(mol, grid=cube, trans=args.trans, trans_k=args.trans_k,
                            trans_a=args.trans_a, denscut=args.denscut)

    # dump file/script for visualizing LOL
    lol.generate_scripts(args.output, isosurf=args.isosurface)
