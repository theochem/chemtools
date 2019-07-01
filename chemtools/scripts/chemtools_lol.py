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


from chemtools import LOL
from chemtools.scripts.common import help_cube, load_molecule_and_grid


__all__ = ['parse_args_lol', 'main_lol']

# description message
lol_dest = """
Visualize Localized Orbital Locator (LOL) using VMD package.

The generated files include:
  output.vmd             The VMD script.
  output-lol.cube        The LOL cube file.
"""


def parse_args_lol(subparser):
    """Parse command-line arguments for computing LOL."""

    # required arguments
    subparser.add_argument(
        'fname',
        help='wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')

    # optional arguments
    subparser.add_argument(
        "-o", "--output",
        default=None,
        type=str,
        help='name of generated output files. By default, it is derived from fname.')

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
        metavar="",
        help=help_cube)

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
    # load molecule & cubic grid
    mol, cube = load_molecule_and_grid(args.fname, args.cube)

    # build model
    lol = LOL.from_molecule(mol, grid=cube, trans=args.trans, trans_k=args.trans_k,
                            trans_a=args.trans_a, denscut=args.denscut)

    # dump cube file/script for visualization
    output = args.output
    if output is None:
        output = args.fname.split(".")[0]
    lol.generate_scripts(output, isosurf=args.isosurface)
