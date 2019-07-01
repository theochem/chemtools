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
"""Electrostatic Potential (ESP) Script."""


import logging

from chemtools import print_vmd_script_isosurface

from chemtools.scripts.common import help_cube, load_molecule_and_grid


logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


__all__ = ['parse_args_esp', 'main_esp']

esp_desp = """
Visualize Electrostatic Potential (ESP) on electron density iso-surface using VMD package.

The generated files include:
  output.vmd          The VMD script.
  output_esp.cube     The ESP cube file.
  output_dens.cube    The density cube file.
"""


def parse_args_esp(subparser):
    """Parse command-line arguments for computing ESP."""
    # required arguments
    subparser.add_argument(
        'fname',
        help='wave-function file.')

    # optional arguments
    subparser.add_argument(
        "-o", "--output",
        default=None,
        type=str,
        metavar="",
        help='name of generated output files. By default, it is derived from fname.')

    # subparser.add_argument(
    #     '--spin',
    #     default='a',
    #     choices=['a', 'b'],
    #     type=str,
    #     help='type of occupied spin orbitals to visualize. [default=%(default)s]')
    #
    # subparser.add_argument(
    #     '--index',
    #     default=None,
    #     type=str,
    #     help='index of spin orbital to visualize represented by comma separated integers.'
    #          'If None, files for generating all occupied molecular orbitals are generated.'
    #          ' [default=%(default)s]')

    subparser.add_argument(
        "-c", "--cube",
        default='0.5,5.0',
        type=str,
        metavar="",
        help=help_cube)

    subparser.add_argument(
        "-i", "--isosurface",
        default=0.002,
        type=float,
        metavar="",
        help="iso-surface value of electron density to visualize. [default=%(default)s]")

    subparser.add_argument(
        '--scalemin',
        default=-0.02,
        type=float,
        metavar="",
        help='minimum value of ESP to color on the electron density iso-surface. '
             '[default=%(default)s]')

    subparser.add_argument(
        '--scalemax',
        default=0.04,
        type=float,
        metavar="",
        help='maximum value of ESP to color on the electron density iso-surface. '
             '[default=%(default)s]')


def main_esp(args):
    """Generate VMD script and cube files for visualizing ESP on electron density iso-surface."""
    # load molecule & cubic grid
    mol, cube = load_molecule_and_grid(args.fname, args.cube)

    # dump cube files & script for visualization
    output = args.output
    if output is None:
        output = args.fname.split(".")[0]
    espname = output + '_esp.cube'
    rhoname = output + '_dens.cube'
    vmdname = output + '.vmd'

    cube.generate_cube(rhoname, mol.compute_density(cube.points))
    cube.generate_cube(espname, mol.compute_esp(cube.points))
    print_vmd_script_isosurface(vmdname, rhoname, colorfile=espname, isosurf=args.isosurface,
                                scalemin=args.scalemin, scalemax=args.scalemax)
