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

from chemtools import Molecule, UniformGrid, print_vmd_script_isosurface

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


__all__ = ['parse_args_esp', 'main_esp']

esp_desp = """
Visualize Electrostatic Potential (ESP) on electron density iso-surface with VMD package.

The generated files include:
  output.vmd          The VMD script.
  output_esp.cube     The esp cube file.
  output_dens.cube    The density cube file.

If VMD is setup on your system, you can visualize iso-surface with the command below:
    $ vmd -e output.vmd
For instruction on how to open the script from the VMD interactive environment,
please refer to ChemTools website.

Note: The output.vmd script requires output_esp.cube & output_dens.cube to visualize ESP
      on electron density iso-surface using VMD software (they files should be all in the
      same directory).
"""


def parse_args_esp(subparser):
    """Parse command-line arguments for computing ESP."""
    # required arguments
    subparser.add_argument(
        'fname',
        help='wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')

    subparser.add_argument(
        'output', help='name of generated cube files and vmd script.')

    # optional arguments
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
        '--cube',
        default='0.5,5.0',
        type=str,
        metavar='N',
        help='specify the cubic grid used for visualizing ESP on density iso-surface. '
        'This can be either a cube file with .cube extension, or a '
        'user-defined cubic grid specified by spacing and extension parameters '
        'separated by a comma. For example, 0.2,5.0 which specifies 0.2 a.u. '
        'distance between grid points, and 5.0 a.u. extension of cubic grid '
        'on each side of the molecule. This cube is used for evaluating '
        'ESP and visualizing it using VMD program. [default=%(default)s]')

    subparser.add_argument(
        '--isosurface',
        default=0.002,
        type=float,
        help='electron density iso-surface value to visualize. [default=%(default)s]')

    subparser.add_argument(
        '--scalemin',
        default=-0.02,
        type=float,
        help='minimum value of ESP to color on the electron density iso-surface. '
             '[default=%(default)s]')

    subparser.add_argument(
        '--scalemax',
        default=0.04,
        type=float,
        help='maximum value of ESP to color on the electron density iso-surface. '
             '[default=%(default)s]')


def main_esp(args):
    """Generate VMD script and cube files for visualizing ESP on electron density iso-surface."""
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

    # dump cube files & script for visualization
    espname = args.output + '_esp.cube'
    rhoname = args.output + '_rho.cube'
    vmdname = args.output + '.vmd'

    cube.generate_cube(rhoname, mol.compute_density(cube.points))
    cube.generate_cube(espname, mol.compute_esp(cube.points))
    print_vmd_script_isosurface(vmdname, rhoname, colorfile=espname, isosurf=args.isosurface)
