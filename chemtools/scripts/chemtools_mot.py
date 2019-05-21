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
"""Molecular Orbital Theory (MOT) Script."""


import logging

from chemtools import Molecule, UniformGrid, MOTBasedTool

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


__all__ = ['parse_args_mot', 'main_mot']

mot_desp = """
Generate a VMD (Visual Molecular Dynamics) script as well as cube file required
for visualizing molecular orbitals with VMD package.

The generated files include:
  output.vmd                 The VMD script.
  output_mo{index}.cube      The MO cube file.

If VMD is setup on your system, you can visualize MO with the command below:
    $ vmd -e output.vmd
For instruction on how to open the script from the VMD interactive environment,
please refer to ChemTools website.

Note: The output.vmd script requires output_mo{index}.cube to plot MO
      in VMD software (they files should be all in the same directory).
"""


def parse_args_mot(subparser):
    """Parse command-line arguments for computing MO."""
    # description message

    # required arguments
    subparser.add_argument(
        'fname',
        help='wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')

    subparser.add_argument(
        'output', help='name of generated cube file and vmd script.')

    # optional arguments
    subparser.add_argument(
        '--spin',
        default='a',
        choices=['a', 'b'],
        type=str,
        help='type of occupied spin orbitals to visualize. [default=%(default)s]')

    subparser.add_argument(
        '--index',
        default=None,
        type=str,
        help='index of spin orbital to visualize represented by comma separated integers.'
             'If None, files for generating all occupied molecular orbitals are generated.'
             ' [default=%(default)s]')

    subparser.add_argument(
        '--cube',
        default='0.2,5.0',
        type=str,
        metavar='N',
        help='specify the cubic grid used for visualizing MO. '
        'This can be either a cube file with .cube extension, or a '
        'user-defined cubic grid specified by spacing and extension parameters '
        'separated by a comma. For example, 0.2,5.0 which specifies 0.2 a.u. '
        'distance between grid points, and 5.0 a.u. extension of cubic grid '
        'on each side of the molecule. This cube is used for evaluating '
        'MO and visualizing it using VMD program. [default=%(default)s]')

    subparser.add_argument(
        '--isosurface',
        default=0.05,
        type=float,
        help='iso-surface value of MO to visualize. [default=%(default)s]')


def main_mot(args):
    """Build MOTBasedTool model and dump VMD script and cube files for visualizing MO."""
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

    if args.index is not None:
        index = [int(item) for item in args.index.split(',')]
        if len(index) == 1:
            index = index[0]
    else:
        index = None

    # build MOT model
    mot = MOTBasedTool.from_molecule(mol)

    # print logging message
    logging.info('')
    logging.info('Initialized : {0}'.format(mot))
    logging.info('# of basis           : {0}'.format(mot.nbasis))
    logging.info('# of a & b electrons : {0}'.format(mot.nelectrons))
    logging.info('')
    logging.info('Index  of a & b LUMO : {0}'.format(mot.lumo_index))
    logging.info('Energy of a & b LUMO : {0}'.format(mot.lumo_energy))
    logging.info('')
    logging.info('Index  of a & b HOMO : {0}'.format(mot.homo_index))
    logging.info('Energy of a & b HOMO : {0}'.format(mot.homo_energy))
    logging.info('')

    # dump file/script for visualizing MOT
    mot.generate_scripts(args.output, spin=args.spin, index=index, grid=cube,
                         isosurf=args.isosurface)
