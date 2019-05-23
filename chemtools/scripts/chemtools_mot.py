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
import numpy as np

from chemtools import Molecule, UniformGrid, MOTBasedTool

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')


__all__ = ['parse_args_mot', 'main_mot']

mot_desp = """
Visualize Molecular Orbitals (MO) using VMD package.

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
        '--output',
        default=None,
        help='name of generated cube file and vmd script. If None, the output name in derived'
             'from fname. [default=%(default)s]')

    # optional arguments
    subparser.add_argument(
        '--info',
        action='store_true',
        default=False,
        help='print basic information on molecule and wave-function. [default=%(default)s]')

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

    hia, hib = np.array(mol.homo_index) - 1
    lia, lib = np.array(mol.lumo_index) - 1
    ea, eb = mol.orbital_energy
    print('')
    print('File: {0}'.format(args.fname))
    # print('Charge      : % 5f' % np.sum(mol.numbers) - np.sum(ne))
    # print('Multiplicity: % 5d' % np.max(ne) - np.min(ne) + 1)
    print('')
    print('Atomic number and coordinates:')
    for index, num in enumerate(mol.numbers):
        coord = mol.coordinates[index, :]
        print('% 2i   %10.6f   %10.6f   %10.6f' % (num, coord[0], coord[1], coord[2]))
    print('')
    print('Information on alpha & beta electrons:')
    print('# electrons  :  % 3.3f       % 3.3f' % mol.nelectrons)
    print('HOMO index   : % 3d        % 5d' % mol.homo_index)
    print('')
    print('LUMO+2 index : %10.6f   %10.6f' % (ea[lia + 2], eb[lib + 2]))
    print('LUMO+1 energy: %10.6f   %10.6f' % (ea[lia + 1], eb[lib + 1]))
    print('LUMO   energy: %10.6f   %10.6f' % mol.lumo_energy)
    print('HOMO   energy: %10.6f   %10.6f' % mol.homo_energy)
    print('HOMO-1 energy: %10.6f   %10.6f' % (ea[hia - 1], eb[hib - 1]))
    print('HOMO-2 energy: %10.6f   %10.6f' % (ea[hia - 2], eb[hib - 2]))
    print('')

    if args.info:
        return

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
    # logging.info('')
    # logging.info('Initialized : {0}'.format(mot))
    # logging.info('# of basis           : {0}'.format(mot.nbasis))
    # logging.info('# of a & b electrons : {0}'.format(mot.nelectrons))
    # logging.info('')
    # logging.info('Index  of a & b LUMO : {0}'.format(mot.lumo_index))
    # logging.info('Energy of a & b LUMO : {0}'.format(mot.lumo_energy))
    # logging.info('')
    # logging.info('Index  of a & b HOMO : {0}'.format(mot.homo_index))
    # logging.info('Energy of a & b HOMO : {0}'.format(mot.homo_energy))
    # logging.info('')

    # dump file/script for visualizing MOT
    if args.output is None:
        args.output = args.fname.rsplit('.')[0]
    mot.generate_scripts(args.output, spin=args.spin, index=index, grid=cube,
                         isosurf=args.isosurface)
