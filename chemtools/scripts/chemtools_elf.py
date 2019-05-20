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
"""Electron Localization Function (ELF) Script."""


from chemtools import Molecule, UniformGrid, ELF


__all__ = ['parse_args_elf', 'main_elf']


def parse_args_elf(subparser):
    """Parse command-line arguments for computing ELF."""
    # description message
    description = """
    Generate a VMD (Visual Molecular Dynamics) script as well as cube file required
    for visualizing Electron Localization Function (ELF) with VMD package.

    The generated files include:
      fname_output.vmd             The VMD script.
      fname_output-elf.cube        The ELF cube file.

    If VMD is setup on your system, you can visualize ELF with the command below:
        $ vmd -e fname_output.vmd
    For instruction on how to open the script from the VMD interactive environment,
    please refer to ChemTools website.

    Note: The fname_output.vmd script requires fname_output-elf.cube to plot ELF
          in VMD software (they files should be all in the same directory).
    """

    # required arguments
    subparser.add_argument(
        'file_wfn',
        help='wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')

    subparser.add_argument(
        'output_name', help='name of generated cube file and vmd script.')

    # optional arguments
    subparser.add_argument(
        '--trans',
        default='rational',
        choices=['rational', 'hyperbolic'],
        type=str,
        help='type of transformation applied to ELF ratio. [default=%(default)s]')

    subparser.add_argument(
        '--trans_k',
        default=2,
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
        help='specify the cubic grid used for visualizing ELF. '
        'This can be either a cube file with .cube extension, or a '
        'user-defined cubic grid specified by spacing and extension parameters '
        'separated by a comma. For example, 0.2,5.0 which specifies 0.2 a.u. '
        'distance between grid points, and 5.0 a.u. extension of cubic grid '
        'on each side of the molecule. This cube is used for evaluating '
        'ELF and visualizing it using VMD program. [default=%(default)s]')

    subparser.add_argument(
        '--isosurface',
        default=0.8,
        type=float,
        help='iso-surface value of ELF to visualize. [default=%(default)s]')

    subparser.add_argument(
        '--denscut',
        default=0.0005,
        type=float,
        help='the ELF value of points with density < denscut is set to zero. '
             '[default=%(default)s]')


def main_elf(args):
    """Build ELF model and dump VMD script and cube files for visualizing ELF."""
    # load molecule
    mol = Molecule.from_file(args.file_wfn)

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

    # build ELF model
    elf = ELF.from_molecule(mol, grid=cube, trans=args.trans, trans_k=args.trans_k,
                            trans_a=args.trans_a, denscut=args.denscut)

    # dump file/script for visualizing ELF
    elf.generate_scripts(args.output_name, isosurf=args.isosurface)
