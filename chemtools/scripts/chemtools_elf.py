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


from chemtools.toolbox.interactions import ELF
from chemtools.scripts.common import help_cube, load_molecule_and_grid


description_elf = """
Visualize Electron Localization Function (ELF) using VMD package.

The generated files include:
  output.vmd             The VMD script.
  output-elf.cube        The ELF cube file.
"""


def parse_args_elf(subparser):
    """Parse command-line arguments for computing ELF."""

    # required arguments
    subparser.add_argument(
        "fname",
        help="wave-function file.")

    # optional arguments
    subparser.add_argument(
        "-o", "--output",
        default=None,
        type=str,
        metavar="",
        help="name of generated output files. By default, it is derived from fname.")

    subparser.add_argument(
        "-t", "--trans",
        default="rational",
        choices=["rational", "hyperbolic"],
        type=str,
        help="type of transformation applied to ELF ratio. [default=%(default)s]")

    subparser.add_argument(
        "-k", "--trans_k",
        default=2,
        type=int,
        metavar="",
        help="parameter k of transformation. [default=%(default)s]")

    subparser.add_argument(
        "-a", "--trans_a",
        default=1,
        type=int,
        metavar="",
        help="parameter a of transformation. [default=%(default)s]")

    subparser.add_argument(
        "-c", "--cube",
        default="0.1,2.0",
        type=str,
        metavar="",
        help=help_cube)

    subparser.add_argument(
        "-i", "--isosurface",
        default=0.8,
        type=float,
        metavar="",
        help="iso-surface value of ELF to visualize. [default=%(default)s]")

    subparser.add_argument(
        "-d", "--denscut",
        default=0.0005,
        type=float,
        metavar="",
        help="the ELF value of points with electron density < denscut is set to zero. "
             "[default=%(default)s]")


def main_elf(args):
    """Build ELF model and dump VMD script and cube files for visualizing ELF."""
    # load molecule & cubic grid
    mol, cube = load_molecule_and_grid(args.fname, args.cube)

    # build model
    elf = ELF.from_molecule(mol, grid=cube, trans=args.trans, trans_k=args.trans_k,
                            trans_a=args.trans_a, denscut=args.denscut)

    # dump files for visualization
    output = args.output
    if output is None:
        output = args.fname.split(".")[0]
    elf.generate_scripts(output, isosurf=args.isosurface)
