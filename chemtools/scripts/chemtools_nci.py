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
"""Non-Covalent Interactions (NCI) Script."""


from chemtools.toolbox.interactions import NCI
from chemtools.scripts.common import help_cube, load_molecule_and_grid


description_nci = """
Visualize Non-Covalent Interactions (NCI) using VMD package.

The generated files include:
  output.vmd             The VMD script.
  output-dens.cube       The signed density cube file.
  output-grad.cube       The reduced density gradient cube file.
  output.png             The reduced density gradient vs. signed density plot.

The values of signed density (density multiplied by the sign of 2nd eigenvalue of
Hessian) are multiplied by 100.0 when being recorded in cube file. Similar to NCIPlot
program, this is used for coloring reduced density gradient iso-surface(s).
The values of reduced density gradient are masked using the given denscut argument
before being recorded in cube file. More specifically, similar to NCIPlot program,
the reduced density gradient value of points for which density > denscut will be
set to 100.0 to have VMD only display reduced density gradient iso-surface(s) for
region with density < denscut.
"""


def parse_args_nci(subparser):
    """Parse command-line arguments for computing NCI."""

    # required arguments
    subparser.add_argument(
        "fname",
        help="wave-function file.")

    subparser.add_argument(
        "-o", "--output",
        default=None,
        type=str,
        help="name of generated output files. By default, it is derived from fname.")

    # optional arguments
    subparser.add_argument(
        "-c", "--cube",
        default="0.1,2.0",
        type=str,
        metavar="",
        help=help_cube)

    subparser.add_argument(
        "--plot",
        default=True,
        action="store_true",
        help="plot reduced density gradient vs. signed density (density"
             "multiplied by sign of hessian's 2nd eigenvalue). This plot is not affected "
             "by the value of denscut argument. [default=%(default)s]")

    subparser.add_argument(
        "-i", "--isosurface",
        default=0.5,
        type=float,
        metavar="",
        help="iso-surface value of reduced density gradient (RDG) to visualize. "
             "[default=%(default)s]")

    subparser.add_argument(
        "-d", "--denscut",
        default=0.05,
        type=float,
        metavar="",
        help="density cutoff used in visualizing reduced density gradient (RDG) "
             "iso-surfaces, and dumping reduced density gradient cube file. "
             "Similar to NCIPlot program, reduced density gradient of points with "
             "density > denscut will be set to 100.0 in the corresponding cube "
             "file. This triggers the VMD to only display reduced density gradient "
             "iso-surface(s) in regions for which density < denscut. "
             "For visualizing all reduced density gradient (RDG) iso-surfaces, "
             "disregarding of density value, set this argument to inf or infinity. "
             "[default=%(default)s]")

    subparser.add_argument(
        "--color",
        default="b",
        type=str,
        help="color of reduced density gradient vs. signed density scatter plot."
             " [default=%(default)s]")


def main_nci(args):
    """Build NCI model and dump VMD script and cube files for visualizing NCI with VMD."""
    # load molecule & cubic grid
    mol, cube = load_molecule_and_grid(args.fname, args.cube)

    # build model
    nci = NCI.from_molecule(mol, grid=cube)

    # dump files for visualization
    output = args.output
    if output is None:
        output = args.fname.split(".")[0]
    nci.generate_scripts(output, isosurf=args.isosurface, denscut=args.denscut)

    # plot reduced density gradient vs. signed density
    if args.plot:
        nci.generate_plot(output, color=args.color)
