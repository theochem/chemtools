#!/usr/bin/env python
# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
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


import argparse
from chemtools import NCI


def parse_args_nci():
    '''
    Parse command-line arguments for computing NCI.
    '''
    parser = argparse.ArgumentParser(prog='chemtools-nci.py',
        description='Generate density & reduced density gradient cube files, as well'
                    'as a VMD (Visual Molecular Dynamics) script to visualize non-covalnet'
                    'inteactions (NCI).\n'
                    'This script will generate 3 files which includes:'
                    'filename_output-dens.cube, filename_output-grad.cube, and filename_output.vmd\n'
                    'If VMD is setup on your system, you can visualize NCI with the command below:\n'
                    '$ vmd -e filename_output.vmd\n'
                    'Note: The filename_output.vmd scripts, requires filename_output-dens.cube & '
                    'filename_output-grad.cube to plot NCI.'
                    'Note: The generated VMD script is the same as the NCIPlot software version 1.0.'
                                     )
    # parser.add_argument('-V', '--version', action='version',
    #     version="%%(prog)s (ChemTools version %s)" % __version__)

    parser.add_argument('filename_wfn',
                        help='Wave-function file. Supported formats: fchk, mkl, molden.input, wfn.')
    parser.add_argument('filename_output',
                        help='Name of generated cube files and vmd script.')

    return parser.parse_args()


def main_nci():
    '''
    Build NCI model using given command-line settings, and dump files/scripts for
    visualizing NCI with VMD.
    '''
    # Parse command-line arguments
    args = parse_args_nci()

    # Build NCI model using default settings
    nci = NCI.from_file(args.filename_wfn)

    # Dump files/scripts for visualizing NCI
    nci.dump_files(args.filename_output)


if __name__ == '__main__':
    main_nci()
