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
# pragma pylint: disable=invalid-name
"""Entry Point Main Script."""


import argparse

from chemtools import __version__
from chemtools.scripts.chemtools_conceptual import (
    main_conceptual_global, main_conceptual_local, parse_args_global,
    parse_args_local)
from chemtools.scripts.chemtools_nci import main_nci, parse_args_nci


__all__ = ['main']

# basic swtich dictionary for storing all main callable function and subparser
SCRIPT_MAIN = {
    'nci': main_nci,
    'lcdft': main_conceptual_local,
    'gcdft': main_conceptual_global,
}


def parse_args_chemtools():
    """Parse entry points arguments for chemtools functionality."""
    description = """ChemTools command-line tools"""
    parser = argparse.ArgumentParser(prog='chemtools', description=description)

    # main parser to handle basic command and hel function
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version="{} (ChemTools version {})".format(parser.prog, __version__))

    # command parser, stored in parser.command
    subparser = parser.add_subparsers(metavar="<Commands>", help='<Functions>', dest='command')

    # sub parser for nci functions
    parser_nci = subparser.add_parser('nci', help='visualize non-covalent interactions')
    parse_args_nci(parser_nci)

    # sub parser for lcdft functions
    parser_nci = subparser.add_parser('lcdft', help='compute local conceptual DFT indicators')
    parse_args_local(parser_nci)

    # sub parser for gcdft functions
    parser_nci = subparser.add_parser('gcdft', help='compute global conceptual DFT indicators')
    parse_args_global(parser_nci)

    return parser.parse_args()


def main():
    """Entry point function for Chemtools."""
    arg = parse_args_chemtools()  # parse all variables for each functions
    main_fun = SCRIPT_MAIN[arg.command]  # call the main executable function
    main_fun(arg)  # run the function
