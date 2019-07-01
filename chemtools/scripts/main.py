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
"""Entry Point Main Script."""

import argparse

from chemtools import __version__
from chemtools.scripts.chemtools_conceptual import (
    main_conceptual_global,
    main_conceptual_local,
    main_conceptual_condensed,
    parse_args_global,
    parse_args_local,
    parse_args_condensed,
    description_global,
)
from chemtools.scripts.chemtools_nci import main_nci, parse_args_nci, description_nci
from chemtools.scripts.chemtools_elf import main_elf, parse_args_elf, description_elf
from chemtools.scripts.chemtools_lol import main_lol, parse_args_lol, description_lol
from chemtools.scripts.chemtools_mot import main_mot, parse_args_mot, description_mot
from chemtools.scripts.chemtools_esp import main_esp, parse_args_esp, description_esp

from argparse import RawDescriptionHelpFormatter


__all__ = ["main"]


# basic swtich dictionary for storing all main callable function and subparser
SCRIPT_MAIN = {
    "mot": main_mot,
    "esp": main_esp,
    "nci": main_nci,
    "elf": main_elf,
    "lol": main_lol,
    "gcdft": main_conceptual_global,
    "lcdft": main_conceptual_local,
    "ccdft": main_conceptual_condensed,
}


def parse_args_chemtools():
    """Parse entry points arguments for chemtools functionality."""
    description = """ChemTools command-line tools"""
    parser = argparse.ArgumentParser(prog="chemtools", description=description)

    # main parser to handle basic command and hel function
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="{} (ChemTools version {})".format(parser.prog, __version__),
    )

    # command parser, stored in parser.command
    subparser = parser.add_subparsers(
        metavar="<Commands>", help="<Functions>", dest="command"
    )

    parser_mot = subparser.add_parser(
        "mot",
        help="Molecular Orbital Theory (MOT).",
        description=description_mot,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_mot(parser_mot)

    parser_esp = subparser.add_parser(
        "esp",
        help="Electrostatic Potential (ESP).",
        description=description_esp,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_esp(parser_esp)

    parser_nci = subparser.add_parser(
        "nci",
        help="Non-Covalent Interactions (NCI).",
        description=description_nci,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_nci(parser_nci)

    parser_elf = subparser.add_parser(
        "elf",
        help="Electron Localization Function (ELF).",
        description=description_elf,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_elf(parser_elf)

    parser_lol = subparser.add_parser(
        "lol",
        help="Localized Orbital Locator (LOL).",
        description=description_lol,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_lol(parser_lol)

    # sub parser for gcdft functions
    parser_g = subparser.add_parser(
        "gcdft",
        help="Global Conceptual Density Functional Theory.",
        description=description_global,
        formatter_class=RawDescriptionHelpFormatter,
    )
    parse_args_global(parser_g)

    # sub parser for lcdft functions
    parser_l = subparser.add_parser(
        "lcdft",
        help="Local Conceptual DFT.",
        description='Local Conceptual Density Functional Theory.',
    )
    parse_args_local(parser_l)

    parser_c = subparser.add_parser(
        "ccdft",
        help="Condensed Conceptual DFT.",
        description='Condensed Conceptual Density Functional Theory.',
    )
    parse_args_condensed(parser_c)

    return parser.parse_args()


def main():
    """Entry point function for Chemtools."""
    arg = parse_args_chemtools()  # parse all variables for each functions
    main_fun = SCRIPT_MAIN[arg.command]  # call the main executable function
    main_fun(arg)  # run the function


if __name__ == "__main__":
    import logging

    try:
        import gooey

        # this is command is equivalent to use decorator
        gooey.Gooey(main)()
    except ImportError:
        logging.warning(
            "Package 'gooey' is needed for GUI command-line tools. Please"
            " install\nit with your package manager."
        )
