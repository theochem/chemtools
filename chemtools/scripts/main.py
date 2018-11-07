"""The entry point script main."""

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
    'ccpt-l': main_conceptual_local,
    'ccpt-g': main_conceptual_global,
}


def parse_args_chemtools():
    """Parse entry points arguments for chemtools functionality."""
    description = ("""ChemTools command-line tools""")
    parser = argparse.ArgumentParser(prog='chemtools', description=description)

    # main parser to handle basic command and hel function
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version="{} (ChemTools version {})".format(parser.prog, __version__))

    # command parser, stored in parser.command
    subparser = parser.add_subparsers(
        metavar="<Commands>", help='<Functions>', dest='command')

    # sub parser for nci functions
    parser_nci = subparser.add_parser(
        'nci', help='visualize non-covalent interactions')
    parse_args_nci(parser_nci)

    # sub parser for ccpt-1 functions
    parser_nci = subparser.add_parser(
        'ccpt-l', help='compute local conceptual DFT indicators')
    parse_args_local(parser_nci)

    # sub parser for ccpt-g functions
    parser_nci = subparser.add_parser(
        'ccpt-g', help='compute global conceptual DFT indicators')
    parse_args_global(parser_nci)

    return parser.parse_args()


def main():
    """Entry point function for Chemtools."""
    arg = parse_args_chemtools()  # parse all variables for each functions
    main_fun = SCRIPT_MAIN[arg.command]  # call the main executable function
    main_fun(arg)  # run the function
