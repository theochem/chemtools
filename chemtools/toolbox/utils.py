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
"""Utility Functions of Toolbox Module."""


from chemtools.utils.molecule import BaseMolecule


__all__ = ["check_arg_molecule"]


def check_arg_molecule(molecule):
    """Return molecule argument after checking.

    Parameters
    ----------
    molecule : BaseMolecule or Sequence of BaseMolecule
        Instance of BaseMolecule class, or sequence of BaseMolecule class instances.

    Returns
    -------
    molecule : BaseMolecule or Sequence of BaseMolecule
        Instance of BaseMolecule or Sequence of BaseMolecule with more than one instance.
    """
    if isinstance(molecule, BaseMolecule):
        return molecule
    if hasattr(molecule, "__iter__") and len(molecule) == 1:
        # sequence of just one molecule
        return molecule[0]
    return molecule
