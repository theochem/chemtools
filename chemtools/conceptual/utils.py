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
"""The Utility Functions of Conceptual Module."""


__all__ = []


def check_dict_energy(dict_energy):
    """Check sanity of number of electrons and corresponding energy values.

    Parameters
    ----------
    dict_energy : dict
        Dictionary of number of electrons (keys) and corresponding energy (values).

    Returns
    -------
    n_ref : float

    energy_m : float

    energy_0 : float

    energy_p : float

    """
    # check length of dictionary & its keys
    if len(dict_energy) != 3 or not all([key >= 0 for key in dict_energy.keys()]):
        raise ValueError('The energy model requires 3 energy values corresponding '
                         'to positive number of electrons!')
    # find reference number of electrons
    n_ref = sorted(dict_energy.keys())[1]
    if n_ref < 1:
        raise ValueError('The n_ref cannot be less than one! Given n_ref={0}'.format(n_ref))
    # check number of electrons differ by one
    if sorted(dict_energy.keys()) != [n_ref - 1, n_ref, n_ref + 1]:
        raise ValueError('Number of electrons should differ by one!')
    # get energy values
    energy_m, energy_0, energy_p = [dict_energy[n] for n in sorted(dict_energy.keys())]
    return n_ref, energy_m, energy_0, energy_p
