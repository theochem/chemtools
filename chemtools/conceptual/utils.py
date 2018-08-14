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


from horton import log


__all__ = ["check_dict_values", "check_number_electrons"]


def check_dict_values(dict_values):
    """Check sanity of number of electrons and corresponding property values.

    Parameters
    ----------
    dict_values : dict
        Dictionary of number of electrons (keys) and corresponding property value (values).
    """
    # check length of dictionary & its keys
    if len(dict_values) != 3 or not all([key >= 0 for key in dict_values.keys()]):
        raise ValueError("The energy model requires 3 keys corresponding to positive "
                         "number of electrons! Given keys={0}".format(dict_values.keys()))
    # find reference number of electrons
    n_ref = sorted(dict_values.keys())[1]
    if n_ref < 1:
        raise ValueError("The n_ref cannot be less than one! Given n_ref={0}".format(n_ref))
    # check that number of electrons differ by one
    if sorted(dict_values.keys()) != [n_ref - 1, n_ref, n_ref + 1]:
        raise ValueError("In current implementation, the number of electrons (keys) should "
                         "differ by one! Given keys={0}".format(dict_values.keys()))
    # check that all values have the same type
    if not all([isinstance(value, type(dict_values[n_ref])) for value in dict_values.values()]):
        raise ValueError("All values in dictionary should be of the same type!")
    # check size of array values are the same
    if hasattr(dict_values[n_ref], "__len__"):
        if not all([value.shape == dict_values[n_ref].shape for value in dict_values.values()]):
            raise ValueError("All array values in dictionary should have the same shape!")
    # get property values
    value_m, value_0, value_p = [dict_values[n] for n in sorted(dict_values.keys())]
    return n_ref, value_m, value_0, value_p


def check_number_electrons(n_elec, n_min, n_max):
    """Check number of electrons to be positive & print warning if outside of interpolation range.

    Parameters
    ----------
    n_elec : float
        Number of electrons.
    n_min : float
        Minimum number of electrons used for interpolation.
    n_max : float
        Maximum number of electrons used for interpolation.
    """
    if not isinstance(n_elec, (int, float)):
        raise ValueError("Number of electrons should be a single number. "
                         "Given n_elec={0}".format(n_elec))
    if n_elec < 0.0:
        raise ValueError("Number of electrons cannot be negative! n_elec={0}".format(n_elec))
    if not n_min <= n_elec <= n_max:
        log.warn("Property evaluated for n_elec={0} outside of interpolation "
                 "region [{1}, {2}].".format(n_elec, n_min, n_max))
