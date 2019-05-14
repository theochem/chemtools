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
# pragma pylint: disable=invalid-name,bad-whitespace,bad-continuation
"""Test chemtools.toolbox.densbased."""


import numpy as np

from numpy.testing import assert_array_almost_equal

from chemtools.toolbox.densbased import DensityLocalTool
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_densbased_from_file_elf_h2o_dimer():
    # load data computed with NCIPLOT by E.R. Johnson and J. Contreras-Garcia
    with path("chemtools.data", "data_elf_nciplot_h2o_dimer_pbe_sto3g.npz") as filename:
        data = np.load(str(filename))
    # test from_file initialization & check ELF
    with path("chemtools.data", "h2o_dimer_pbe_sto3g.fchk") as filename:
        tool = DensityLocalTool.from_file(filename, spin='ab', index=None, points=data["points"])
    # assert_array_almost_equal(tool.electron_localization_function, data["elf"], decimal=5)
