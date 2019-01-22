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
"""Test chemtools.analysis.elf."""


from chemtools.toolbox.elf import ELF
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_h2o_b3lyp_sto3g_elf():
    with path('chemtools.data', 'water_b3lyp_sto3g.fchk') as file_path:
        elf = ELF.from_file(file_path)
    elf.dump_isosurface_files('h2o', isosurf=0.8)
