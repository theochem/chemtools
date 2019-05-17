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
"""Test chemtools.toolbox.interactions."""


import numpy as np
from numpy.testing import assert_allclose, assert_raises
from chemtools.toolbox.interactions import ELF, LOL
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


# def test_h2o_b3lyp_sto3g_elf():
#     with path('chemtools.data', 'water_b3lyp_sto3g.fchk') as file_path:
#         elf = ELF.from_file(file_path)
#     elf.generate_scripts('h2o', isosurf=0.8)


def test_elf_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    elf = ELF(data['nuc_dens'], data['nuc_grad'], data['nuc_ked_pd'])
    assert_allclose(elf.value, data['nuc_elf'], rtol=1.e-6, atol=1.e-6)
    # check raises
    dens, grad, ked = data['nuc_dens'], data['nuc_grad'], data['nuc_ked_pd']
    assert_raises(ValueError, ELF, dens, grad, ked, trans_k=-1)
    assert_raises(ValueError, ELF, dens, grad, ked, trans_a=2, trans_k=-1)
    assert_raises(ValueError, ELF, dens, grad, ked, trans_k=0)
    assert_raises(ValueError, ELF, dens, grad, ked, trans_a=0)
    assert_raises(ValueError, ELF, dens, grad, ked, trans='inverse_rational')


def test_lol_h2o_nuclei():
    # test against multiwfn 3.6 dev src
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    lol = LOL(data['nuc_dens'], data['nuc_grad'], data['nuc_ked_pd'])
    assert_allclose(lol.value, data['nuc_lol'], rtol=1.e-6, atol=1.e-6)
    # check raises
    dens, grad, ked = data['nuc_dens'], data['nuc_grad'], data['nuc_ked_pd']
    assert_raises(ValueError, LOL, dens, grad, ked, trans_k=-1)
    assert_raises(ValueError, LOL, dens, grad, ked, trans_a=2, trans_k=-1)
    assert_raises(ValueError, LOL, dens, grad, ked, trans_k=0)
    assert_raises(ValueError, LOL, dens, grad, ked, trans_a=0)
    assert_raises(ValueError, LOL, dens, grad, ked, trans='rational')
