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


import sys
import numpy as np

from numpy.testing import assert_allclose

from chemtools.toolbox.densbased import DensityLocalTool

if sys.version_info.major == 2:
    from chemtools.wrappers2.molecule import Molecule

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_densbased_from_file_h2o():
    # test against multiwfn
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        tool = DensityLocalTool.from_file(fname, data['coords'], spin='ab', index=None)
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.gradient_norm, data['nuc_grad_norm'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.laplacian, data['nuc_lap'], rtol=1.e-7, atol=0.)
    ked_tf = 0.3 * (3.0 * np.pi**2.0)**(2.0 / 3.0) * data['nuc_dens'] ** (5. / 3.)
    assert_allclose(tool.ked_thomas_fermi, ked_tf, rtol=1.e-7, atol=0.)


def test_densbased_from_molecule_h2o():
    # test against multiwfn
    with path('chemtools.data', 'data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz') as fname:
        data = np.load(str(fname))
    with path('chemtools.data', 'h2o_q+0_ub3lyp_ccpvtz.fchk') as fname:
        mol = Molecule.from_file(fname)
        tool = DensityLocalTool.from_molecule(mol, data['coords'])
    assert_allclose(tool.density, data['nuc_dens'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.gradient, data['nuc_grad'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.gradient_norm, data['nuc_grad_norm'], rtol=1.e-7, atol=0.)
    assert_allclose(tool.laplacian, data['nuc_lap'], rtol=1.e-7, atol=0.)
    ked_tf = 0.3 * (3.0 * np.pi**2.0)**(2.0 / 3.0) * data['nuc_dens'] ** (5. / 3.)
    assert_allclose(tool.ked_thomas_fermi, ked_tf, rtol=1.e-7, atol=0.)
