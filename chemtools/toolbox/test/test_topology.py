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
"""Test chemtools.toolbox.topology."""


import numpy as np

from chemtools.toolbox.topology import TopologicalTool
from chemtools.wrappers.molecule import Molecule
from chemtools.utils.cube import UniformGrid

try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


def test_critical_point_h2o():
    # test against multiwfn 3.6 dev src
    with path("chemtools.data", "data_multiwfn36_fchk_h2o_q+0_ub3lyp_ccpvtz.npz") as fname:
        data = np.load(str(fname))
    nna, bcp = data["nna_coords"], data["bcp_coords"]
    # find critical points
    with path("chemtools.data", "h2o_q+0_ub3lyp_ccpvtz.fchk") as fpath:
        mol = Molecule.from_file(fpath)
    cub = UniformGrid.from_molecule(mol, spacing=0.15, extension=0.1, rotate=False)
    top = TopologicalTool.from_molecule(mol, points=cub.points)
    # check NA
    assert len(top.nna) == 3
    assert sum([np.allclose(top.nna[0].point, point, rtol=0.0) for point in nna]) == 1
    assert sum([np.allclose(top.nna[1].point, point, rtol=0.0) for point in nna]) == 1
    assert sum([np.allclose(top.nna[2].point, point, rtol=0.0) for point in nna]) == 1
    # check BCP
    assert len(top.bcp) == 2
    print(sum([np.allclose(top.bcp[1].point, point, rtol=0.0) for point in bcp]))
    assert sum([np.allclose(top.bcp[0].point, point, rtol=0.0) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[1].point, point, rtol=0.0) for point in bcp]) == 1
    # check total number of CP
    assert len(top.rcp) == 0
    assert len(top.ccp) == 0
    assert len(top.cps) == 5
    assert top.poincare_hopf_equation


def test_critical_point_li2():
    # test against multiwfn 3.6 dev src
    with path("chemtools.data", "data_multiwfn36_fchk_li2_q+0_ub3lyp_ccpvdz.npz") as fname:
        data = np.load(str(fname))
        nna, bcp = data["nna_coords"], data["bcp_coords"]
    # find critical points
    with path("chemtools.data", "li2_q+0_ub3lyp_ccpvdz.fchk") as fpath:
        mol = Molecule.from_file(fpath)
    cub = UniformGrid.from_molecule(mol, spacing=0.25, extension=0.1, rotate=False)
    top = TopologicalTool.from_molecule(mol, points=cub.points)
    # check NA
    assert len(top.nna) == 3
    assert sum([np.allclose(top.nna[0].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[1].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[2].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    # check BCP
    assert len(top.bcp) == 2
    assert sum([np.allclose(top.bcp[0].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[1].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    # check total number of CP
    assert len(top.rcp) == 0
    assert len(top.ccp) == 0
    assert len(top.cps) == 5
    assert top.poincare_hopf_equation


def test_critical_point_cyclopropenium_oxide():
    # test against multiwfn 3.6 dev src
    with path("chemtools.data", "data_multiwfn36_fchk_c3h2o_q+0_ub3lyp_augccpvdz.npz") as fname:
        data = np.load(str(fname))
        nna, bcp, rcp = data["nna_coords"], data["bcp_coords"], data["rcp_coords"]
    # find critical points
    with path("chemtools.data", "c3h2o_q+0_ub3lyp_augccpvdz.fchk") as fpath:
        mol = Molecule.from_file(fpath)
    cub = UniformGrid.from_molecule(mol, spacing=0.25, extension=0.1, rotate=False)
    top = TopologicalTool.from_molecule(mol, points=cub.points)
    # check NA
    assert len(top.nna) == 6
    assert sum([np.allclose(top.nna[0].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[1].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[2].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[3].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[4].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[5].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    # check BCP
    assert len(top.bcp) == 6
    assert sum([np.allclose(top.bcp[0].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[1].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[2].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[3].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[4].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[5].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    # check RCP
    assert len(top.rcp) == 1
    assert sum([np.allclose(top.rcp[0].point, point, rtol=0.0, atol=1.0e-5) for point in rcp]) == 1
    # check total number of CP
    assert len(top.ccp) == 0
    assert len(top.cps) == 13
    assert top.poincare_hopf_equation


def test_critical_point_c4h4():
    # test against multiwfn 3.6 dev src
    with path("chemtools.data", "data_multiwfn36_fchk_c4h4_ub3lyp_ccpvdz.npz") as fname:
        data = np.load(str(fname))
        nna, bcp = data["nna_coords"], data["bcp_coords"]
        rcp, ccp = data["rcp_coords"], data["ccp_coords"]
    # find critical points
    with path("chemtools.data", "c4h4_ub3lyp_ccpvdz.fchk") as fpath:
        mol = Molecule.from_file(fpath)
    cub = UniformGrid.from_molecule(mol, spacing=0.25, extension=0.1, rotate=False)
    top = TopologicalTool.from_molecule(mol, points=cub.points)
    # check NA
    assert len(top.nna) == 8
    assert sum([np.allclose(top.nna[0].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[1].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[2].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[3].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[4].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[5].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[6].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    assert sum([np.allclose(top.nna[7].point, point, rtol=0.0, atol=1.0e-5) for point in nna]) == 1
    # check BCP
    assert len(top.bcp) == 10
    assert sum([np.allclose(top.bcp[0].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[2].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[3].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[4].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[5].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[6].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[7].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[8].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    assert sum([np.allclose(top.bcp[9].point, point, rtol=0.0, atol=1.0e-5) for point in bcp]) == 1
    # check RCP
    assert len(top.rcp) == 4
    assert sum([np.allclose(top.rcp[0].point, point, rtol=0.0, atol=1.0e-5) for point in rcp]) == 1
    assert sum([np.allclose(top.rcp[1].point, point, rtol=0.0, atol=1.0e-5) for point in rcp]) == 1
    assert sum([np.allclose(top.rcp[2].point, point, rtol=0.0, atol=1.0e-5) for point in rcp]) == 1
    assert sum([np.allclose(top.rcp[3].point, point, rtol=0.0, atol=1.0e-5) for point in rcp]) == 1
    # check CCP
    assert len(top.ccp) == 1
    assert sum([np.allclose(top.ccp[0].point, point, rtol=0.0, atol=1.0e-5) for point in ccp]) == 1
    # check total number of CP
    assert len(top.cps) == 23
    assert top.poincare_hopf_equation
