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
"""Test chemtools.toolbox.utils."""


import numpy as np

from numpy.testing import assert_raises

from chemtools import context, CubeGen
from chemtools.utils.wrappers import HortonMolecule
from chemtools.toolbox.utils import get_matching_attr, get_molecular_grid
from chemtools.toolbox.utils import get_dict_energy, get_dict_density, get_dict_population


def test_get_matching_attr_raises():
    # check molecule
    filename = [context.get_fn("test/ch4_uhf_ccpvdz.fchk"),
                context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")]
    assert_raises(ValueError, get_matching_attr, filename, "numbers")
    assert_raises(ValueError, get_matching_attr, filename, "coordinates")
    # check matching attribute
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_matching_attr, molecule, "numbers")
    assert_raises(ValueError, get_matching_attr, molecule, "coordinates")


def test_get_molecular_grid_raises():
    # check atomic numbers shape
    molecule = HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    grid = CubeGen.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))
    assert_raises(ValueError, get_molecular_grid, molecule, grid)
    assert_raises(ValueError, get_molecular_grid, [molecule], grid)
    assert_raises(ValueError, get_molecular_grid, [molecule, molecule], grid)
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_molecular_grid, molecule, grid)
    # check atomic numbers (to be added)

    # check atomic coordinate
    molecule = HortonMolecule.from_file(context.get_fn("test/water_b3lyp_sto3g.fchk"))
    grid = CubeGen.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))
    assert_raises(ValueError, get_molecular_grid, molecule, grid)
    assert_raises(ValueError, get_molecular_grid, molecule, grid)
    assert_raises(ValueError, get_molecular_grid, [molecule], grid)
    assert_raises(ValueError, get_molecular_grid, [molecule, molecule], grid)
    molecule = [HortonMolecule.from_file(context.get_fn("test/water_b3lyp_sto3g.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_molecular_grid, molecule, grid)
    molecule = [HortonMolecule.from_file(context.get_fn("test/water_b3lyp_sto3g.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+1_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_molecular_grid, molecule, grid)


def test_get_dict_energy_raises():
    # check molecule
    filename = [context.get_fn("test/ch4_uhf_ccpvdz.fchk"),
                context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")]
    assert_raises(ValueError, get_dict_energy, filename)
    assert_raises(ValueError, get_dict_energy, context.get_fn("test/ch4_uhf_ccpvdz.fchk"))
    assert_raises(ValueError, get_dict_energy, "gibberish")
    # check repeated molecules
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))]
    assert_raises(ValueError, get_dict_energy, molecule)
    # check molecules with the same number of molecules
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_dict_energy, molecule)


def test_get_dict_density_raises():
    # check molecule
    filename = [context.get_fn("test/ch4_uhf_ccpvdz.fchk"),
                context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")]
    assert_raises(ValueError, get_dict_density, filename, np.array([[0., 0., 0.]]))
    assert_raises(ValueError, get_dict_density, "gibberish", np.array([[0., 0., 0.]]))
    # check repeated molecules
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk"))]
    assert_raises(ValueError, get_dict_density, molecule, np.array([[0., 0., 0.]]))
    # check molecules with the same number of molecules
    molecule = [HortonMolecule.from_file(context.get_fn("test/ch4_uhf_ccpvdz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_dict_density, molecule, np.array([[0., 0., 0.]]))


def test_get_dict_population_raises():
    # check molecule
    assert_raises(ValueError, get_dict_population, "gibberish", "RMF", "hi")
    # check number of molecules
    molecule = [HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+1_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_dict_population, molecule, "RMF", "esp")
    # check condensing approach
    molecule = HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk"))
    assert_raises(ValueError, get_dict_population, molecule, "fm", "h")
    assert_raises(ValueError, get_dict_population, molecule, "rm", "hi")
    assert_raises(ValueError, get_dict_population, molecule, "gibberish", "esp")
    # check scheme
    molecule = [HortonMolecule.from_file(context.get_fn("test/h2o_q+0_ub3lyp_ccpvtz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q+1_ub3lyp_ccpvtz.fchk")),
                HortonMolecule.from_file(context.get_fn("test/h2o_q-1_ub3lyp_ccpvtz.fchk"))]
    assert_raises(ValueError, get_dict_population, molecule, "rmf", "gibberish")
