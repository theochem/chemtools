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
# pylint: skip-file
""" Tests chemtools.utils.output
"""

from nose.tools import assert_raises
from chemtools.utils import output


def test_vmd_script_start():
    """ Tests output._vmd_script_start
    """
    assert output._vmd_script_start() == ('#!/usr/local/bin/vmd\n'
                                          '# VMD script written by save_state $Revision: 1.41 $\n'
                                          '# VMD version: 1.8.6\n'
                                          'set viewplist\n'
                                          'set fixedlist\n'
                                          '#\n'
                                          '# Display settings\n'
                                          'display projection Orthographic\n'
                                          'display nearclip set 0.000000\n'
                                          '#\n')

def test_vmd_script_molecule():
    """ Tests output._vmd_script_molecule
    """
    assert_raises(ValueError, output._vmd_script_molecule)
    assert_raises(TypeError, output._vmd_script_molecule, 'example.log')
    assert ('# load new molecule\n'
            'mol new test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'
            '\n#\n'
            '# representation of the atoms\n'
            'mol delrep 0 top\n'
            'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
            'mol color Name\n'
            'mol selection {{all}}\n'
            'mol material Opaque\n'
            'mol addrep top\n'
            '#\n') == output._vmd_script_molecule('test.xyz')
    assert ('# load new molecule\n'
            'mol new test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'
            '\n'
            'mol addfile test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor'
            ' all\n#\n'
            '# representation of the atoms\n'
            'mol delrep 0 top\n'
            'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
            'mol color Name\n'
            'mol selection {{all}}\n'
            'mol material Opaque\n'
            'mol addrep top\n'
            '#\n') == output._vmd_script_molecule('test.xyz', 'test.xyz')
    assert ('# load new molecule\n'
            'mol new test.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'
            '\n'
            'mol addfile test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor'
            ' all\n#\n'
            '# representation of the atoms\n'
            'mol delrep 0 top\n'
            'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
            'mol color Name\n'
            'mol selection {{all}}\n'
            'mol material Opaque\n'
            'mol addrep top\n'
            '#\n') == output._vmd_script_molecule('test.cube', 'test.xyz')
