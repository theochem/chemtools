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
"""Test chemtools.utils.output."""

import shutil
import tempfile
from contextlib import contextmanager
from nose.tools import assert_raises
import numpy as np
from chemtools.utils import output


@contextmanager
def tmpdir(name):
    """Create temporary directory that gets deleted after accessing it."""
    dn = tempfile.mkdtemp(name)
    try:
        yield dn
    finally:
        shutil.rmtree(dn)

header = ('#!/usr/local/bin/vmd\n'
          '# VMD script written by save_state $Revision: 1.41 $\n'
          '# VMD version: 1.8.6\n'
          'set viewplist\n'
          'set fixedlist\n'
          '#\n'
          '# Display settings\n'
          'display projection Orthographic\n'
          'display nearclip set 0.000000\n'
          'color Name {C} gray\n'
          '#\n')


def test_vmd_script_start():
    """Test output._vmd_script_start."""
    assert output._vmd_script_start() == header


def test_vmd_script_molecule():
    """Test output._vmd_script_molecule."""
    assert_raises(ValueError, output._vmd_script_molecule)
    assert_raises(TypeError, output._vmd_script_molecule, 'example.log')
    assert output._vmd_script_molecule('test.xyz') == \
        ('# load new molecule\n'
         'mol new test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
         '#\n'
         '# representation of the atoms\n'
         'mol delrep 0 top\n'
         'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
         'mol color Name\n'
         'mol selection {{all}}\n'
         'mol material Opaque\n'
         'mol addrep top\n'
         '#\n')
    assert output._vmd_script_molecule('test.xyz', 'test.xyz') == \
        ('# load new molecule\n'
         'mol new test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
         'mol addfile test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor '
         'all\n'
         '#\n'
         '# representation of the atoms\n'
         'mol delrep 0 top\n'
         'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
         'mol color Name\n'
         'mol selection {{all}}\n'
         'mol material Opaque\n'
         'mol addrep top\n'
         '#\n')
    assert output._vmd_script_molecule('test.cube', 'test.xyz') == \
        ('# load new molecule\n'
         'mol new test.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n'
         'mol addfile test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor '
         'all\n'
         '#\n'
         '# representation of the atoms\n'
         'mol delrep 0 top\n'
         'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
         'mol color Name\n'
         'mol selection {{all}}\n'
         'mol material Opaque\n'
         'mol addrep top\n'
         '#\n')


def test_vmd_script_isosurface():
    """Test output._vmd_script_isosurface."""
    assert_raises(TypeError, output._vmd_script_isosurface, isosurf=None)
    assert_raises(TypeError, output._vmd_script_isosurface, isosurf=1)
    assert_raises(TypeError, output._vmd_script_isosurface, index=1.0)
    assert_raises(TypeError, output._vmd_script_isosurface, index=None)
    assert_raises(TypeError, output._vmd_script_isosurface, show_type='boxes')
    assert_raises(TypeError, output._vmd_script_isosurface, show_type=None)
    assert_raises(TypeError, output._vmd_script_isosurface, draw_type='lskdfj')
    assert_raises(TypeError, output._vmd_script_isosurface, draw_type=None)
    assert_raises(TypeError, output._vmd_script_isosurface, material='lksjdf')
    assert_raises(TypeError, output._vmd_script_isosurface, material=None)
    assert_raises(TypeError, output._vmd_script_isosurface, scalemin=1)
    assert_raises(TypeError, output._vmd_script_isosurface, scalemin=None)
    assert_raises(TypeError, output._vmd_script_isosurface, scalemax=1)
    assert_raises(TypeError, output._vmd_script_isosurface, scalemax=None)
    assert_raises(TypeError, output._vmd_script_isosurface, colorscheme=-1)
    assert_raises(TypeError, output._vmd_script_isosurface, colorscheme=1057)
    assert_raises(TypeError, output._vmd_script_isosurface, colorscheme='asdfasdf')
    assert_raises(TypeError, output._vmd_script_isosurface, colorscheme=None)

    assert output._vmd_script_isosurface() == ('# add representation of the surface\n'
                                               'mol representation Isosurface 0.50000 0 0 0 1 1\n'
                                               'mol color Volume 0\n'
                                               'mol selection {all}\n'
                                               'mol material Opaque\n'
                                               'mol addrep top\n'
                                               'mol selupdate 1 top 0\n'
                                               'mol colupdate 1 top 0\n'
                                               'mol scaleminmax top 1 -0.050000 0.050000\n'
                                               'mol smoothrep top 1 0\n'
                                               'mol drawframes top 1 {now}\n'
                                               'color scale method RGB\n'
                                               '#\n')
    assert output._vmd_script_isosurface(colorscheme=1) == \
        ('# add representation of the surface\n'
         'mol representation Isosurface 0.50000 0 0 0 1 1\n'
         'mol color ColorID 1\n'
         'mol selection {all}\n'
         'mol material Opaque\n'
         'mol addrep top\n'
         'mol selupdate 1 top 0\n'
         'mol colupdate 1 top 0\n'
         'mol scaleminmax top 1 -0.050000 0.050000\n'
         'mol smoothrep top 1 0\n'
         'mol drawframes top 1 {now}\n'
         'color scale method RGB\n'
         '#\n')


def test_vmd_script_vector_field():
    """Test output._vmd_script_vector_field."""
    centers = np.array([[1, 2, 3]])
    unit_vecs = np.array([[1, 0, 0]])
    weights = np.array([1])
    assert_raises(ValueError, output._vmd_script_vector_field, centers, unit_vecs, np.array([1, 2]))
    assert_raises(ValueError, output._vmd_script_vector_field, centers, np.array([[1, 2, 3]]),
                  weights)

    assert_raises(TypeError, output._vmd_script_vector_field, np.array([1, 2, 3]),
                  unit_vecs, weights)
    assert_raises(TypeError, output._vmd_script_vector_field, np.array([[1, 2, 3, 4]]),
                  unit_vecs, weights)
    assert_raises(TypeError, output._vmd_script_vector_field, [[1, 2, 3]],
                  unit_vecs, weights)

    assert_raises(TypeError, output._vmd_script_vector_field, centers, np.array([1, 2, 3]), weights)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, np.array([[1, 2, 3, 4]]),
                  weights)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, [[1, 2, 3]],
                  weights)

    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, np.array([[1]]))
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, [1])

    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  has_shadow=None)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  has_shadow=0)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  material='lksjdf')
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  material=None)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  color=-1)
    assert_raises(TypeError, output._vmd_script_vector_field, centers, unit_vecs, weights,
                  color=1057)
    assert output._vmd_script_vector_field(centers, unit_vecs, weights) == \
        ('# Add function for vector field\n'
         'proc vmd_draw_arrow {mol center unit_dir cyl_radius cone_radius length} {\n'
         'set start [vecsub $center [vecscale [vecscale 0.5 $length] $unit_dir]]\n'
         'set end [vecadd $start [vecscale $length $unit_dir]]\n'
         'set middle [vecsub $end [vecscale [vecscale 1.732050808 $cone_radius] $unit_dir]]\n'
         'graphics $mol cylinder $start $middle radius $cyl_radius\n'
         'graphics $mol cone $middle $end radius $cone_radius\n'
         '}\n'
         '#\n'
         'draw materials on\n'
         'draw material Transparent\n'
         'draw color 0\n'
         'draw arrow {1 2 3} {1 0 0} 0.08 0.15 0.7\n'
         '#\n')
    assert output._vmd_script_vector_field(centers, unit_vecs, weights, has_shadow=False) == \
        ('# Add function for vector field\n'
         'proc vmd_draw_arrow {mol center unit_dir cyl_radius cone_radius length} {\n'
         'set start [vecsub $center [vecscale [vecscale 0.5 $length] $unit_dir]]\n'
         'set end [vecadd $start [vecscale $length $unit_dir]]\n'
         'set middle [vecsub $end [vecscale [vecscale 1.732050808 $cone_radius] $unit_dir]]\n'
         'graphics $mol cylinder $start $middle radius $cyl_radius\n'
         'graphics $mol cone $middle $end radius $cone_radius\n'
         '}\n'
         '#\n'
         'draw materials off\n'
         'draw material Transparent\n'
         'draw color 0\n'
         'draw arrow {1 2 3} {1 0 0} 0.08 0.15 0.7\n'
         '#\n')
    assert output._vmd_script_vector_field(centers, unit_vecs, np.array([1e-2])) == \
        ('# Add function for vector field\n'
         'proc vmd_draw_arrow {mol center unit_dir cyl_radius cone_radius length} {\n'
         'set start [vecsub $center [vecscale [vecscale 0.5 $length] $unit_dir]]\n'
         'set end [vecadd $start [vecscale $length $unit_dir]]\n'
         'set middle [vecsub $end [vecscale [vecscale 1.732050808 $cone_radius] $unit_dir]]\n'
         'graphics $mol cylinder $start $middle radius $cyl_radius\n'
         'graphics $mol cone $middle $end radius $cone_radius\n'
         '}\n'
         '#\n'
         'draw materials on\n'
         'draw material Transparent\n'
         'draw color 0\n'
         '#\n')


def test_print_vmd_script_isosurface():
    """Test print_vmd_script_isosurface."""
    # check TypeError:
    assert_raises(TypeError, output.print_vmd_script_isosurface, 'test.vmd', 'iso.cube',
                  colorscheme=[1], negative=True)

    with tmpdir('chemtools.utils.test.test_base.test_vmd_script_isosurface') as dn:
        vmd = '%s/%s' % (dn, 'test.vmd')

        output.print_vmd_script_isosurface(vmd, 'iso.cube')

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new iso.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 '
                 'waitfor all\n'
                 '#\n'
                 '# representation of the atoms\n'
                 'mol delrep 0 top\n'
                 'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
                 'mol color Name\n'
                 'mol selection {{all}}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 '#\n'
                 '# add representation of the surface\n'
                 'mol representation Isosurface 0.50000 0 0 0 1 1\n'
                 'mol color ColorID 0\n'
                 'mol selection {all}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 'mol selupdate 1 top 0\n'
                 'mol colupdate 1 top 0\n'
                 'mol scaleminmax top 1 -0.050000 0.050000\n'
                 'mol smoothrep top 1 0\n'
                 'mol drawframes top 1 {now}\n'
                 'color scale method RGB\n'
                 '#\n')

        output.print_vmd_script_isosurface(vmd, 'iso.cube', colorfile='col.cube',
                                           isosurf=0.6, material='Transparent',
                                           scalemin=-0.06, scalemax=0.08)

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new col.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 '
                 'waitfor all\n'
                 'mol addfile iso.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 '
                 'waitfor all\n'
                 '#\n'
                 '# representation of the atoms\n'
                 'mol delrep 0 top\n'
                 'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
                 'mol color Name\n'
                 'mol selection {{all}}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 '#\n'
                 '# add representation of the surface\n'
                 'mol representation Isosurface 0.60000 1 0 0 1 1\n'
                 'mol color Volume 0\n'
                 'mol selection {all}\n'
                 'mol material Transparent\n'
                 'mol addrep top\n'
                 'mol selupdate 1 top 0\n'
                 'mol colupdate 1 top 0\n'
                 'mol scaleminmax top 1 -0.060000 0.080000\n'
                 'mol smoothrep top 1 0\n'
                 'mol drawframes top 1 {now}\n'
                 'color scale method RGB\n'
                 '#\n')

        output.print_vmd_script_isosurface(vmd, 'iso.cube', colorscheme=[0, 1], negative=True)

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new iso.cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 '
                 'waitfor all\n'
                 '#\n'
                 '# representation of the atoms\n'
                 'mol delrep 0 top\n'
                 'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
                 'mol color Name\n'
                 'mol selection {{all}}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 '#\n'
                 '# add representation of the surface\n'
                 'mol representation Isosurface 0.50000 0 0 0 1 1\n'
                 'mol color ColorID 0\n'
                 'mol selection {all}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 'mol selupdate 1 top 0\n'
                 'mol colupdate 1 top 0\n'
                 'mol scaleminmax top 1 -0.050000 0.050000\n'
                 'mol smoothrep top 1 0\n'
                 'mol drawframes top 1 {now}\n'
                 'color scale method RGB\n'
                 '#\n'
                 '# add representation of the surface\n'
                 'mol representation Isosurface -0.50000 0 0 0 1 1\n'
                 'mol color ColorID 1\n'
                 'mol selection {all}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 'mol selupdate 1 top 0\n'
                 'mol colupdate 1 top 0\n'
                 'mol scaleminmax top 1 -0.050000 0.050000\n'
                 'mol smoothrep top 1 0\n'
                 'mol drawframes top 1 {now}\n'
                 'color scale method RGB\n'
                 '#\n')


def test_print_vmd_script_multiple_cube():
    """Test print_vmd_script_multiple_cube."""
    # check TypeError and ValueError:
    assert_raises(TypeError, output.print_vmd_script_multiple_cube, 'test.vmd', 'iso.cube')
    assert_raises(ValueError, output.print_vmd_script_multiple_cube, 'test.vmd',
                  ['iso.cube', 'iso.wrong_end'])

    ratom = ('# representation of the atoms\n'
             'mol delrep 0 top\n'
             'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
             'mol color Name\n'
             'mol selection {{all}}\n'
             'mol material Opaque\n'
             'mol addrep top\n'
             '#\n')

    def rsurf(iso, n, c):
        """Generate representation of the surface."""
        return ('# add representation of the surface\n'
                'mol representation Isosurface {0} {1} 0 0 1 1\n'.format(iso, n) +
                'mol color ColorID {0}\n'.format(c) +
                'mol selection {all}\n'
                'mol material Opaque\n'
                'mol addrep top\n'
                'mol selupdate 1 top 0\n'
                'mol colupdate 1 top 0\n'
                'mol scaleminmax top 1 -0.050000 0.050000\n'
                'mol smoothrep top 1 0\n'
                'mol drawframes top 1 {now}\n'
                'color scale method RGB\n'
                '#\n')

    with tmpdir('chemtools.utils.test.test_base.test_vmd_script_multiple_cube') as dn:
        vmd = '%s/%s' % (dn, 'test.vmd')
        c1 = '%s/%s' % (dn, 'iso1.cube')
        c2 = '%s/%s' % (dn, 'iso2.cube')
        open(c1, 'a').close()
        open(c2, 'a').close()

        output.print_vmd_script_multiple_cube(vmd, [c1, c2])

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c1) +
                 'autobonds 1 waitfor all\n'
                 'mol addfile {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c2) +
                 'autobonds 1 waitfor all\n'
                 '#\n' + ratom + rsurf('0.50000', '0', '0') + rsurf('0.50000', '1', '1'))

        output.print_vmd_script_multiple_cube(vmd, [c1, c2], isosurfs=0.6)

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c1) +
                 'autobonds 1 waitfor all\n'
                 'mol addfile {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c2) +
                 'autobonds 1 waitfor all\n'
                 '#\n' + ratom + rsurf('0.60000', '0', '0') + rsurf('0.60000', '1', '1'))

        output.print_vmd_script_multiple_cube(vmd, [c1, c2], isosurfs=[0.6, 0.8], colors=[3, 4])

        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c1) +
                 'autobonds 1 waitfor all\n'
                 'mol addfile {0} type cube first 0 last -1 step 1 filebonds 1 '.format(c2) +
                 'autobonds 1 waitfor all\n'
                 '#\n' + ratom + rsurf('0.60000', '0', '3') + rsurf('0.80000', '1', '4'))

        # check TypeError and ValueError:
        assert_raises(TypeError, output.print_vmd_script_multiple_cube, vmd, [c1, c2],
                      isosurfs=[0.6, 0.8, 0.4])
        assert_raises(TypeError, output.print_vmd_script_multiple_cube, vmd, [c1, c2],
                      isosurfs=[0.6, 'error'])
        assert_raises(TypeError, output.print_vmd_script_multiple_cube, vmd, [c1, c2],
                      isosurfs=[0.6, 0.8], colors=[3, 4, 8])
        assert_raises(ValueError, output.print_vmd_script_multiple_cube, vmd, [c1, c2],
                      isosurfs=[0.6, 0.8], colors=[3, 1060])


def test_print_vmd_script_vector_field():
    """Test output.print_vmd_script_vector_field."""
    centers = np.array([[1, 2, 3]])
    vecs = np.array([[1, 0, 0]])
    with tmpdir('chemtools.utils.test.test_base.test_print_vmd_script_vector_field') as dn:
        vmd = '%s/%s' % (dn, 'test.vmd')
        output.print_vmd_script_vector_field(vmd, 'test.xyz', centers, vecs)
        with open(vmd, 'r') as content_file:
            assert content_file.read() == \
                (header +
                 '# load new molecule\n'
                 'mol new test.xyz type {xyz} first 0 last -1 step 1 filebonds 1 autobonds 1 '
                 'waitfor all\n'
                 '#\n'
                 '# representation of the atoms\n'
                 'mol delrep 0 top\n'
                 'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
                 'mol color Name\n'
                 'mol selection {{all}}\n'
                 'mol material Opaque\n'
                 'mol addrep top\n'
                 '#\n'
                 '# Add function for vector field\n'
                 'proc vmd_draw_arrow {mol center unit_dir cyl_radius cone_radius length} {\n'
                 'set start [vecsub $center [vecscale [vecscale 0.5 $length] $unit_dir]]\n'
                 'set end [vecadd $start [vecscale $length $unit_dir]]\n'
                 'set middle [vecsub $end [vecscale [vecscale 1.732050808 $cone_radius] '
                 '$unit_dir]]\n'
                 'graphics $mol cylinder $start $middle radius $cyl_radius\n'
                 'graphics $mol cone $middle $end radius $cone_radius\n'
                 '}\n'
                 '#\n'
                 'draw materials on\n'
                 'draw material Transparent\n'
                 'draw color 0\n'
                 'draw arrow {1 2 3} {1.0 0.0 0.0} 0.08 0.15 0.7\n'
                 '#\n')
