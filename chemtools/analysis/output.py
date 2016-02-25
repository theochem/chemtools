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
'''Output methods for Analyze Module.'''

def _print_vmd(file, densfile, rdgfile, isosurf, denscut):
    # print out the .vmd file as created with NCIPLOT by E.R. Johnson and J. Contreras-Garcia.
    with open(file, 'w') as f:

        print >> f, '#!/usr/local/bin/vmd'
        print >> f, '# VMD script written by save_state $Revision: 1.41 $'
        print >> f, '# VMD version: 1.8.6               '
        print >> f, 'set viewplist                      '
        print >> f, 'set fixedlist                      '
        print >> f, '#'
        print >> f, '# Display settings                 '
        print >> f, 'display projection   Orthographic  '
        print >> f, 'display nearclip set 0.000000      '
        print >> f, '#'
        print >> f, '# load new molecule                '
        print >> f, 'mol new {0}  type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'.format(densfile)
        print >> f, 'mol addfile {0}  type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'.format(rdgfile)
        print >> f, '#'
        print >> f, '# representation of the atoms'
        print >> f, 'mol delrep 0 top'
        print >> f, 'mol representation CPK 1.000000 0.300000 118.000000 131.000000'
        print >> f, 'mol color Name'
        print >> f, 'mol selection {all}'
        print >> f, 'mol material Opaque'
        print >> f, 'mol addrep top'
        print >> f, '#'
        print >> f, '# add representation of the surface'
        print >> f, 'mol representation Isosurface {0:.5f} 1.0 0.0 0.0 1 1'.format(isosurf)
        print >> f, 'mol color Volume 0'
        print >> f, 'mol selection {all}'
        print >> f, 'mol material Opaque'
        print >> f, 'mol addrep top'
        print >> f, 'mol selupdate 1 top 0'
        print >> f, 'mol colupdate 1 top 0'
        print >> f, 'mol scaleminmax top 1 {0:.6f} {1:.6f}'.format(-denscut, denscut)
        print >> f, 'mol smoothrep top 1 0'
        print >> f, 'mol drawframes top 1 {now}'
        print >> f, 'color scale method BGR'
        print >> f, 'set colorcmds {{color Name {C} gray}}'
        print >> f, '#some more'

