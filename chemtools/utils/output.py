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
"""Output Module.

   This module contains functions for generating scripts for visualizing
   purposes using VMD, GaussView, etc.
"""
import os
import numpy as np


def _vmd_script_start():
    """ Generates part of the beginning part of the VMD script
    """
    return ('#!/usr/local/bin/vmd\n'
            '# VMD script written by save_state $Revision: 1.41 $\n'
            '# VMD version: 1.8.6\n'
            'set viewplist\n'
            'set fixedlist\n'
            '#\n'
            '# Display settings\n'
            'display projection Orthographic\n'
            'display nearclip set 0.000000\n'
            '#\n'
           )


def _vmd_script_molecule(*mol_files):
    """ Generates part of the VMD script that loads the molecule information

    Parameters
    ----------
    mol_files : str
        Names of the input files that represent the moelcule
        .xyz and .cube files are supported
        Cube files can correspond to the density, reduced density gradient, isosurface, color of
        isosurface, etc

    Returns
    -------
    Part of the VMD script that constructs the molecule

    Raises
    ------
    TypeError
        If unsupported file type (i.e. not xyz or cube)
    """
    output = '# load new molecule\n'
    for i, mol in enumerate(mol_files):
        if i == 0:
            mol_type = 'new'
        else:
            mol_type = 'addfile'

        ext = os.path.splitext(mol)[1]
        if  ext == '.xyz':
            file_type = 'xyz'
        elif ext == '.cube':
            file_type = 'cube'
        else:
            raise TypeError('Unsupported file type, {0}'.format(ext))

        output += ('mol {0} {1} type {2} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'
                   '\n'.format(mol_type, mol, file_type))
    output += ('#\n'
               '# representation of the atoms\n'
               'mol delrep 0 top\n'
               'mol representation CPK 1.000000 0.300000 118.000000 131.000000\n'
               'mol color Name\n'
               'mol selection {{all}}\n'
               'mol material Opaque\n'
               'mol addrep top\n'
               '#\n'
              )
    return output


def _vmd_script_isosurface(isosurf=0.5, index=0, show_type='isosurface', draw_type='solid surface',
                           material='Opaque', scalemin=-0.05, scalemax=0.05, colorscheme='RGB'):
    """ Generates part of the VMD script that configures the isosurface

    Parameters
    ----------
    isosurf : float
        Isovalue at which the isosurface is plotted
        Default is 0.5
    index : int
        Index of the file that contains the isosurface data (in the order loaded by
        `_vmd_script_molecule`)
        Default is 0
    show_type : str
        Option that controls what will be shown
        One of 'isosurface', 'box', and 'box+isosurface'
    draw_type : str
        Option that controls how the isosurface will be drawn
        One of 'solid surface', 'wireframe', 'points', 'shaded points'
    material : str
        The material setting of the isosurface
        One of 'Opaque', 'Transparent', 'BrushedMetal', 'Diffuse', 'Ghost', 'Glass1', 'Glass2',
        'Glass3', 'Glossy', 'HardPlastic', 'MetallicPastel', 'Steel', 'Translucent', 'Edgy',
        'EdgyShiny', 'EdgyGlass', 'Goodsell', 'AOShiny', 'AOChalky', 'AOEdgy', 'BlownGlass',
        'GlassBubble', 'RTChrome'.
        Default is 'Opaque'
    scalemin : float
        Smallest value to color on the isosurface
        Default is -0.05
    scalemax : float
        Largest value to color on the isosurface
        Default is 0.05
    colorscheme : str, int
        Color scheme used for the isosurface.
        If str, then appropriate gradient is used to color the isosurface. It must be one of the
        following options
            =======  =====================================
            Options  Description
            =======  =====================================
            'RGB'    small=red, middle=green, large=blue
            'BGR'    small=blue, middle=green, large=red
            'RWB'    small=red, middle=white, large=blue
            'BWR'    small=blue, middle=white, large=red
            'RWG'    small=red, middle=white, large=green
            'GWR'    small=green, middle=white, large=red
            'GWB'    small=green, middle=white, large=blue
            'BWG'    small=blue, middle=white, large=green
            'BlkW'   small=black, large=white
            'WBlk'   small=white, large=black
            =======  =====================================
        If int, isosurface is colored with just one color. There are 32 colors are available in
        VMD, check the program or website for more details.
        Default is 'RGB'

    Returns
    -------
    Part of the VMD script that constructs the isosurface

    Raises
    ------
    TypeError
        If `isosurf` is not a float
        If `index` is not an integer
        If `show_type` is not one of 'isosurface', 'box', or 'box+isosurface'
        If `draw_type` is not one of 'solid surface', 'wireframe', 'points', or 'shaded points'
        If `material` is not one of 'Opaque', 'Transparent', 'BrushedMetal', 'Diffuse', 'Ghost',
        'Glass1', 'Glass2', 'Glass3', 'Glossy', 'HardPlastic', 'MetallicPastel', 'Steel',
        'Translucent', 'Edgy', 'EdgyShiny', 'EdgyGlass', 'Goodsell', 'AOShiny', 'AOChalky',
        'AOEdgy', 'BlownGlass', 'GlassBubble', 'RTChrome'
        If `scalemin` is not float
        If `scalemax` is not float
        If `colorscheme` is not one of 'RGB', 'BGR', 'RWB', 'BWR', 'RWG', 'GWR', 'GWB', 'BWG',
        'BlkW', 'WBlk'
    """
    if not isinstance(isosurf, float):
        raise TypeError('`isosurf` must be a float')

    if not isinstance(index, int):
        raise TypeError('`index` must be an integer')

    if show_type not in ['isosurface', 'box', 'box+isosurface']:
        raise TypeError('Unsupported `show_type`. Must be one of {0}, {1}, or {2}'
                        ''.format('box', 'isosurface', 'box+isosurface'))
    show_type = {'isosurface':0, 'box':1, 'box+isosurface':2}[show_type]

    if draw_type not in ['solid surface', 'wireframe', 'points', 'shaded points']:
        raise TypeError('Unsupported `draw_type`. Must be one of {0}, {1}, {2} or {3}'
                        ''.format('solid surface', 'wireframe', 'points', 'shaded points'))
    draw_type = {'solid surface':0, 'wireframe':1, 'points':2, 'shaded points':3}[draw_type]

    allowed_materials = ['Opaque', 'Transparent', 'BrushedMetal', 'Diffuse', 'Ghost', 'Glass1',
                         'Glass2', 'Glass3', 'Glossy', 'HardPlastic', 'MetallicPastel', 'Steel',
                         'Translucent', 'Edgy', 'EdgyShiny', 'EdgyGlass', 'Goodsell', 'AOShiny',
                         'AOChalky', 'AOEdgy', 'BlownGlass', 'GlassBubble', 'RTChrome']
    if material not in allowed_materials:
        raise TypeError('Unsupported `material`. Must be one of {0}'.format(allowed_materials))

    if not isinstance(scalemin, float):
        raise TypeError('`scalemin` must be a float')
    if not isinstance(scalemax, float):
        raise TypeError('`scalemax` must be a float')

    if (isinstance(colorscheme, str) and
            colorscheme not in ['RGB', 'BGR', 'RWB', 'BWR', 'RWG', 'GWR', 'GWB', 'BWG', 'BlkW',
                                'WBlk']):
        raise TypeError('Unsupported colorscheme, {0}'.format(colorscheme))
    if isinstance(colorscheme, int) and not (0 <= colorscheme < 33):
        raise TypeError('Unsupported colorscheme, {0}'.format(colorscheme))

    output = '# add representation of the surface\n'
    # mol representation Isosurface {0} {1} {2} {3} {4} {5}
    # {0} = isovalue
    # {1} = file index that will be plotted
    # {2} = what to show
    #       0 -> isosurface
    #       1 -> box
    #       2 -> box and isosurface
    # {3} = how to draw
    #       0 -> solid surface
    #       1 -> wireframe
    #       2 -> points
    #       3 -> shaded points
    # {4} = (integer) step size (distance between plot objects)
    # {5} = (integer) width of plot object (e.g. wire, point)
    output += ('mol representation Isosurface {isosurf:.5f} {index} {show} {draw} 1 1\n'
               ''.format(isosurf=isosurf, index=index, show=show_type, draw=draw_type))

    if isinstance(colorscheme, int):
        output += 'mol color ColorID {0}\n'.format(colorscheme)
    else:
        output += 'mol color Volume 0\n'

    output += ('mol selection {{all}}\n'
               'mol material {material}\n'
               'mol addrep top\n'
               'mol selupdate 1 top 0\n'
               'mol colupdate 1 top 0\n'
               'mol scaleminmax top 1 {scalemin:.6f} {scalemax:.6f}\n'
               'mol smoothrep top 1 0\n'
               'mol drawframes top 1 {{now}}\n'.format(material=material, scalemin=scalemin,
                                                       scalemax=scalemax)
              )

    if isinstance(colorscheme, str):
        output += 'color scale method {0}\n'.format(colorscheme)
    else:
        output += 'color scale method RGB\n'

    output += 'set colorcmds {{{{color Name {{C}} gray}}}}\n'
    output += '#\n'
    return output


def _vmd_script_vector_field(centers, unit_directions, weights):
    """ Generates part of the VMD script that constructs the vector field

    Parameters
    ----------
    centers : np.ndarray(N, 3)
        Coordinates of the centers of each vector
    unit_directions : np.ndarray(N, 3)
        Unit direction of each vector
    weights : np.ndarray(N)
        Weights that determine the size (length and/or thickness) of each vector
    """
    # check unit directions
    if not np.allclose(np.linalg.norm(unit_directions, axis=1), 1):
        raise ValueError('Given direction vectors are not unit vectors')

    # vmd/tcl function for constructing arrow
    output = '# Add function for vector field\n'
    output += ('proc vmd_draw_arrow {mol center unit_dir cyl_radius cone_radius length} {\n'
               'set start [vecsub $center [vecscale [vecscale 0.5 $length] $unit_dir]]\n'
               'set end [vecadd $start [vecscale $length $unit_dir]]\n'
               'set middle [vecsub $end [vecscale [vecscale 1.732050808 $cone_radius] $unit_dir]]\n'
               'graphics $mol cylinder $start $middle radius $cyl_radius\n'
               'graphics $mol cone $middle $end radius $cone_radius\n'
               '}\n'
               '#\n')

    def decompose_weight(weight):
        """ Decomposes a weight to the corresponding cylinder radius, cone radius and length

        Parameters
        ----------
        weight : float
            Weight of a vector

        Returns
        -------
        cyl_radius : float
            Radius of cylinder in vector
        cone_radius : float
            Radius of cone in vector
        length : float
            Length of vector
        """
        # FIXME: needs to be played around with
        return (0.3, 0.8, weight)

    for (center_x, center_y, center_z), (unit_x, unit_y, unit_z), weight in zip(centers,
                                                                                unit_directions,
                                                                                weights):
        output += ('draw arrow {{{0} {1} {2}}} {{{3} {4} {5}}} {6} {7} {8}\n'
                   ''.format(center_x, center_y, center_z, unit_x, unit_y, unit_z,
                             *decompose_weight(weight)))
    return output


def print_vmd_script_nci(scriptfile, densfile, rdgfile, isosurf=0.5, denscut=0.05):
    r"""
    Generate VMD (Visual Molecular Dynamics) script for visualizing NCI
    (non-covalent interactions) isosurfaces subject to the constraint of
    density(r) < denscut, i.e. low-density, and colored based on the
    sign(:math:`\lambda_2`) :math:`\rho`.

    Parameters
    ----------
    scriptfile : str
        Name of VMD script file to generate.
    densfile : str
        Name of density cube file.
    rdgfile : str
        Name of reduced density gradient cube file.
    isosurf : float, default=0.5
        Reduced density gradient isosurface to visualize.
    denscut : float, default=0.05
        Density cutoff used in creating reduced density gradient cube file.
        Similar to NCIPlot program, reduced density gradient of points with
        density > denscut will be set to 100.0 to display reduced density gradient
        isosurface subject to the constraint of low density.

    Note
    ----
    The script is the same as the one generated by NCIPlot software version 1.0.
    """
    output = _vmd_script_start()
    output += _vmd_script_molecule(densfile, rdgfile)
    output += _vmd_script_isosurface(isosurf=isosurf, scalemin=-denscut, scalemax=denscut, index=1)
    with open(scriptfile, 'w') as f:
        f.write(output)


def print_vmd_script_isosurface(scriptfile, isofile, colorfile=None, isosurf=0.5, material='Opaque',
                                scalemin=-0.05, scalemax=0.05, colorscheme='RGB', negative=False):
    """
    Generate VMD (Visual Molecular Dynamics) script for visualizing the isosurface based on
    one cube file when coloring by the value of another cube file on the isosurface.

    Parameters
    ----------
    scriptfile : str
        Name of VMD script file to generate.
    isofile : str
        Name of cube file used in VMD script for visualizing the isosurface.
    colorfile : str, default=None
        Name of cube file used in VMD script for coloring the isosurface.
        If None, the isofile is used for coloring.
    isosurf : float, default=0.5
        The value of the isosurface to visualize used in VMD script.
    material : str, default='Opaque'
        The material setting of the isosurface used in VMD script.

            Options: 'Opaque', 'Transparent', 'BrushedMetal', 'Diffuse', 'Ghost',
            'Glass1', 'Glass2', 'Glass3', 'Glossy', 'HardPlastic', 'MetallicPastel',
            'Steel', 'Translucent', 'Edgy', 'EdgyShiny', 'EdgyGlass', 'Goodsell',
            'AOShiny', 'AOChalky', 'AOEdgy', 'BlownGlass', 'GlassBubble', 'RTChrome'.

    scalemin : float, default=-0.05
        Smallest value to color on the isosurface used in VMD script.
    scalemax : float, default=0.05
        Largest value to color on the isosurface used in VMD script.
    colorscheme : str or int, default='RGB'
        Color scheme used in VMD script for coloring the isosurface.

            For color-scale, a sting is needed specifying the color scheme. The default is 'RGB',
            ans more options are available in table below.
            Alternatively, the isosurface can be colored with just one color.
            For this option, an integer is needed specifying the color.
            1057 colors are available in VMD, check the program or website for more details.
            If you want to plot the negative isosurface as well, a list of two colors can be given
            to specify the color of positive and negative iso-surfaces. e.g. [0, 1]

            =======  =====================================
            Options  Description
            =======  =====================================
            'RGB'    small=red, middle=green, large=blue
            'BGR'    small=blue, middle=green, large=red
            'RWB'    small=red, middle=white, large=blue
            'BWR'    small=blue, middle=white, large=red
            'RWG'    small=red, middle=white, large=green
            'GWR'    small=green, middle=white, large=red
            'GWB'    small=green, middle=white, large=blue
            'BWG'    small=blue, middle=white, large=green
            'BlkW'   small=black, large=white
            'WBlk'   small=white, large=black
            =======  =====================================
    negative : bool, default=False
        Determines if you want to plot the negative of the isosurface as well. The default is false.
    """
    if not negative and isinstance(colorscheme, (int, str)):
        pos_color, neg_color = colorscheme, colorscheme
    elif negative and hasattr(colorscheme, '__iter__') and len(colorscheme) == 2:
        pos_color, neg_color = colorscheme
    else:
        raise TypeError('To plot the negative of the isosurface with integer color scheme, '
                        'two integers must be given. The first integer corresponds to the '
                        'ColorID for the positive surface, and the second for the negative '
                        'surface.'
                        'If the negative is not plotted, then colorscheme must be given as an '
                        'integer or a string')

    output = _vmd_script_start()
    if colorfile is not None:
        output += _vmd_script_molecule(colorfile, isofile)
        file_index = 1
    else:
        output += _vmd_script_molecule(isofile)
        file_index = 0

    output += _vmd_script_isosurface(isosurf=isosurf, index=file_index, material=material,
                                     scalemin=scalemin, scalemax=scalemax, colorscheme=pos_color)

    if negative:
        output += _vmd_script_isosurface(isosurf=-isosurf, index=file_index, material=material,
                                         scalemin=scalemin, scalemax=scalemax,
                                         colorscheme=neg_color)

    with open(scriptfile, 'w') as f:
        f.write(output)


def print_vmd_script_multiple_cube(scriptfile, cubes, isosurfs=None, material='Opaque',
                                   scalemin=-0.05, scalemax=0.05, colors=None):
    """ Generate VMD (Visual Molecular Dynamics) script for visualizing multiple cube files
    simultaneously where data from each cube file is colored differently

    Parameters
    ----------
    scriptfile : str
        Name of VMD script file to generate.
    cubes : list of str
        Names of cube files to plot
    isosurfs : float, list of float
        Isovalue at which the plot (isosurface) is generated
        If a float is given, then this is the value of isosurface for all cube files
        Default value is 0.5 for all isosurfaces
    material : str
        The material setting of the isosurface used in VMD script.
        One of 'Opaque', 'Transparent', 'BrushedMetal', 'Diffuse', 'Ghost',
        'Glass1', 'Glass2', 'Glass3', 'Glossy', 'HardPlastic', 'MetallicPastel',
        'Steel', 'Translucent', 'Edgy', 'EdgyShiny', 'EdgyGlass', 'Goodsell',
        'AOShiny', 'AOChalky', 'AOEdgy', 'BlownGlass', 'GlassBubble', 'RTChrome'.
        Default is 'Opaque'
    scalemin : float
        Smallest value to color on the isosurface used in VMD script.
        Default is -0.05
    scalemax : float
        Largest value to color on the isosurface used in VMD script.
        Default is 0.05
    colors : list of int
        Colors of each cube file data
        Each integer corresponds to a color. See VMD program or manual for details.
        Default selects random color for each cube file

    Note
    ----
    Not quite sure what happens when the number of cube files exceeds 33 (possiblly the  maximum
    number of ColorID's in VMD)

    Raises
    ------
    TypeError
        If cube files are not provided as a list or tuple
        If colors are not provided as a list or tuple of the same length as the cube files
    ValueError
        If any of the cube files cannot be found
        If any of the colors are not an integer between 0 and 32
    """
    if not isinstance(cubes, (list, tuple)):
        raise TypeError('The cube files must be given as a list or tuple')
    if not all(os.path.isfile(cube) for cube in cubes):
        raise ValueError('Cannot find at least one of the cube files')

    if isosurfs is None:
        isosurfs = [0.5 for i in cubes]
    elif isinstance(isosurfs, float):
        isosurfs = [isosurfs for i in cubes]
    if not (isinstance(isosurfs, (list, tuple)) and len(isosurfs) != len(cubes)):
        raise TypeError('The isosurfs must be provided as a list or tuple of same length as the '
                        'number of cube files')
    elif not all(isinstance(isosurf, float) for isosurf in isosurfs):
        raise TypeError('Each isosurface value must be a float')

    if colors is None:
        colors = range(len(cubes))
    elif not (isinstance(colors, (list, tuple)) and len(colors) != len(cubes)):
        raise TypeError('The colors must be provided as a list or tuple of the same length as the '
                        'number of cube files')
    elif not all(isinstance(color, int) and 0 <= color < 33 for color in colors):
        raise ValueError('Each color must be given as an integer between 0 and 32')

    output = _vmd_script_start()
    output += _vmd_script_molecule(*cubes)
    for i, (cube, isosurf, color) in enumerate(zip(cubes, isosurfs, colors)):
        output += _vmd_script_isosurface(isosurf=isosurf, index=i, material=material,
                                         scalemin=scalemin, scalemax=scalemax, colorscheme=color)

    with open(scriptfile, 'w') as f:
        f.write(output)
