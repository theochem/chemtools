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
"""Module for constructing different mesh (grid coordinates)."""


import numpy as np


def plane_mesh(points, spacing, extension):
    """Return the grid points on the plane spanned by the given three coordinates.

    Parameters
    ----------
    points : np.ndarray(3, 3)
        Points that are on the plane.
        Rows corespond to the different points, columns correspond to the x, y, and z components.
        The first set of coordinate will be used to establish the new vertical (i.e. y) direction in
        the plane.
    spacing : float
        Upper bound to the spacing between adjacent grid points. This means that the spacing between
        adjacent points will be less than the given value.
    extension : float
        Distance from the center of the three points that will define the edges of the mesh.

    Returns
    -------
    grid_points : np.ndarray(N_v, N_h, 3)
        Points on the plane spanned by the given three points.
        First index corresponds to the position along the vertical direction of the plane mesh. This
        direction is the same as the direction of the first coordinate wrt origin. `N_v` is the
        number of points vertically along the the mesh.
        Second index corresponds to the position along the horizontal direction of the plane mesh.
        `N_v` is the number of points horizontally along the the mesh.

    Raises
    ------
    TypeError
        If `plane_points` is not two-dimensional numpy array of shape (3, 3).
        If `spacing` is not a float.
        If `extension` is not int or float.
    ValueError
        If `spacing` is less than or equal to zero.
        If `extension` is less than or equal to zero.
        If the three points on the plane are on a line.

    """
    if not (isinstance(points, np.ndarray) and points.ndim == 2 and points.shape == (3, 3)):
        raise TypeError("Arguments points must be given as a 2D numpy array of shape (3, 3).")
    if not isinstance(spacing, float) or spacing <= 0:
        raise TypeError("Argument spacing must be a float greater than 0.")
    if not isinstance(extension, (int, float)) or extension < 0:
        raise TypeError("Argument extension must be int or float greater than 0.")

    center = np.average(points, axis=0)
    vec = points - center
    length_vec = (np.sum(vec ** 2, axis=1) ** 0.5)[:, None]
    if np.any(length_vec == 0):
        raise ValueError("Three points on the plane cannot be in a line.")
    unit_vec = vec / length_vec
    if np.unique(unit_vec, axis=0).shape[0] < 3:
        raise ValueError("Three points on the plane cannot be in a line.")

    # edges of the cube
    edges = unit_vec * extension + points
    # length of the edges from the center
    length_edges = (np.sum((edges - center) ** 2, axis=1) ** 0.5)[:, None]
    # use first point as a reference
    angle_01 = np.arccos(unit_vec[0].dot(unit_vec[1]))
    angle_02 = np.arccos(unit_vec[0].dot(unit_vec[2]))
    height_1 = length_edges[1] * np.sin(angle_01 - np.pi / 2)
    width_1 = length_edges[1] * np.cos(angle_01 - np.pi / 2)
    height_2 = length_edges[2] * np.sin(angle_02 - np.pi / 2)
    width_2 = length_edges[2] * np.cos(angle_02 - np.pi / 2)

    height = max(height_1, height_2) + length_edges[0]
    width = width_1 + width_2
    unit_vertical = unit_vec[0]
    origin = center + unit_vertical * length_edges[0]
    if height_1 > height_2:
        vec_horizontal = edges[1] - (center - height_1 * unit_vertical)
    else:
        vec_horizontal = edges[2] - (center - height_2 * unit_vertical)
    origin += vec_horizontal
    unit_horizontal = vec_horizontal / np.sum(vec_horizontal ** 2) ** 0.5

    output = np.linspace(
        origin, origin - unit_horizontal * width, num=int(width / spacing) + 2, endpoint=True
    )
    output = np.linspace(
        output, output - unit_vertical * height, num=int(height / spacing) + 2, endpoint=True
    )
    return output
