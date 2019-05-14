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
"""Simple Plotting Module."""


import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from matplotlib import rcParams


__all__ = ['plot_scatter']


def plot_scatter(x, y, fname, color='b', xlabel=None, ylabel=None, xlim=None, ylim=None):
    r"""Scatter plot of y versus x.

    Parameters
    ----------
    x : 1-D array or sequence.
        Array or sequence containing data on x axis.
    y : 1-D array or sequence.
        Array or sequence containing data on y axis.
    fname : str
        A string representing the path to a filename for storing the plot.
        If the given filename does not have a proper extension, the 'png' format is used
        by default, i.e. plot is saved as filename.png.
        Supported formats, which can be specified as filename extensions, include:

        - 'svgz' or 'svg' (Scalable Vector Graphics)
        - 'tif' or 'tiff' (Tagged Image File Format)
        - 'raw' (Raw RGBA bitmap)
        - 'png' (Portable Network Graphics)
        - 'ps' (Postscript)
        - 'eps' (Encapsulated Postscript)
        - 'rgba' (Raw RGBA bitmap)
        - 'pdf' (Portable Document Format)

    color : str, optional
        Color of plot. To customize color, see http://matplotlib.org/users/colors.html
    xlabel : str, optional
        The x axis label.
    ylabel : str, optional
        The y axis label.
    xlim : 1-D array or sequence of length 2, optional
        The lower and higher limit of x axis.
    ylim : 1-D array or sequence of length 2, optional
        The lower and higher limit of y axis.

    """
    # set font
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = ['Times New Roman']
    rcParams['mathtext.fontset'] = 'stix'
    # create figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # scatter plot
    if len(x) != len(y):
        raise ValueError('Length of x & y does not match! {0}!={1}'.format(len(x), len(y)))
    plt.scatter(x, y, marker='o', color=color)
    # set axis range
    if xlim:
        if len(xlim) != 2:
            raise ValueError('Argument xlim={0} should have a length 2!'.format(len(xlim)))
    plt.xlim(*xlim)
    if ylim:
        if len(ylim) != 2:
            raise ValueError('Argument ylim={0} should have a length 2!'.format(len(ylim)))
        plt.ylim(*ylim)
    # set axis label
    if xlabel:
        plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    if ylabel:
        plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    # hide the right, top and bottom spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
    # save plot ('.png' extension is added by default, if filename is not a supported format)
    plt.savefig(fname, dpi=800)
