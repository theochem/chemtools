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
# pragma pylint: disable=wrong-import-position
"""Module for Non-Covalent Interactions (NCI) analysis of Quantum Chemistry Output Files."""


import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from matplotlib import rcParams

from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densbased import DensityBasedTool
from chemtools.utils.cube import CubeGen
from chemtools.outputs.output_vmd import print_vmd_script_nci


__all__ = ['NCI']


class NCI(object):
    """Non-Covalent Interactions (NCI) Class."""

    def __init__(self, density, rdgradient, cube, hessian=None):
        """Initialize class using density, reduced density gradient and `CubeGen` instance.

        Parameters
        ----------
        density : np.array
            Density evaluated on grid points of `cube`.
        rdgradient : np.array
            Reduced density gradient evaluated on grid points of `cube`
        cube : instance of `CubeGen`, optional
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        hessian : np.array, optional
            Hessian of density evaluated on grid points of `cube`. This is a array with shape
            (n, 6) where n is the number of grid points of `cube`.
        """
        if density.shape != (len(cube.points),):
            raise ValueError('Shape of density argument {0} does not match '
                             'expected ({1},) shape.'.format(density.shape, len(cube.points)))
        if rdgradient.shape != (len(cube.points),):
            raise ValueError('Shape of rdgradient argument {0} does not '
                             'match expected ({1},) shape.'.format(density.shape, len(cube.points)))

        if hessian is not None:
            if hessian.shape != (len(cube.points), 6):
                raise ValueError('Shape of hessian argument {0} does not match expected ({1}, 6)'
                                 ' shape.'.format(hessian.shape, len(cube.points)))

            # convert the (n, 6) shape to (n, 3, 3) to calculate eigenvalues.
            hestri = np.zeros((len(cube.points), 3, 3))
            tmp = np.zeros((3, 3))
            for i in range(0, len(cube.points)):
                tmp[np.triu_indices(3)] = hessian[i, :]
                hestri[i, :] = tmp

            # compute hessian and its eigenvalues on cubic grid
            eigvalues = np.linalg.eigvalsh(hestri, UPLO='U')

            # use sign of second eigenvalue to distinguish interaction types
            sdens = np.sign(eigvalues[:, 1]) * density

            self._signed_density = sdens
            self._eigvalues = eigvalues

        else:
            self._signed_density = None
            self._eigvalues = None

        self._density = density
        self._rdgrad = rdgradient
        self._cube = cube

    @classmethod
    def from_file(cls, filename, cube=None):
        """Initialize class using wave-function file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        cube : instance of `CubeGen`, optional
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        molecule = Molecule.from_file(filename)
        return cls.from_molecule(molecule, cube)

    @classmethod
    def from_molecule(cls, molecule, cube=None):
        """Initialize class from ``Molecule`` object.

        Parameters
        ----------
        molecule : ``Molecule``
            Instance of ``Molecule``.
        cube : instance of `CubeGen`, optional
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        # generate or check cubic grid
        if cube is None:
            cube = CubeGen.from_molecule(molecule.numbers, molecule.pseudo_numbers,
                                         molecule.coordinates, spacing=0.1, threshold=2.0)
        elif not isinstance(cube, CubeGen):
            raise ValueError('Argument cube should be an instance of CubeGen!')
        # compute density, gradient & hessian on cubic grid
        dens = molecule.compute_density(cube.points)
        grad = molecule.compute_gradient(cube.points)
        hess = molecule.compute_hessian(cube.points)
        # compute reduced gradient
        rdgrad = DensityBasedTool(dens, grad).reduced_density_gradient
        return cls(dens, rdgrad, cube, hessian=hess)

    @property
    def signed_density(self):
        r"""Signed electron density.

        Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on a grid,
        signed by the second eigenvalue of the Hessian at that point, i.e.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.
        """
        return self._signed_density

    @property
    def eigvalues(self):
        r"""Eigenvalues of Hessian."""
        return self._eigvalues

    def plot(self, filename, color='b'):
        r"""Plot reduced density gradient.

        Reduced density gradient vs.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.

        Parameters
        ----------
        filename : str
            Name of generated 2D plot.

            If the given filename does not have a proper extension (representing its format),
            the 'png' format is used by default (i.e. plot is saved as filename.png).

            Supported formats (which should be specified as filename extensions) include:

            - 'svgz' or 'svg' (Scalable Vector Graphics)
            - 'tif' or 'tiff' (Tagged Image File Format)
            - 'raw' (Raw RGBA bitmap)
            - 'png' (Portable Network Graphics)
            - 'ps' (Postscript)
            - 'eps' (Encapsulated Postscript)
            - 'rgba' (Raw RGBA bitmap)
            - 'pdf' (Portable Document Format)

        color : str, optional
            Color of plot. Default is blue specified with 'b'.
            For details on specifying colors, please refer to
            http://matplotlib.org/users/colors.html
        """
        # set font
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Times New Roman']
        rcParams['mathtext.fontset'] = 'stix'
        # create figure
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        # scatter plot
        plt.scatter(self._signed_density, self._rdgrad, marker='o', color=color)
        # set axis range and label
        plt.xlim(-0.2, 0.2)
        plt.ylim(0.0, 2.0)
        plt.xlabel(r'sgn$\mathbf{(\lambda_2)}$ $\times$ $\mathbf{\rho(r)}$ (a.u)',
                   fontsize=12, fontweight='bold')
        plt.ylabel('Reduced Density Gradient', fontsize=12, fontweight='bold')
        # hide the right, top and bottom spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        # save plot ('.png' extension is added by default, if filename is not a supported format)
        plt.savefig(filename, dpi=800)

    def dump_files(self, filename, isosurf=0.50, denscut=0.05):
        r"""Generate cube files and VMD script to visualize non-covalent interactions (NCI).

        Generate density and reduced density gradient cube files, as well as a VMD (Visual
        Molecular Dynamics) script to visualize non-covalent interactions (NCI).

        Parameters
        ----------
        filename : str
            Name of generated cube files and vmd script.
        isosurf : float, optional
            Value of reduced density gradient (RDG) iso-surface used in VMD script.
        denscut : float, optional
            Density cutoff used in creating reduced density gradient cube file.
            Similar to NCIPlot program, reduced density gradient of points with
            density > denscut will be set to 100.0 to display reduced density gradient
            iso-surface subject to the constraint of low density.
            To visualize all reduced density gradient iso-surfaces, disregarding of the
            corresponding density value, set this argument equal to infity using `float('inf')`.
        Note
        ----
        The generated cube files and script imitate the NCIPlot software version 1.0.
        """
        # similar to NCIPlot program, reduced density gradient of points with
        # density > cutoff will be set to 100.0 before generating cube file to
        # display reduced density gradient iso-surface subject to the constraint
        # of low density, i.e. density < denscut.
        cutrdg = np.array(self._rdgrad, copy=True)
        cutrdg[abs(self._density) > denscut] = 100.0

        # similar to NCIPlot program, sign(hessian second eigenvalue)*density is
        # multiplied by 100.0 before generating cube file used for coloring the
        # reduced density gradient iso-surface.
        if self._signed_density is not None:
            dens = 100.0 * self._signed_density
        else:
            dens = 100.0 * self._density

        # name of output files
        densfile = filename + '-dens.cube'    # density cube file
        rdgfile = filename + '-grad.cube'     # reduced density gradient cube file
        vmdfile = filename + '.vmd'           # vmd script file
        # dump density & reduced density gradient cube files
        self._cube.dump_cube(densfile, dens)
        self._cube.dump_cube(rdgfile, cutrdg)
        # write VMD scripts
        print_vmd_script_nci(vmdfile, densfile, rdgfile, isosurf, denscut * 100.0)
