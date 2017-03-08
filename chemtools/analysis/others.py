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
'''Module for Density-Based and Orbital-Based Analysis of Quantum Chemistry Output Files.

   This modules contains wrappers which take outputs of quantum chemistry softwares and
   compute various descriptive tools based on the density and orbital information.
'''

import numpy as np
from horton import IOData
from chemtools.tool.densitytool import DensityLocalTool
from chemtools.utils import CubeGen
from chemtools.analysis.output import _print_vmd_script_nci
import matplotlib.pyplot as plt


class NCI(object):
    '''
    Class for density-based analysis of one quantum chemistry output file.
    '''
    def __init__(self, iodata):
        '''
        Parameters
        ----------
        iodata : IOData
            The molecule's IOData.
        '''
        if not isinstance(iodata, IOData):
            raise ValueError('Argument cube should be an instance of IOData!')

        self._mol = iodata

    @classmethod
    def from_file(cls, filename):
        '''
        Initialize class from files.
        '''
        # case of one file not given as a list
        if isinstance(filename, (str, unicode)):
            return cls(IOData.from_file(filename))
        # case of list of file(s)
        for file in filename:
            raise ValueError('Multiple files are not suported')

    @classmethod
    def from_iodata(cls, iodata):
        '''
        Initialize class from `IOData` objects.
        '''
        # This classmethod is redundant, but just included for uniformity
        return cls(iodata)

    def dump_files(self, filename, isosurf=0.50, denscut=0.05, cube=None, plot=False):
        r'''
        Generate density and reduced density gradient cube files, as well as a VMD (Visual
        Molecular Dynamics) script to visualize non-covalnet inteactions (NCI).

        Parameters
        ----------
        filename : str
            Name of generated cube files and vmd script.
        isosurf : float, default=0.5
            Value of reduced density gradient (RDG) iso-surface used in VMD script.
        denscut : float, default=0.05
            Density cutoff used in creating reduced density gradient cube file.
            Similar to NCIPlot program, reduced density gradient of points with
            density > denscut will be set to 100.0 to display reduced density gradient
            iso-surface subject to the constraint of low density.
        cube : instance of `CubeGen`, default=None
            Cubic grid used for calculating and visualizating the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        plot : boolean, default=False
            Generate a plot of reduced density gradient RDG vs. sign(:math:`\lambda_2`):math:`\rho`.

        Note
        ----
        The generated cube files and script imitate the NCIPlot software version 1.0.
        '''
        # Generate or check cubic grid
        if cube is None:
            cube = CubeGen.from_molecule(self._mol.numbers, self._mol.pseudo_numbers,
                                         self._mol.coordinates, spacing=0.1, threshold=2.0)
        elif not isinstance(cube, CubeGen):
            raise ValueError('Argument cube should be an instance of CubeGen!')

        # Compute density & gradient on cubic grid
        dm_full = self._mol.get_dm_full()
        dens = self._mol.obasis.compute_grid_density_dm(dm_full, cube.points)
        grad = self._mol.obasis.compute_grid_gradient_dm(dm_full, cube.points)
        # Calculate reduced density gradient on cubic grid
        localtool = DensityLocalTool(dens, grad, hessian=None)
        rdgrad = localtool.reduced_density_gradient
        # Compute hessian and its eigenvalues on cubuc grid
        hess = _compute_hessian(self._mol, cube.points)
        eigvalues = np.linalg.eigvalsh(hess)
        # Use sign of second eigenvalue to distinguish interaction types
        dens = np.sign(eigvalues[:, 1]) * dens

        # Plot reduced density gradient vs. sign(lambda_2)*density
        if plot:
            plt.scatter(dens, rdgrad, lw = 0.5)
            plt.xlim(-0.2, 0.2)
            plt.ylim(0.0, 2.0)
            plt.xlabel(r'sign($\lambda_2$)$\rho$ (a.u)')
            plt.ylabel('Reduced Density Gradient')
            plt.show()

        # Similar to NCIPlot program, reduced density gradient of points with
        # density > cutoff will be set to 100.0 before generating cube file to
        # display reduced density gradient iso-surface subject to the constraint
        # of low density, i.e. density < denscut.
        rdgrad[abs(dens) > denscut] = 100.0
        # Similar to NCIPlot program, sign(hessian second eigenvalue)*density is
        # multiplied by 100.0 before genetaing cube file used for coloring the
        # reduced density gradient iso-surface.
        dens *= 100.0

        # Name of output files:
        densfile = filename + '-dens.cube'    # density cube file
        rdgfile  = filename + '-grad.cube'    # reduced density gradient cube file
        vmdfile  = filename + '.vmd'          # reduced vmd script file
        # Dump density & reduced density gradient cube files
        cube.dump_cube(densfile, dens)
        cube.dump_cube(rdgfile, rdgrad)
        # Make VMD scripts for visualization
        _print_vmd_script_nci(vmdfile, densfile, rdgfile, isosurf, denscut*100.0)


def _compute_hessian(mol, points):
    '''
    Compute hessian of electron density defined as the second-order partial
    derivative of electron density w.r.t. coordinates.

    Parameters
    ----------
    mol : instance of `IOData`
        An `IOData` object.
    points : np.ndarray
        Coordinates of grid points for calculating hessian.

    Note
    ----
    This finite difference implementation is temporary until hessian is implemented in HORTON.
    '''
    # Make small change in coordinates along x, y, & z directions
    eps = 1.0e-5
    dx = points + np.array([[eps, 0.0, 0.0]])
    dy = points + np.array([[0.0, eps, 0.0]])
    dz = points + np.array([[0.0, 0.0, eps]])
    # Calculate gradient on original grid points
    dm = mol.get_dm_full()
    gpnt = mol.obasis.compute_grid_gradient_dm(dm, points)
    # Calculate hessian using finite difference
    hess = np.zeros((points.shape[0], 3, 3), float)
    hess[:,:,0] = (mol.obasis.compute_grid_gradient_dm(dm, dx) - gpnt) / eps
    hess[:,:,1] = (mol.obasis.compute_grid_gradient_dm(dm, dy) - gpnt) / eps
    hess[:,:,2] = (mol.obasis.compute_grid_gradient_dm(dm, dz) - gpnt) / eps
    # Make Hessian symmetric to avoid complex eigenvalues
    hess = 0.5 * (hess + np.transpose(hess, axes=(0, 2, 1)))
    return hess
