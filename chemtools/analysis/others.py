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
"""Module for Density-Based and Orbital-Based Analysis of Quantum Chemistry Output Files.

   This modules contains wrappers which take outputs of quantum chemistry software and
   compute various descriptive tools based on the density and orbital information.
"""

import numpy as np
from horton import IOData
from chemtools.tool.densitytool import DensityLocalTool
from chemtools.utils import CubeGen
from chemtools.tool.orbitaltool import OrbitalLocalTool
from chemtools.analysis.output import print_vmd_script_nci
import matplotlib.pyplot as plt


class NCI(object):
    """
    Class for the Non-Covalent Interactions (NCI).
    """
    def __init__(self, density, rdgradient, cube, hessian=None):
        """
        Parameters
        ----------
        cube : instance of `CubeGen`, default=None
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        if hessian is not None:
            # Compute hessian and its eigenvalues on cubuc grid
            eigvalues = np.linalg.eigvalsh(hessian)

            # Use sign of second eigenvalue to distinguish interaction types
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
        """
        Initialize class from file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        cube : instance of `CubeGen`, default=None
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        # case of one file not given as a list
        if isinstance(filename, (str, unicode)):
            return cls.from_iodata(IOData.from_file(filename), cube)
        # case of list of file(s)
        for _ in filename:
            raise ValueError('Multiple files are not supported')

    @classmethod
    def from_iodata(cls, iodata, cube=None):
        """
        Initialize class from ``IOData`` object from ``HORTON`` library.

        Parameters
        ----------
        iodata : ``IOData``
            Instance of ``IOData``.
        cube : instance of `CubeGen`, default=None
            Cubic grid used for calculating and visualizing the NCI.
            If None, it is constructed from molecule with spacing=0.1 and threshold=2.0
        """
        # Generate or check cubic grid
        if cube is None:
            cube = CubeGen.from_molecule(iodata.numbers, iodata.pseudo_numbers,
                                         iodata.coordinates, spacing=0.1, threshold=2.0)
        elif not isinstance(cube, CubeGen):
            raise ValueError('Argument cube should be an instance of CubeGen!')

        # Compute density & gradient on cubic grid
        dm_full = iodata.get_dm_full()
        dens = iodata.obasis.compute_grid_density_dm(dm_full, cube.points)
        grad = iodata.obasis.compute_grid_gradient_dm(dm_full, cube.points)
        # Compute hessian on cubic grid
        hess = _compute_hessian(iodata, cube.points)

        # initialize DensityLocalTool & compute reduced gradient
        temp = DensityLocalTool(dens, grad)
        rdgrad = temp.reduced_density_gradient

        return cls(dens, rdgrad, cube, hessian=hess)

    @property
    def signed_density(self):
        r"""
        Electron density :math:`\rho\left(\mathbf{r}\right)` evaluated on a grid,
        signed by the second eigenvalue of the Hessian at that point, i.e.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.
        """
        return self._signed_density

    @property
    def eigvalues(self):
        r"""
        Eigenvalues of Hessian.
        """
        return self._eigvalues

    def plot(self, filename):
        r"""
        Plot reduced density gradient vs.
        :math:`\text{sgn}\left(\lambda_2\right) \times \rho\left(\mathbf{r}\right)`.
        """
        plt.scatter(self._signed_density, self._rdgrad, lw=0.5)
        plt.xlim(-0.2, 0.2)
        plt.ylim(0.0, 2.0)
        plt.xlabel(r'sgn($\lambda_2$)$\rho$ (a.u)')
        plt.ylabel('Reduced Density Gradient')
        plt.savefig(filename)

    def dump_files(self, filename, isosurf=0.50, denscut=0.05):
        r"""
        Generate density and reduced density gradient cube files, as well as a VMD (Visual
        Molecular Dynamics) script to visualize non-covalent interactions (NCI).

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

        Note
        ----
        The generated cube files and script imitate the NCIPlot software version 1.0.
        """
        # Similar to NCIPlot program, reduced density gradient of points with
        # density > cutoff will be set to 100.0 before generating cube file to
        # display reduced density gradient iso-surface subject to the constraint
        # of low density, i.e. density < denscut.
        cutrdg = np.array(self._rdgrad, copy=True)
        cutrdg[abs(self._density) > denscut] = 100.0
        # Similar to NCIPlot program, sign(hessian second eigenvalue)*density is
        # multiplied by 100.0 before generating cube file used for coloring the
        # reduced density gradient iso-surface.
        if self._signed_density is not None:
            dens = 100.0 * self._signed_density
        else:
            dens = 100.0 * self._density

        # Name of output files:
        densfile = filename + '-dens.cube'    # density cube file
        rdgfile = filename + '-grad.cube'    # reduced density gradient cube file
        vmdfile = filename + '.vmd'          # vmd script file
        # Dump density & reduced density gradient cube files
        self._cube.dump_cube(densfile, dens)
        self._cube.dump_cube(rdgfile, cutrdg)
        # Make VMD scripts for visualization
        print_vmd_script_nci(vmdfile, densfile, rdgfile, isosurf, denscut * 100.0)


def _compute_hessian(mol, points):
    """
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
    """
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
    hess[:, :, 0] = (mol.obasis.compute_grid_gradient_dm(dm, dx) - gpnt) / eps
    hess[:, :, 1] = (mol.obasis.compute_grid_gradient_dm(dm, dy) - gpnt) / eps
    hess[:, :, 2] = (mol.obasis.compute_grid_gradient_dm(dm, dz) - gpnt) / eps
    # Make Hessian symmetric to avoid complex eigenvalues
    hess = 0.5 * (hess + np.transpose(hess, axes=(0, 2, 1)))
    return hess


class OrbitalAnalysis(OrbitalLocalTool):
    """
    Class for orbital-based analysis.
    """
    def __init__(self, points, obasis, exp_alpha, exp_beta=None):
        super(OrbitalAnalysis, self).__init__(points, obasis, exp_alpha, exp_beta=exp_beta)

    @classmethod
    def from_file(cls, filename, points):
        """
        Initialize class from file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        # case of one file not given as a list
        if isinstance(filename, (str, unicode)):
            mol = IOData.from_file(filename)
            if hasattr(mol, 'exp_beta'):
                return cls(points, mol.obasis, mol.exp_alpha, mol.exp_beta)
            else:
                return cls(points, mol.obasis, mol.exp_alpha)
        # case of list of file(s)
        for _ in filename:
            raise ValueError('Multiple files are not supported')

    @classmethod
    def from_iodata(cls, iodata, points):
        """
        Initialize class from `IOData` object.

        Parameters
        ----------
        iodata : `IOData`
            Instance of `IOData`.
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        # check if iodata has exp_beta
        if hasattr(iodata, 'exp_beta'):
            return cls(points, iodata.obasis, iodata.exp_alpha, iodata.exp_beta)
        else:
            return cls(points, iodata.obasis, iodata.exp_alpha)
