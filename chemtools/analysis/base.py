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
'''Analyze Quantum Chemistry Output Files Module.'''


import numpy as np
from horton import IOData
from chemtools.tool.globaltool import LinearGlobalTool, QuadraticGlobalTool, ExponentialGlobalTool, RationalGlobalTool
from chemtools.tool.densitytool import DensityLocalTool
from chemtools.utils import CubeGen
from chemtools.analysis.output import _print_vmd
import matplotlib.pyplot as plt

class Analyze_1File(object):
    '''
    Class for analyzing one quantum chemistry output file.
    '''
    def __init__(self, molecule_filename, model='quadratic', energy_expr=None):
        '''
        Parameters
        ----------
        molecule_filename : str
            The path to the molecule's file.
        model : str, default='quadratic'
            Energy model used to calculate descriptive tools.
            The available models include:
            * 'linear'; refer to :py:class:`chemtools.tool.globaltool.LinearGlobalTool` for more information.
            * 'quadratic'; refer to :py:class:`chemtools.tool.globaltool.QuadraticGlobalTool` for more information.
            * 'exponential'; refer to :py:class:`chemtools.tool.globaltool.ExponentialGlobalTool` for more information.
            * 'rational'; refer to :py:class:`chemtools.tool.globaltool.RationalGlobalTool` for more information.
            * 'general'; refer to :py:class:`chemtools.tool.globaltool.GeneralGlobalTool` for more information.
            If 'general' model is selected, an energy expression should be given.
        energy_expr : ``Sympy.expr``, default=None
            Energy expresion used, if 'general' model is selected.
        '''
        mol = IOData.from_file(molecule_filename)
        self._mol = mol
        if model not in ['linear', 'quadratic', 'exponential', 'rational', 'general']:
            raise ValueError('Argument model={0} is not supported.'.format(model))
        if model is 'general' and energy_expr is None:
            raise ValueError('Argument energy_expr is required when model=\'general\'.')
        self._model = model
        self.energy_expr = energy_expr

        # TODO: Some attributes of the self._mol should become the class attribute
        #       like coordinates, numbers, energy, etc.

        # Get E(HOMO), E(LUMO) & number of electrons
        homo_energy = self._mol.exp_alpha.homo_energy
        lumo_energy = self._mol.exp_alpha.lumo_energy
        n_elec = int(np.sum(self._mol.exp_alpha.occupations))

        # HACK: self._mol might not have the exp_beta attribute, making it crash
        if hasattr(self._mol, 'exp_beta'):
            if self._mol.exp_beta is not None:
                n_elec += int(np.sum(self._mol.exp_beta.occupations))
                if self._mol.exp_beta.homo_energy > homo_energy:
                    homo_energy = self._mol.exp_beta.homo_energy
                if self._mol.exp_beta.lumo_energy < lumo_energy:
                    lumo_energy = self._mol.exp_beta.lumo_energy

        # HACK: wfn files don't have an energy attribute...
        if hasattr(self._mol, 'energy'):
            # Compute E(N), E(N+1), & E(N-1)
            energy_zero = self._mol.energy
            energy_plus = energy_zero - lumo_energy
            energy_minus = energy_zero - homo_energy

            # Define global tool
            if model == 'linear':
                self._globaltool = LinearGlobalTool(energy_zero, energy_plus, energy_minus, n_elec)
            elif model == 'quadratic':
                self._globaltool = QuadraticGlobalTool(energy_zero, energy_plus, energy_minus, n_elec)
            elif model == 'exponential':
                self._globaltool = ExponentialGlobalTool(energy_zero, energy_plus, energy_minus, n_elec)
            elif model == 'rational':
                self._globaltool = RationalGlobalTool(energy_zero, energy_plus, energy_minus, n_elec)
            elif model == 'general':
                pass

        # ???: do we really want to compute the density here? if so, on what grid? Becke vs Cube?
        # Compute electron density (of the N electron system)
        # density = None

        # Compute gradient of electron density
        # gradient = None

        # Compute Hessian of electron density
        # hessian = None

        # Define density-based local tool
        # self._localtool = DensityLocalTool(density, gradient, hessian)


    @property
    def model(self):
        '''
        Energy model used to calculate descriptive tools.
        '''
        return self._model

    @property
    def globaltool(self):
        '''
        Instance of one of the gloabl reactivity tool classes.
        '''
        return self._globaltool

    # @property
    # def localtool(self):
    #     '''
    #     Instance of one of the local reactivity tools classes.
    #     '''
    #     return self._localtool

    def compute_nci(self, filename, isosurf=0.50, denscut=0.05, cube=None, plot=False):
        r'''
        Generate all files to plot NCI.
        
        Parameters
        ----------
        filename : str
            Name for the outputfiles: *-dens.cube, *-grad.cube and *.vmd.
        isosurf : float, default=0.5
            Value of reduced density gradient (RDG) isosurface used in VMD script.
        denscut : float, default=0.5
            Density cutoff used in creating cube RDG file. 
            RDG of points with density > cutoff will be set to 100.
        cube : instance of CubeGen, default=None
            Cubic grid used to calculate NCI. 
            By deafault it is constructed from the molecule with spacing=0.1 and threshold=2.0
        plot : boolean, default=False
            Generate a plot of RDG vs. sign(:math:`\lambda_2`):math:`\rho`.
        '''
        #
        # the files used:
        #
        densfile = filename + '-dens.cube'    # name of the density cube file
        rdgfile  = filename + '-grad.cube'    # name of the reduced density gradient cube file
        vmdfile  = filename + '.vmd'          # name of the reduced vmd script file

        if cube is None:
            cube = CubeGen.from_molecule(self._mol.numbers, self._mol.pseudo_numbers, self._mol.coordinates, spacing=0.1,threshold=2.0)
        else:
            if not isinstance(cube, CubeGen):
                raise ValueError('cube should be an instance of CubeGen!')

        dm_full = self._mol.get_dm_full()

        dens = self._mol.obasis.compute_grid_density_dm(dm_full, cube.gridpoints)
        grad = self._mol.obasis.compute_grid_gradient_dm(dm_full, cube.gridpoints)
        dloctool = DensityLocalTool(dens, grad)
        rdgrad = dloctool.reduced_density_gradient

        hess = _compute_hessian(self._mol, cube.gridpoints)
        dens = np.sign(np.linalg.eigvalsh(hess)[:,1])*dens

        if plot:
            # plotting the reduced density gradient vs. sign(lambda_2)density
            # FIXME: this might need to go, but I would like to keep it for testing & demo for now.
            plt.xlim(-0.2, 0.2)
            plt.ylim(0.0, 2.0)
            plt.scatter(dens, rdgrad, lw = 0.5)
            plt.xlabel(r"sign($\lambda_2$)$\rho$ (a.u)")
            plt.ylabel(r"Reduced gradient")
            plt.show()

        #mask rdgrad:
        rdgrad[abs(dens) > denscut] = 100.0

        cube.dump_cube(densfile ,dens*100.0)
        cube.dump_cube(rdgfile,rdgrad)
        _print_vmd(vmdfile, densfile, rdgfile, isosurf, denscut*100.0)


def _compute_hessian(mol, x):
    #untill the electronic Hessian is included in HORTON, we do it by finite difference:
    dm = mol.get_dm_full()
    eps = 1.0e-5
    dx = x + np.array([[eps, 0.0, 0.0]])
    dy = x + np.array([[0.0, eps, 0.0]])
    dz = x + np.array([[0.0, 0.0, eps]])
    hess = np.zeros((x.shape[0], 3, 3), float)
    gpnt = mol.obasis.compute_grid_gradient_dm(dm, x)
    hess[:,:,0] = (mol.obasis.compute_grid_gradient_dm(dm, dx) - gpnt) / eps
    hess[:,:,1] = (mol.obasis.compute_grid_gradient_dm(dm, dy) - gpnt) / eps
    hess[:,:,2] = (mol.obasis.compute_grid_gradient_dm(dm, dz) - gpnt) / eps

    return 0.5*(hess + np.transpose(hess, axes=(0, 2, 1)))
