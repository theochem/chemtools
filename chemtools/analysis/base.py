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

    def compute_nci(self, filename, isosurf=0.50000, denscut=0.10, cube=None):
        '''
        Generate all files to plot NCI.
        '''
        # FIXME: I still have to tweak the cutoffs
        #
        # the files used:
        #
        densfile = filename + '-dens.cube'    # name of the density cube file
        rdgfile  = filename + '-grad.cube'    # name of the reduced density gradient cube file
        vmdfile  = filename + '.vmd'          # name of the reduced vmd script file

        if cube is None:
            cube = CubeGen.from_molecule(self._mol.numbers, self._mol.pseudo_numbers, self._mol.coordinates, spacing=0.1,threshold=2.0)

        dm_full = self._mol.get_dm_full()

        dens = self._mol.obasis.compute_grid_density_dm(dm_full, cube.gridpoints)
        grad = self._mol.obasis.compute_grid_gradient_dm(dm_full, cube.gridpoints)
        dloctool = DensityLocalTool(dens, grad)
        rdgrad = dloctool.reduced_density_gradient

        #mask rdgrad:
        rdgrad[abs(dens) > denscut] = 100.0

        hess = self._compute_hessian(self._mol, cube.gridpoints)
        dens = np.sign(np.linalg.eigvalsh(hess)[:,1])*dens

        cube.dump_cube(densfile ,dens*100.0)
        cube.dump_cube(rdgfile,rdgrad)
        _print_vmd(vmdfile, densfile, rdgfile, isosurf, denscut*100.0)

        plt.xlim(-0.1, 0.1)
        plt.ylim(0.0, 2.0)
        plt.scatter(dens,rdgrad)
        plt.show()

    # ???: ok to put it here?

def _compute_hessian(mol, x):
    #untill the electronic Hessian is included in HORTON, we do it by finite difference:
    dm = mol.get_dm_full()
    eps = 1.0e-5
    dx = x + np.array([[eps, 0.0, 0.0]])
    dy = x + np.array([[0.0, eps, 0.0]])
    dz = x + np.array([[0.0, 0.0, eps]])
    hess = np.zeros((x.shape[0], 3, 3), float)
    gpnt = mol.obasis._compute_gradient(dm, x)
    hess[:,:,0] = (mol.obasis._compute_gradient(dm, dx) - gpnt) / eps
    hess[:,:,1] = (mol.obasis._compute_gradient(dm, dy) - gpnt) / eps
    hess[:,:,2] = (mol.obasis._compute_gradient(dm, dz) - gpnt) / eps

    return 0.5*(hess + np.transpose(hess, axes=(0, 2, 1)))


def _print_vmd(file, densfile, rdgfile, isosurf, denscut):
    # print out the .vmd file...
    # FIXEME: The latest version of NCIPLOT uses a slightly different script. I have to see what changed and why
    with open(file, 'w') as f:
        print >> f, '#!/usr/local/bin/vmd'
        print >> f, '# VMD script written by save_state $Revision: 1.10 $'
        print >> f, '# VMD version: 1.8.6              '
        print >> f, 'set viewplist                     '
        print >> f, 'set fixedlist                     '
        print >> f, '# Display settings                '
        print >> f, 'display projection   Orthographic '
        print >> f, 'display nearclip set 0.000000     '
        print >> f, '# load new molecule               '
        print >> f, 'mol new {0}  type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'.format(densfile)
        print >> f, 'mol addfile {0}  type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all'.format(rdgfile)
        print >> f, '#'
        print >> f, '# representation of the atoms'
        print >> f, 'mol delrep 0 top'
        print >> f, 'mol representation CPK 1.000000 0.300000 118.000000 131.000000'
        print >> f, 'mol color Name'
        print >> f, 'mol selection {{all}}'
        print >> f, 'mol material Opaque'
        print >> f, 'mol addrep top'
        print >> f, '#'
        print >> f, '# add representation of the surface'
        print >> f, 'mol representation Isosurface {0:.5f} 1 0 0 1 1'.format(isosurf)
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
