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
        if self._mol.exp_beta is not None:
            n_elec += int(np.sum(self._mol.exp_beta.occupations))
            if self._mol.exp_beta.homo_energy > homo_energy:
                homo_energy = self._mol.exp_beta.homo_energy
            if self._mol.exp_beta.lumo_energy < lumo_energy:
                lumo_energy = self._mol.exp_beta.lumo_energy

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

    # def make_scripts_nci(self):
    #     '''
    #     Genrate scripts to plot NCI through VMD.
    #     '''
    #     # uses self._localtool.compute_nci() and then write scripts.
    #     pass
