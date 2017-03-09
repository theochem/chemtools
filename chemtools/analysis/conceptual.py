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
'''Module for Conceptual Density Functional Theory Analysis of Quantum Chemistry Output Files.

   This modules contains wrappers which take outputs of quantum chemistry softwares and
   compute various conceptual density functional theory (DFT) descriptive tools.
'''

import numpy as np
from horton import log
from horton import IOData, BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes
from chemtools.tool.globaltool import LinearGlobalTool, QuadraticGlobalTool, ExponentialGlobalTool, RationalGlobalTool, GeneralGlobalTool
from chemtools.tool.localtool import LinearLocalTool, QuadraticLocalTool
from chemtools.tool.condensedtool import CondensedTool
from chemtools.utils import CubeGen



class GlobalConceptualDFT(object):
    '''
    Class for global conceptual density functional theory (DFT) analysis of molecular quantum
    chemistry output file(s). If only one molecule is provided, the frontiner molecular orbital
    (FMO) approach is invoked, otherwise finite difference (FD) approach is taken.
    '''
    def __init__(self, dict_values, model='quadratic', **kwargs):
        '''
        Parameters
        ----------
        dict_values : dict
            Dictionary of number_electron:energy
        model : str, default='quadratic'
            Energy model used to calculate global descriptive tools.
            The available models include: 'linear', 'quadratic', 'exponential', 'rational',
            and 'general'. Please see '' for more information.
        '''
        # available models for global tools
        select_tool = {'linear': LinearGlobalTool, 'quadratic': QuadraticGlobalTool,
                       'exponential': ExponentialGlobalTool, 'rational': RationalGlobalTool,
                       'general': GeneralGlobalTool,
                      }
        if model.lower() not in select_tool.keys():
            raise ValueError('Model={0} is not available!'.format(model.lower()))
        self._model = model.lower()

        if self._model != 'general':
            # get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            # get energy values of sorted number of electrons
            energies = [dict_values[key] for key in nelectrons]

            if len(nelectrons) != 3:
                raise NotImplementedError('For model={0}, three energy values are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons are required! #electrons={1}'.format(self._model, nelectrons))

            # check consecutive number of elctrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of electrons should differ by 1! #electrons={1}'.format(self._model, nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for global tool
            args = [energies[1], energies[2], energies[0], reference]
            self._tool = select_tool[model](*args)
        else:
            self._tool = select_tool[model](*args, **kwargs)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    @property
    def model(self):
        '''Energy model.'''
        return self._model

    def __getattr__(self, attr):
        '''
        Return class attribute.
        '''
        if not hasattr(self._tool, attr):
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    def __repr__(self):
        '''
        '''
        available = dir(self._tool)
        attrs = [attr for attr in available if not callable(getattr(self, attr)) and not attr.startswith('_')]
        attrs.sort()
        methods = [attr for attr in available if callable(getattr(self, attr)) and not attr.startswith('_')]
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 40)
        content += '\n'.join(attrs)
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 40)
        content += '\n'.join(methods)
        return content

    def _log_init(self):
        '''
        Print an initial informative message.
        '''
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference Energy', self._tool.energy_zero),
                         ('Reference #Electrons', self._tool.n0),
                         ('Energy Values', [getattr(self._tool, attr) for attr in ['_energy_minus', 'energy_zero', '_energy_plus']]),
                         #('Parameters', self._tool.energy_params)
                        ])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, **kwargs):
        '''
        Initialize class from files.
        '''
        # case of one file not given as a list
        if isinstance(filenames, (str, unicode)):
            return cls.from_iodata(IOData.from_file(filenames), model, **kwargs)
        iodatas = []
        # case of list of file(s)
        for filename in filenames:
            iodatas.append(IOData.from_file(filename))
        return cls.from_iodata(iodatas, model, **kwargs)

    @classmethod
    def from_iodata(cls, iodatas, model, **kwargs):
        '''
        Initialize class from `IOData` objects.
        '''
        # case of one IOData object not given as a list
        if isinstance(iodatas, IOData):
            iodatas = [iodatas]

        if len(iodatas) == 1:
            # Frontiner Molecular Orbital (FMO) Approach
            mol = iodatas[0]
            # get homo & lumo energy of alpha electrons
            homo_energy = mol.exp_alpha.homo_energy
            lumo_energy = mol.exp_alpha.lumo_energy
            # get number of alpha electrons
            nelec = int(np.sum(mol.exp_alpha.occupations))
            if hasattr(mol, 'exp_beta') and mol.exp_beta is not None:
                # add number of beta electrons
                nelec += int(np.sum(mol.exp_beta.occupations))
                # use homo energy of beta electrons, if it has higher energy
                if mol.exp_beta.homo_energy > homo_energy:
                    homo_energy = mol.exp_beta.homo_energy
                # use lumo energy of beta electrons, if it has lower energy
                if mol.exp_beta.lumo_energy < lumo_energy:
                    lumo_energy = mol.exp_beta.lumo_energy
            else:
                nelec *= 2
            if not hasattr(mol, 'energy'):
                raise ValueError('Molecule object does not contain energy value!')
            # store number of electron and energy in a dictionary
            dict_values = dict([(nelec, mol.energy),
                                (nelec + 1, mol.energy + lumo_energy),
                                (nelec - 1, mol.energy - homo_energy)])
        else:
            # Finite Difference (FD) Approach
            dict_values = {}
            for iodata in iodatas:
                # get the number of electrons
                nelec = int(np.sum(iodata.exp_alpha.occupations))
                if hasattr(iodata, 'exp_beta') and iodata.exp_beta is not None:
                    nelec += int(np.sum(iodata.exp_beta.occupations))
                else:
                    nelec *= 2
                if nelec in dict_values.keys():
                    raise ValueError('Two molecules have {0} electrons!'.format(nelec))
                # get the energy
                if not hasattr(iodata, 'energy'):
                    raise ValueError('Molecule object does not contain energy value!')
                # store number of electron and energy in a dictionary
                dict_values[nelec] = iodata.energy

        return cls(dict_values, model, **kwargs)


class LocalConceptualDFT(object):
    '''
    Class for local conceptual density functional theory (DFT) analysis of molecular quantum
    chemistry output file(s). If only one molecule is provided, the frontiner molecular orbital
    (FMO) approach is invoked, otherwise finite difference (FD) approach is taken.
    '''
    def __init__(self, dict_values, model='quadratic', **kwargs):
        '''
        Parameters
        ----------
        dict_value : dict
            Dictionary of number_electron:density
        model : str, default='quadratic'
            Energy model used to calculate local descriptive tools.
            The available models include: 'linear', 'quadratic', and 'general'.
            Please see '' for more information.
        '''
        # available models for local tools
        select_tool = {'linear':LinearLocalTool, 'quadratic':QuadraticLocalTool}

        if model not in select_tool.keys():
            raise ValueError('Model={0} is not available!'.format(model))
        self._model = model

        if self._model != 'general':
            # get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            # get energy values of sorted number of electrons (from small to big)
            densities = [dict_values[key] for key in nelectrons]

            if len(nelectrons) != 3:
                raise NotImplementedError('For model={0}, three densities are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons are required! #electrons={1}'.format(self._model, nelectrons))

            # check consecutive number of elctrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of electrons should differ by 1! #electrons={1}'.format(self._model, nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for local tool
            args = [densities[1], densities[2], densities[0], reference]
            self._tool = select_tool[model](*args)
        else:
            self._tool = select_tool[model](*args)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    def __getattr__(self, attr):
        '''
        Return class attribute.
        '''
        if not hasattr(self._tool, attr):
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    def _log_init(self):
        '''
        Print an initial informative message.
        '''
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference #Electrons', self._tool.n0),
                        ])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, points=None, **kwargs):
        '''
        Initialize class from files.
        '''
        # case of one file not given as a list
        if isinstance(filenames, (str, unicode)):
            return cls.from_iodata(IOData.from_file(filenames), model, points, **kwargs)
        # case of list of file(s)
        iodatas = []
        for filename in filenames:
            iodatas.append(IOData.from_file(filename))
        return cls.from_iodata(iodatas, model, points, **kwargs)

    @classmethod
    def from_iodata(cls, iodatas, model, points=None, **kwargs):
        '''
        Initialize class from `IOData` objects.
        '''
        # case of one IOData object not given as a list
        if isinstance(iodatas, IOData):
            iodatas = [iodatas]

        if len(iodatas) == 1:
            # Frontiner Molecular Orbital (FMO) Approach
            mol = iodatas[0]
            # get homo & lumo index & energy of alpha electrons
            homo_index = mol.exp_alpha.get_homo_index()
            lumo_index = mol.exp_alpha.get_lumo_index()
            homo_energy = mol.exp_alpha.homo_energy
            lumo_energy = mol.exp_alpha.lumo_energy
            # set homo & lumo expasion to alpha orbitals
            homo_exp, lumo_exp = 'exp_alpha', 'exp_alpha'

            # get number of alpha electrons
            nelec = int(np.sum(mol.exp_alpha.occupations))
            if hasattr(mol, 'exp_beta') and mol.exp_beta is not None:
                # add number of beta electrons
                nelec += int(np.sum(mol.exp_beta.occupations))
                # use homo energy of beta electrons, if it has higher energy
                if mol.exp_beta.homo_energy > homo_energy:
                    homo_energy = mol.exp_beta.homo_energy
                    homo_index = mol.exp_beta.get_homo_index()
                    homo_exp = 'exp_beta'
                # use lumo energy of beta electrons, if it has lower energy
                if mol.exp_beta.lumo_energy < lumo_energy:
                    lumo_energy = mol.exp_beta.lumo_energy
                    lumo_index = mol.exp_beta.get_lumo_index()
                    lumo_exp = 'exp_beta'
            else:
                nelec *= 2
            # compute homo & lumo density
            homo_dens = mol.obasis.compute_grid_orbitals_exp(getattr(mol, homo_exp), points,
                                                             np.array([homo_index]))**2
            lumo_dens = mol.obasis.compute_grid_orbitals_exp(getattr(mol, lumo_exp), points,
                                                             np.array([lumo_index]))**2
            homo_dens = homo_dens.flatten()
            lumo_dens = lumo_dens.flatten()
            # store number of electron and density in a dictionary
            dens = mol.obasis.compute_grid_density_dm(mol.get_dm_full(), points)
            dict_values = dict([(nelec, dens),
                                (nelec + 1, dens + lumo_dens),
                                (nelec - 1, dens - homo_dens)])
        else:
            # Finite Difference (FD) Approach
            dict_values = {}
            for iodata in iodatas:
                # get number of electrons
                nelec = int(np.sum(iodata.exp_alpha.occupations))
                if hasattr(iodata, 'exp_beta') and iodata.exp_beta is not None:
                    nelec += int(np.sum(iodata.exp_beta.occupations))
                else:
                    nelec *= 2
                if nelec in dict_values.keys():
                    raise ValueError('Two molecules have {0} electrons!'.format(nelec))
                # compute density
                dens = iodata.obasis.compute_grid_density_dm(iodata.get_dm_full(), points)
                # store number of electron and density in a dictionary
                dict_values[nelec] = dens

        return cls(dict_values, model, **kwargs)
