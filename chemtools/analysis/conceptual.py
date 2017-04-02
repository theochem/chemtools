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
"""Module for Conceptual Density Functional Theory Analysis of Quantum Chemistry Output Files.

   This modules contains wrappers which take outputs of quantum chemistry software and
   compute various conceptual density functional theory (DFT) descriptive tools.
"""

import numpy as np
from horton import log
from horton import IOData, BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes
from chemtools.tool.globaltool import LinearGlobalTool, QuadraticGlobalTool
from chemtools.tool.globaltool import ExponentialGlobalTool, RationalGlobalTool, GeneralGlobalTool
from chemtools.tool.localtool import LinearLocalTool, QuadraticLocalTool
from chemtools.utils import CubeGen


class GlobalConceptualDFT(object):
    """
    Class for global conceptual density functional theory (DFT) analysis of molecular quantum
    chemistry output file(s). If only one molecule is provided, the frontier molecular orbital
    (FMO) approach is invoked, otherwise finite difference (FD) approach is taken.
    """
    def __init__(self, dict_values, model='quadratic'):
        """
        Parameters
        ----------
        dict_values : dict
            Dictionary of number_electron:energy
        model : str, default='quadratic'
            Energy model used to calculate global descriptive tools.
            The available models include: 'linear', 'quadratic', 'exponential', 'rational',
            and 'general'. Please see '' for more information.
        """
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
                raise NotImplementedError(
                    'For model={0}, three energy values are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons '.format(self._model) +
                                 'are required! #electrons={0}'.format(nelectrons))

            # check consecutive number of electrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of'.format(self._model) +
                                 'electrons should differ by 1! #electrons={0}'.format(nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for global tool
            args = [energies[1], energies[2], energies[0], reference]
            self._tool = select_tool[self._model](*args)
        else:
            # self._tool = select_tool[model](*args, **kwargs)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    @property
    def model(self):
        """Energy model."""
        return self._model

    def __getattr__(self, attr):
        """
        Return class attribute.
        """
        # if not hasattr(self._tool, attr):
        value = getattr(self._tool, attr, 'error')
        if value == 'error':
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    def __repr__(self):
        """
        Print table of available class attributes and methods.
        """
        available = dir(self._tool)
        is_public = lambda item: not item.startswith('_')
        attrs = [atr for atr in available if not callable(getattr(self, atr)) and is_public(atr)]
        attrs.sort()
        methods = [attr for attr in available if callable(getattr(self, attr)) and is_public(attr)]
        methods.sort()
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        for attr in attrs:
            if getattr(self._tool, attr) is not None:
                content += '\n%s   % .6f' % (attr.ljust(25), getattr(self._tool, attr))
            else:
                content += '\n%s   %s' % (attr.ljust(25), ' ---')
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """
        Print an initial informative message.
        """
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference Energy', self._tool.energy_zero),
                         ('Reference #Electrons', self._tool.n0),
                         ('Energy Values', [getattr(self._tool, attr) for attr in
                                            ['_energy_minus', 'energy_zero', '_energy_plus']]),
                         # ('Parameters', self._tool.params)
                        ])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, **kwargs):
        """
        Initialize class from files.

        Parameters
        ----------
        filenames : list,tuple
            List/tuple containing the path to molecule's files.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        points : np.array
            Points on which the local descriptive tools is evaluated.
        """
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
        """
        Initialize class from `IOData` objects.

        Parameters
        ----------
        iodatas : list,tuple
            List/tuple containing `IOData` objects.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        points : np.array
            Points on which the local descriptive tools is evaluated.
        """
        # case of one IOData object not given as a list
        if isinstance(iodatas, IOData):
            iodatas = [iodatas]

        if len(iodatas) == 1:
            # Frontier Molecular Orbital (FMO) Approach
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
            for index, iodata in enumerate(iodatas):
                # check atomic numbers
                if index == 0:
                    atomic_numbers = iodata.numbers
                elif not np.all(abs(atomic_numbers - iodata.numbers) < 1.e-6):
                    raise ValueError(
                        'Molecule 1 & {0} have different atomic numbers!'.format(index + 1))
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
    """
    Class for local conceptual density functional theory (DFT) analysis of molecular quantum
    chemistry output file(s). If only one molecule is provided, the frontier molecular orbital
    (FMO) approach is invoked, otherwise finite difference (FD) approach is taken.

    Note: If FD approach is taken, the geometries of all molecules need to be the same.
    """
    def __init__(self, dict_values, model='quadratic', coordinates=None, numbers=None):
        """
        Parameters
        ----------
        dict_value : dict
            Dictionary of number_electron:density
        model : str, default='quadratic'
            Energy model used to calculate local descriptive tools.
            The available models include: 'linear', 'quadratic', and 'general'.
            Please see '' for more information.
        """
        # available models for local tools
        select_tool = {'linear': LinearLocalTool, 'quadratic': QuadraticLocalTool}

        if model.lower() not in select_tool.keys():
            raise ValueError('Model={0} is not available!'.format(model))
        self._model = model.lower()

        # check shape of coordinates
        if coordinates is not None and coordinates.shape[1] != 3:
            raise ValueError('Argument coordinate should be a 2D-array with 3 columns!' +
                             ' shape={0}'.format(coordinates.shape))

        # check number of atoms given by numbers and coordinates match
        if len(numbers) != len(coordinates):
            raise ValueError('Numbers & coordinates should represent same number of ' +
                             'atoms! {0}!={1}'.format(len(numbers), len(coordinates)))

        self._coordinates = coordinates
        self._numbers = numbers

        # check densities shape
        for index, value in enumerate(dict_values.values()):
            if index == 0:
                shape = value.shape
            elif value.shape != shape:
                raise ValueError(
                    'Densities should have the same shape! {1}!={2}'.format(shape, value.shape))

        if self._model != 'general':
            # get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            # get energy values of sorted number of electrons (from small to big)
            densities = [dict_values[key] for key in nelectrons]

            if len(nelectrons) != 3:
                raise NotImplementedError(
                    'For model={0}, three densities are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons '.format(self._model) +
                                 'are required! #electrons={0}'.format(nelectrons))

            # check consecutive number of electrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of'.format(self._model) +
                                 'electrons should differ by 1! #electrons={0}'.format(nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for local tool
            args = [densities[1], densities[2], densities[0], reference]
            self._tool = select_tool[self._model](*args)
        else:
            # self._tool = select_tool[model](*args)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    @property
    def model(self):
        """Energy model."""
        return self._model

    @property
    def n0(self):
        """
        Reference number of electrons, i.e. :math:`N_0`.
        """
        return self._tool.n0

    @property
    def coordinates(self):
        """Cartesian coordinates of atoms."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic numbers."""
        return self._numbers

    def __getattr__(self, attr):
        """
        Return class attribute.
        """
        # if not hasattr(self._tool, attr):
        value = getattr(self._tool, attr, None)
        if value is None:
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    def __repr__(self):
        """
        Print table of available class attributes and methods.
        """
        available = dir(self._tool) + dir(self)
        is_public = lambda item: not item.startswith('_')
        attrs = [atr for atr in available if not callable(getattr(self, atr)) and is_public(atr)]
        attrs.sort()
        # remove n0, because it is both an attr of self._tool and self (duplicate)
        attrs.remove('n0')
        methods = [attr for attr in available if callable(getattr(self, attr)) and is_public(attr)]
        methods.sort()
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(attrs)
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """
        Print an initial informative message.
        """
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference #Electrons', self._tool.n0),
                        ])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, points):
        """
        Initialize class from files.

        Parameters
        ----------
        filenames : list,tuple
            List/tuple containing the path to molecule's files.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        points : np.array
            Points on which the local descriptive tools is evaluated.
        """
        # case of one file not given as a list
        if isinstance(filenames, (str, unicode)):
            return cls.from_iodata(IOData.from_file(filenames), model, points)
        # case of list of file(s)
        iodatas = []
        for filename in filenames:
            iodatas.append(IOData.from_file(filename))
        return cls.from_iodata(iodatas, model, points)

    @classmethod
    def from_iodata(cls, iodatas, model, points):
        """
        Initialize class from `IOData` objects.

        Parameters
        ----------
        iodatas : list,tuple
            List/tuple containing `IOData` objects.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        points : np.array
            Points on which the local descriptive tools is evaluated.
        """
        # check points array
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError('Argument points should be a 2D array with 3 columns! ' +
                             'given shape={0}'.format(points.shape))

        # case of one IOData object not given as a list
        if isinstance(iodatas, IOData):
            iodatas = [iodatas]

        if len(iodatas) == 1:
            # Frontier Molecular Orbital (FMO) Approach
            mol = iodatas[0]
            # get homo & lumo index & energy of alpha electrons
            homo_index = mol.exp_alpha.get_homo_index()
            lumo_index = mol.exp_alpha.get_lumo_index()
            homo_energy = mol.exp_alpha.homo_energy
            lumo_energy = mol.exp_alpha.lumo_energy
            # set homo & lumo expansion to alpha orbitals
            homo_exp, lumo_exp = 'exp_alpha', 'exp_alpha'

            # get number of alpha electrons
            nelec = int(np.sum(mol.exp_alpha.occupations))
            if hasattr(mol, 'exp_beta') and mol.exp_beta is not None:
                # add number of beta electrons
                nelec += int(np.sum(mol.exp_beta.occupations))
                # use homo energy of beta electrons, if it has higher energy
                if mol.exp_beta.homo_energy > homo_energy:
                    # homo_energy = mol.exp_beta.homo_energy
                    homo_index = mol.exp_beta.get_homo_index()
                    homo_exp = 'exp_beta'
                # use lumo energy of beta electrons, if it has lower energy
                if mol.exp_beta.lumo_energy < lumo_energy:
                    # lumo_energy = mol.exp_beta.lumo_energy
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
            # get coordinates of molecule
            coordinates = mol.coordinates
            # get atomic numbers of molecule
            atomic_numbers = mol.numbers

        else:
            # Finite Difference (FD) Approach
            dict_values = {}
            for index, iodata in enumerate(iodatas):
                # check atomic numbers
                if index == 0:
                    atomic_numbers = iodata.numbers
                elif not np.all(abs(atomic_numbers - iodata.numbers) < 1.e-6):
                    raise ValueError(
                        'Molecule 1 & {0} have different atomic numbers!'.format(index + 1))
                # check coordinates
                if index == 0:
                    # get coordinates of molecule
                    coordinates = iodata.coordinates
                elif not np.all(abs(coordinates - iodata.coordinates) < 1.e-4):
                    raise ValueError(
                        'Molecule 1 & {0} have different geometries!'.format(index + 1))
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

        return cls(dict_values, model, coordinates, atomic_numbers)


class CondensedConceptualDFT(object):
    """
    Class for condensed conceptual density functional theory (DFT) analysis of molecular quantum
    chemistry output file(s). If only one molecule is provided, the frontier molecular orbital
    (FMO) approach is invoked, otherwise finite difference (FD) approach is taken.

    Note: If FD approach is taken, the geometries of molecules can be different only for RMF
          approach.
    """
    def __init__(self, dict_values, model='quadratic'):
        """
        Parameters
        ----------
        dict_value : dict
            Dictionary of number_electron:atomic_populations
        model : str, default='quadratic'
            Energy model used to calculate local descriptive tools.
            The available models include: 'linear', 'quadratic', and 'general'.
            Please see '' for more information.
        """
        # available models for local tools
        select_tool = {'linear': LinearLocalTool, 'quadratic': QuadraticLocalTool}

        if model.lower() not in select_tool.keys():
            raise ValueError('Model={0} is not available!'.format(model))
        self._model = model.lower()

        if self._model != 'general':
            # get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            # get energy values of sorted number of electrons (from small to big)
            populations = [dict_values[key] for key in nelectrons]

            if len(nelectrons) != 3:
                raise NotImplementedError(
                    'For model={0}, three densities are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons'.format(self._model) +
                                 ' are required! #electrons={0}'.format(nelectrons))

            # check consecutive number of electrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of '.format(self._model) +
                                 'electrons should differ by 1! #electrons={0}'.format(nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for condense tool
            args = [populations[1], populations[2], populations[0], reference]
            self._tool = select_tool[self._model](*args)
        else:
            # self._tool = select_tool[model](*args)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    def __getattr__(self, attr):
        """
        Return class attribute.
        """
        # if not hasattr(self._tool, attr):
        value = getattr(self._tool, attr, None)
        if value is None:
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    def __repr__(self):
        """
        Print table of available class attributes and methods.
        """
        available = dir(self._tool)
        is_public = lambda item: not item.startswith('_')
        attrs = [atr for atr in available if not callable(getattr(self, atr)) and is_public(atr)]
        attrs.sort()
        attrs.remove('n0')
        methods = [attr for attr in available if callable(getattr(self, attr)) and is_public(attr)]
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(attrs)
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """
        Print an initial informative message.
        """
        if log.do_medium:
            log('Initialize: Condensed Class')
            # log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference #Electrons', self._tool.n0),
                        ])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, approach='FMR', scheme='h', grid=None, proatomdb=None):
        """
        Initialize class from files.

        Parameters
        ----------
        filenames : list,tuple
            List/tuple containing the path to molecule's files.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        approach : str
            Choose between 'FMR' (fragment of molecular response) or 'RMF'
            (response of molecular fragment).
        scheme: str
            Partitioning scheme. Options: 'h', 'hi', 'mbis'.
        grid: instance of ``BeckeMolGrid``
            Grid used for partitioning
        proatomdb: instance of ``ProAtomDB``
            Proatom database used for partitioning. Only 'h' and 'hi' requires that.
        """
        # case of one file not given as a list
        if isinstance(filenames, (str, unicode)):
            return cls.from_iodata(IOData.from_file(filenames), model, grid, approach,
                                   scheme, proatomdb)
        # case of list of file(s)
        iodatas = []
        for filename in filenames:
            iodatas.append(IOData.from_file(filename))
        return cls.from_iodata(iodatas, model, grid, approach, scheme, proatomdb)

    @classmethod
    def from_iodata(cls, iodatas, model, grid=None, approach='FMR', scheme='mbis',
                    proatomdb=None, **kwargs):
        """
        Initialize class from `IOData` objects.

        Parameters
        ----------
        iodatas : list,tuple
            List/tuple containing `IOData` instances.
        model : str
            Energy model used to calculate local descriptive tools.
            Available models are 'linear' and 'quadratic'.
        approach : str
            Choose between 'FMR' (fragment of molecular response) or 'RMF'
            (response of molecular fragment).
        scheme: str
            Partitioning scheme. Options: 'h', 'hi', 'mbis'.
        grid: instance of ``BeckeMolGrid``
            Grid used for partitioning
        proatomdb: instance of ``ProAtomDB``
            Proatom database used for partitioning. Only 'h' and 'hi' requires that.
        """
        if grid is not None:
            if not isinstance(grid, BeckeMolGrid):
                raise ValueError('Currently, only BeckeMolGrid is supported for condensing!')

        # case of one IOData object not given as a list
        if isinstance(iodatas, IOData):
            iodatas = [iodatas]

        if len(iodatas) == 1:
            # Frontier Molecular Orbital (FMO) Approach
            mol = iodatas[0]
            if grid is None:
                grid = BeckeMolGrid(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                    agspec='fine', random_rotate=False, mode='keep')
            points = grid.points
            # get homo & lumo index & energy of alpha electrons
            homo_index = mol.exp_alpha.get_homo_index()
            lumo_index = mol.exp_alpha.get_lumo_index()
            homo_energy = mol.exp_alpha.homo_energy
            lumo_energy = mol.exp_alpha.lumo_energy
            # set homo & lumo expansion to alpha orbitals
            homo_exp, lumo_exp = 'exp_alpha', 'exp_alpha'

            # get number of alpha electrons
            nelec = int(np.sum(mol.exp_alpha.occupations))
            if hasattr(mol, 'exp_beta') and mol.exp_beta is not None:
                # add number of beta electrons
                nelec += int(np.sum(mol.exp_beta.occupations))
                # use homo energy of beta electrons, if it has higher energy
                if mol.exp_beta.homo_energy > homo_energy:
                    # homo_energy = mol.exp_beta.homo_energy
                    homo_index = mol.exp_beta.get_homo_index()
                    homo_exp = 'exp_beta'
                # use lumo energy of beta electrons, if it has lower energy
                if mol.exp_beta.lumo_energy < lumo_energy:
                    # lumo_energy = mol.exp_beta.lumo_energy
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

            # compute density of reference molecule
            dens = mol.obasis.compute_grid_density_dm(mol.get_dm_full(), points)
            # partitioning
            if scheme.lower() not in wpart_schemes:
                raise ValueError('Partitioning scheme={0} not supported! ' +
                                 'Select from: {1}'.format(scheme, wpart_schemes))

            wpart = wpart_schemes[scheme]
            if isinstance(grid, BeckeMolGrid):
                # compute population of reference system
                if scheme.lower() not in ['mbis', 'b']:
                    if proatomdb is None:
                        proatomdb = ProAtomDB.from_refatoms(mol.numbers)
                    kwargs = {'proatomdb': proatomdb}
                else:
                    kwargs = {}
                part0 = wpart(mol.coordinates, mol.numbers,
                              mol.pseudo_numbers, grid, dens, **kwargs)
                part0.do_all()
                pops0 = part0['populations']
                # compute population of N+1 and N-1 electron system
                if approach.lower() == 'fmr':
                    # fragment of molecular response
                    popsp = cls.condense_to_atoms(dens + lumo_dens, part0)
                    popsm = cls.condense_to_atoms(dens - homo_dens, part0)
                elif approach.lower() == 'rmf':
                    # response of molecular fragment
                    partp = wpart(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                  grid, dens + lumo_dens, **kwargs)
                    partp.do_all()
                    popsp = partp['populations']
                    partm = wpart(mol.coordinates, mol.numbers, mol.pseudo_numbers,
                                  grid, dens - homo_dens, **kwargs)
                    partm.do_all()
                    popsm = partm['populations']
            elif isinstance(grid, CubeGen):
                raise NotImplementedError(
                    'Condensing with CubGen will be implemented in near future!')
            else:
                raise NotImplementedError('Condensing does not support grid={0}!'.format(grid))

            # Store number of electron and populations in a dictionary
            dict_values = dict([(nelec, pops0),
                                (nelec + 1, popsp),
                                (nelec - 1, popsm)])
        else:
            # Finite Difference (FD) Approach
            dict_values = {}
            same_coordinates = True
            for index, iodata in enumerate(iodatas):
                # check atomic numbers
                if index == 0:
                    atomic_numbers = iodata.numbers
                elif not np.all(abs(atomic_numbers - iodata.numbers) < 1.e-6):
                    raise ValueError(
                        'Molecule 1 & {0} have different atomic numbers!'.format(index + 1))
                # check coordinates of grid and molecule match
                if grid is not None:
                    if not np.all(abs(grid.centers - iodata.coordinates) < 1.e-4):
                        print abs(grid.centers - iodata.coordinates)
                        raise ValueError(
                            'Coordinates of grid and molecule {0} should match!'.format(index))
                if index == 0:
                    coordinates = iodata.coordinates
                elif not np.all(abs(coordinates - iodata.coordinates) < 1.e-4):
                    print coordinates - iodata.coordinates
                    same_coordinates = False
                    # if geometries are not the same, only rmf can be done.
                    if approach.lower() == 'fmr':
                        raise ValueError('When geometries of molecules are different, ' +
                                         'only approach=RMF is possible!')

                # get number of electrons
                nelec = int(np.sum(iodata.exp_alpha.occupations))
                if hasattr(iodata, 'exp_beta') and iodata.exp_beta is not None:
                    nelec += int(np.sum(iodata.exp_beta.occupations))
                else:
                    nelec *= 2
                if nelec in dict_values.keys():
                    raise ValueError('Two molecules have {0} electrons!'.format(nelec))
                # store number of electron and iodata in a dictionary
                dict_values[nelec] = iodata

            # Get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            if len(nelectrons) != 3:
                raise ValueError('Condensed conceptual DFT within FD approach, ' +
                                 'currently only works for 3 molecules!')

            reference = nelectrons[1]
            # get energy values of sorted number of electrons (from small to big)
            molm, mol0, molp = [dict_values[key] for key in nelectrons]

            # compute density of reference molecule
            if grid is None:
                grid = BeckeMolGrid(mol0.coordinates, mol0.numbers, mol0.pseudo_numbers,
                                    agspec='fine', random_rotate=False, mode='keep')
            points = grid.points
            dens0 = mol0.obasis.compute_grid_density_dm(mol0.get_dm_full(), points)

            # partitioning
            if scheme.lower() not in wpart_schemes:
                raise ValueError('Partitioning scheme={0} not supported!'.format(scheme) +
                                 'Select from: {0}'.format(wpart_schemes))

            wpart = wpart_schemes[scheme]
            if isinstance(grid, BeckeMolGrid):
                # compute population of reference system
                if scheme.lower() not in ['mbis', 'b']:
                    if proatomdb is None:
                        proatomdb = ProAtomDB.from_refatoms(mol0.numbers)
                    kwargs = {'proatomdb': proatomdb}
                else:
                    kwargs = {}
                part0 = wpart(mol0.coordinates, mol0.numbers,
                              mol0.pseudo_numbers, grid, dens0, **kwargs)
                part0.do_all()
                pops0 = part0['populations']
                # compute population of N+1 and N-1 electron system
                if approach.lower() == 'fmr':
                    # fragment of molecular response
                    assert same_coordinates
                    densp = molp.obasis.compute_grid_density_dm(molp.get_dm_full(), points)
                    densm = molm.obasis.compute_grid_density_dm(molm.get_dm_full(), points)
                    popsp = cls.condense_to_atoms(densp, part0)
                    popsm = cls.condense_to_atoms(densm, part0)
                elif approach.lower() == 'rmf':
                    # response of molecular fragment
                    if not same_coordinates:
                        grid = BeckeMolGrid(molp.coordinates, molp.numbers, molp.pseudo_numbers,
                                            agspec='fine', random_rotate=False, mode='keep')
                        points = grid.points
                    densp = molp.obasis.compute_grid_density_dm(molp.get_dm_full(), points)
                    partp = wpart(molp.coordinates, molp.numbers,
                                  molp.pseudo_numbers, grid, densp, **kwargs)
                    partp.do_all()
                    popsp = partp['populations']
                    if not same_coordinates:
                        grid = BeckeMolGrid(molm.coordinates, molm.numbers, molm.pseudo_numbers,
                                            agspec='fine', random_rotate=False, mode='keep')
                        points = grid.points
                    densm = molm.obasis.compute_grid_density_dm(molm.get_dm_full(), points)
                    partm = wpart(molm.coordinates, molm.numbers,
                                  molm.pseudo_numbers, grid, densm, **kwargs)
                    partm.do_all()
                    popsm = partm['populations']
            elif isinstance(grid, CubeGen):
                raise NotImplementedError(
                    'Condensing with CubGen will be implemented in near future!')
            else:
                raise NotImplementedError('Condensing does not support grid={0}!'.format(grid))

            # Store number of electron and populations in a dictionary
            dict_values = dict([(reference, pops0),
                                (reference + 1, popsp),
                                (reference - 1, popsm)])
        return cls(dict_values, model)

    @staticmethod
    def condense_to_atoms(local_property, part):
        r"""
        Return condensed values of the local descriptor
        :math:`p_{\text{local}}\left(\mathbf{r}\right)` into atomic contribution :math:`P_A`
        defined as,

        .. math::
           P_A = \int \omega_A\left(\mathbf{r}\right)
                      p_{\text{local}}\left(\mathbf{r}\right) d\mathbf{r}

        Parameters
        ----------
        local_property : np.ndarray
            Local descriptor evaluated on grid.
        part : part instance
            Instance of `HORTON` partitioning calss.
        """
        condensed = np.zeros(part.natom)
        for index in xrange(part.natom):
            at_grid = part.get_grid(index)
            at_weight = part.cache.load('at_weights', index)
            wcor = part.get_wcor(index)
            local_prop = part.to_atomic_grid(index, local_property)
            condensed[index] = at_grid.integrate(at_weight, local_prop, wcor)
        return condensed
