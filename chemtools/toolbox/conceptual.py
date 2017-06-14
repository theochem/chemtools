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
"""
Module for Conceptual Density Functional Theory Analysis of Quantum Chemistry Output Files.

This modules contains wrappers which take outputs of quantum chemistry software and
compute various conceptual density functional theory (DFT) descriptive tools.
"""

import numpy as np
from horton import log
from horton import BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes
from chemtools.utils.molecule import BaseMolecule
from chemtools.toolbox.molecule import make_molecule
from chemtools.conceptual.linear import LinearGlobalTool, LinearLocalTool
from chemtools.conceptual.quadratic import QuadraticGlobalTool, QuadraticLocalTool
from chemtools.conceptual.exponential import ExponentialGlobalTool
from chemtools.conceptual.rational import RationalGlobalTool
from chemtools.conceptual.general import GeneralGlobalTool

__all__ = ['GlobalConceptualDFT', 'LocalConceptualDFT', 'CondensedConceptualDFT']


class BaseConceptualDFT(object):
    """Base class for conceptual density functional theory (DFT) analysis."""

    def __init__(self, dict_values, dict_models, model, coordinates, numbers):
        """
        Initialize class.

        Parameters
        ----------
        dict_values : dict
            Dictionary of number_electron:property_value.
        dict_models : dict
            Dictionary of energy_model:energy_object.
        model : str
            Energy model.
        coordinates : np.ndarray
            Coordinates of atomic centers.
        numbers : np.ndarray
            Atomic number of atomic centers.
        """
        if model.lower() not in dict_models.keys():
            raise ValueError('Model={0} is not available!'.format(model.lower()))
        self._model = model.lower()

        # check shape of coordinates
        if coordinates is not None and (coordinates.ndim != 2 or coordinates.shape[1] != 3):
            raise ValueError('Argument coordinate should be a 2D-array with 3 columns! '
                             'shape={0}'.format(coordinates.shape))

        # check number of atoms given by numbers and coordinates match
        if numbers is not None and coordinates is not None and len(numbers) != len(coordinates):
            raise ValueError('Numbers & coordinates should represent same number of atoms! '
                             '{0}!={1}'.format(len(numbers), len(coordinates)))

        self._coordinates = coordinates
        self._numbers = numbers

        if self._model != 'general':
            # get sorted number of electrons
            nelectrons = sorted(dict_values.keys())
            # get value for sorted number of electrons
            values = [dict_values[key] for key in nelectrons]

            if len(nelectrons) != 3:
                raise NotImplementedError(
                    'For model={0}, three values are required!'.format(self._model))

            # check number of electrons are integers
            if not all(isinstance(item, int) or
                       (isinstance(item, float) and (item).is_integer()) for item in nelectrons):
                raise ValueError('For model={0}, integer number of electrons are required! '
                                 '#electrons={1}'.format(self._model, nelectrons))

            # check consecutive number of electrons change by one
            if not all([y - x == 1 for x, y in zip(nelectrons[:-1], nelectrons[1:])]):
                raise ValueError('For model={0}, consecutive number of electrons should differ by '
                                 '1! #electrons={1}'.format(self._model, nelectrons))

            # obtain reference number of electrons
            if len(nelectrons) % 2 == 0:
                raise NotImplementedError('Even number of molecules is not implemented yet!')
            else:
                # consider the middle system as reference
                reference = nelectrons[len(nelectrons) // 2]

            # make a list of arguments for global tool
            args = [values[1], values[2], values[0], reference]
            self._tool = dict_models[self._model](*args)
        else:
            # self._tool = select_tool[model](*args, **kwargs)
            raise NotImplementedError('Model={0} is not covered yet!'.format(self._model))

        # print screen information
        self._log_init()

    def __getattr__(self, attr):
        """Return class attribute."""
        value = getattr(self._tool, attr, 'error')
        if isinstance(value, str) and value == 'error':
            raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return getattr(self._tool, attr)

    @property
    def model(self):
        """Energy model."""
        return self._model

    @property
    def coordinates(self):
        """Cartesian coordinates of atoms."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic numbers."""
        return self._numbers

    @staticmethod
    def load_file(filenames):
        """
        Return `BaseMolecule` instances corresponding to filenames.

        Parameters
        ----------
        filenames : str, list, tuple
            String specifying the path to molecule's file, or list/tuple of strings specifying
            path to molecule files.
        """
        if isinstance(filenames, (str, unicode)):
            # case of one file not given as a list
            molecule = make_molecule(filenames)
        else:
            # case of multiple file(s)
            molecule = [make_molecule(filename) for filename in filenames]
        return molecule

    @staticmethod
    def load_molecule(molecules, points=None, part=None):
        """
        Initialize class from `BaseMolecule` objects.

        Parameters
        ----------
        molecules : instance, list, tuple
            Instance of `BaseMolecule` object or list/tuple of `BaseMolecule` objects
        model : str
            Energy model.
        points : np.array
            Points on which the local properties are evaluated.
        part : dict
            Dictionary specifying condensing details. This includes:
                'scheme' : str
                    Name of partitioning scheme
                'approach' : str
                    Condensing approach, select 'fmo' or 'rmf'.
                'grid' : instance
                    Instance of BeckeMolGrid molecular grid.
                'proatomdb' : instance
                    Instance of ProatomDB proatom database.
                'kwargs' : list
                    Extra keyword arguments.
        """
        # check points array
        if points is not None and (points.ndim != 2 or points.shape[1] != 3):
            raise ValueError('Argument points should be a 2D array with 3 columns! '
                             'given shape={0}'.format(points.shape))
        if isinstance(molecules, BaseMolecule):
            # case of one BaseMolecule object
            return BaseConceptualDFT.load_frontier_molecular_orbital(molecules, points, part)
        elif len(molecules) == 1:
            # case of one BaseMolecule object given in a list
            return BaseConceptualDFT.load_frontier_molecular_orbital(molecules[0], points, part)
        else:
            # case of multiple BaseMolecule objets
            return BaseConceptualDFT.load_finite_difference(molecules, points, part)

    @staticmethod
    def load_frontier_molecular_orbital(molecule, points=None, part=None):
        """
        Load molecule and return its properties within frontier molecular orbital theory approach.

        Parameters
        ----------
        molecules : list,tuple
            List/tuple containing `BaseMolecule` objects.
        model : str
            Energy model.
        points : np.array
            Points on which the local properties are evaluated.
        part : dict
            Dictionary specifying condensing details. This includes:
                'scheme' : str
                    Name of partitioning scheme
                'approach' : str
                    Condensing approach, select 'fmo' or 'rmf'.
                'grid' : instance
                    Instance of BeckeMolGrid molecular grid.
                'proatomdb' : instance
                    Instance of ProatomDB proatom database.
                'kwargs' : list
                    Extra keyword arguments.
        """
        # get homo & lumo energy and spin
        homo_energy = molecule.homo_energy[0]
        homo_spin = 'a'
        if molecule.homo_energy[1] > homo_energy:
            homo_energy = molecule.homo_energy[1]
            homo_spin = 'b'
        lumo_energy = molecule.lumo_energy[0]
        lumo_spin = 'a'
        if molecule.lumo_energy[1] > lumo_energy:
            lumo_energy = molecule.lumo_energy[1]
            lumo_spin = 'b'
        nelec = sum(molecule.nelectrons)

        # store number of electron and energy in a dictionary
        dict_energy = dict([(nelec, molecule.energy),
                            (nelec + 1, molecule.energy + lumo_energy),
                            (nelec - 1, molecule.energy - homo_energy)])

        # compute and record densities on given points in a dictionary
        dict_dens = {}
        if points is not None:
            # compute homo & lumo density
            spin_to_index = {'a': 0, 'b': 1}
            homo_dens = molecule.compute_density(points, homo_spin,
                                                 molecule.homo_index[spin_to_index[homo_spin]])
            lumo_dens = molecule.compute_density(points, lumo_spin,
                                                 molecule.lumo_index[spin_to_index[lumo_spin]])
            # store number of electron and density in a dictionary
            dens = molecule.compute_density(points, 'ab', None)
            dict_dens = dict([(nelec, dens),
                              (nelec + 1, dens + lumo_dens),
                              (nelec - 1, dens - homo_dens)])

        # compute and record populations using given grid in a dictionary
        if part is not None:
            # get molecular grid
            if part['grid'] is not None:
                if not isinstance(part['grid'], BeckeMolGrid):
                    raise ValueError('Only BeckeMolGrid is supported for condensing!')
                else:
                    grid = part['grid']
                    # check coordinates of grid and molecule match
                    if not np.all(abs(grid.centers - molecule.coordinates) < 1.e-4):
                        raise ValueError('Coordinates of grid and molecule do not match!')
                    if not np.all(abs(grid.numbers - molecule.numbers) < 1.e-4):
                        raise ValueError('Atomic numbers of grid and molecule do not match!')
            else:
                # make default grid
                grid = BeckeMolGrid(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                                    agspec='fine', random_rotate=False, mode='keep')
            # compute density of reference system on grid point
            dens = molecule.compute_density(grid.points, 'ab', None)
            # compute homo & lumo density on grid points
            spin_to_index = {'a': 0, 'b': 1}
            homo_dens = molecule.compute_density(grid.points, homo_spin,
                                                 molecule.homo_index[spin_to_index[homo_spin]])
            lumo_dens = molecule.compute_density(grid.points, lumo_spin,
                                                 molecule.lumo_index[spin_to_index[lumo_spin]])
            # get partitioning class
            if part['scheme'].lower() not in wpart_schemes:
                raise ValueError('Partitioning scheme={0} not supported! '
                                 'Select from: {1}'.format(part['scheme'], wpart_schemes.keys()))
            wpart = wpart_schemes[part['scheme']]
            # make proatom database
            kwargs = {}
            if part['scheme'].lower() not in ['mbis', 'b']:
                if part['proatomdb'] is None:
                    proatomdb = ProAtomDB.from_refatoms(molecule.numbers)
                else:
                    proatomdb = part['proatomdb']
                kwargs['proatomdb'] = proatomdb
            # compute population of reference system
            part0 = wpart(molecule.coordinates, molecule.numbers,
                          molecule.pseudo_numbers, grid, dens, **kwargs)
            part0.do_all()
            pops0 = part0['populations']
            # compute population of (N + 1) and (N - 1) electron systems
            if part['approach'].lower() == 'fmr':
                # fragment of molecular response
                popsp = CondensedConceptualDFT.condense_to_atoms(dens + lumo_dens, part0)
                popsm = CondensedConceptualDFT.condense_to_atoms(dens - homo_dens, part0)
            elif part['approach'].lower() == 'rmf':
                # response of molecular fragment
                partp = wpart(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                              grid, dens + lumo_dens, **kwargs)
                partp.do_all()
                popsp = partp['populations']
                partm = wpart(molecule.coordinates, molecule.numbers, molecule.pseudo_numbers,
                              grid, dens - homo_dens, **kwargs)
                partm.do_all()
                popsm = partm['populations']
            # Store number of electron and populations in a dictionary
            dict_pops = dict([(nelec, pops0),
                              (nelec + 1, popsp),
                              (nelec - 1, popsm)])

        # get coordinates of molecule
        coordinates = molecule.coordinates
        # get atomic numbers of molecule
        atomic_numbers = molecule.numbers

        if points is None and part is None:
            return coordinates, atomic_numbers, dict_energy
        elif part is None:
            return coordinates, atomic_numbers, dict_energy, dict_dens
        elif points is None:
            return coordinates, atomic_numbers, dict_energy, dict_pops
        # else:
        #     return coordinates, atomic_numbers, dict_energy, dict_dens, dict_pops

    @staticmethod
    def load_finite_difference(molecules, points=None, part=None):
        """
        Load molecules and return their properties within finite difference theory approach.

        Parameters
        ----------
        molecules : list,tuple
            List/tuple containing `BaseMolecule` objects.
        model : str
            Energy model.
        points : np.array
            Points on which the local properties are evaluated.
        part : dict
            Dictionary specifying condensing details. This includes:
                'scheme' : str
                    Name of partitioning scheme
                'approach' : str
                    Condensing approach, select 'fmo' or 'rmf'.
                'grid' : instance
                    Instance of BeckeMolGrid molecular grid.
                'proatomdb' : instance
                    Instance of ProatomDB proatom database.
                'kwargs' : list
                    Extra keyword arguments.
        """
        dict_energy, dict_moleculea = {}, {}
        for index, molecule in enumerate(molecules):
            # check atomic numbers
            if index == 0:
                atomic_numbers = molecule.numbers
            else:
                reference, numbers = atomic_numbers, molecule.numbers
                if reference.shape != numbers.shape or not np.all(abs(reference - numbers) < 1.e-6):
                    raise ValueError(
                        'Molecule 1 & {0} have different atomic numbers!'.format(index + 1))

            # get number of electrons
            nelec = sum(molecule.nelectrons)
            if nelec in dict_energy.keys():
                raise ValueError('Two molecules have {0} electrons!'.format(nelec))

            # store number of electron and iodata in a dictionary
            dict_moleculea[nelec] = molecule

            # get and store energy
            # if not hasattr(molecule, 'energy'):
            #     raise ValueError('Molecule object does not contain energy value!')
            # store number of electron and energy in a dictionary
            dict_energy[nelec] = molecule.energy

            # compute and record densities on given points in a dictionary
            if points is not None:
                # check coordinates
                if index == 0:
                    dict_dens = {}
                    coordinates = molecule.coordinates
                else:
                    if not np.all(abs(coordinates - molecule.coordinates) < 1.e-4):
                        raise ValueError(
                            'Molecule 1 & {0} have different geometries!'.format(index + 1))
                # computer & store density for a given number of electron
                dens = molecule.compute_density(points, 'ab', None)
                dict_dens[nelec] = dens

            if part is not None:
                # make molecular grid
                if part['grid'] is not None:
                    if not isinstance(part['grid'], BeckeMolGrid):
                        raise ValueError('Only BeckeMolGrid is supported for condensing!')
                    else:
                        grid = part['grid']
                        # check coordinates of grid and molecule match
                        ref, cord = grid.centers, molecule.coordinates
                        if ref.shape != cord.shape or not np.all(abs(ref - cord) < 1.e-4):
                            raise ValueError('Coordinates of grid and molecule {0} '
                                             'should match!'.format(index + 1))
                else:
                    # default grid
                    grid = BeckeMolGrid(molecule.coordinates, molecule.numbers,
                                        molecule.pseudo_numbers, agspec='fine',
                                        random_rotate=False, mode='keep')
                # check geometries
                if index == 0:
                    coordinates = molecule.coordinates
                    same_coordinates = True
                elif not np.all(abs(coordinates - molecule.coordinates) < 1.e-4):
                    same_coordinates = False
                    # if geometries are not the same, only rmf can be done.
                    if part['approach'].lower() == 'fmr':
                        raise ValueError('When geometries of molecules are different, only '
                                         'approach=RMF is possible!')
        # Get sorted number of electrons
        nelectrons = sorted(dict_moleculea.keys())
        if len(nelectrons) != 3:
            raise ValueError('Conceptual DFT within FD approach, currently '
                             'only works for 3 molecules!')
        # number of electrons of reference system
        reference = nelectrons[1]
        # coordinates and atomic numbers
        coordinates = dict_moleculea[reference].coordinates
        atomic_numbers = dict_moleculea[reference].numbers

        if part is not None:
            # get partitioning class
            if part['scheme'].lower() not in wpart_schemes:
                raise ValueError('Partitioning scheme={0} not supported!'
                                 'Select from: {1}'.format(part['scheme'], wpart_schemes))
            wpart = wpart_schemes[part['scheme']]
            kwargs = part['kwargs']
            if part['scheme'].lower() not in ['mbis', 'b']:
                if part['proatomdb'] is None:
                    proatomdb = ProAtomDB.from_refatoms(atomic_numbers)
                kwargs['proatomdb'] = proatomdb

            # compute population of reference molecule
            mol0 = dict_moleculea[reference]
            dens0 = mol0.compute_density(grid.points, 'ab', None)
            part0 = wpart(mol0.coordinates, mol0.numbers,
                          mol0.pseudo_numbers, grid, dens0, **kwargs)
            part0.do_all()
            pops0 = part0['populations']
            # record poputation of reference system
            dict_pops = dict([(reference, pops0)])

            # compute and record populations given grid in a dictionary
            for nelec, molecule in dict_moleculea.iteritems():
                # skip reference system
                if nelec == reference:
                    continue
                # compute population
                if part['approach'].lower() == 'fmr':
                    # fragment of molecular response
                    assert same_coordinates
                    dens = molecule.compute_density(grid.points, 'ab', None)
                    pops = CondensedConceptualDFT.condense_to_atoms(dens, part0)
                elif part['approach'].lower() == 'rmf':
                    # response of molecular fragment
                    if not same_coordinates:
                        grid = BeckeMolGrid(molecule.coordinates, molecule.numbers,
                                            molecule.pseudo_numbers, agspec='fine',
                                            random_rotate=False, mode='keep')
                        points = grid.points
                    dens = molecule.compute_density(grid.points, 'ab', None)
                    parts = wpart(molecule.coordinates, molecule.numbers,
                                  molecule.pseudo_numbers, grid, dens, **kwargs)
                    parts.do_all()
                    pops = parts['populations']
                else:
                    raise ValueError(
                        'Partitioning approach={0} is not recognized!'.format(part['approach']))
                # Store number of electron and populations in a dictionary
                dict_pops[nelec] = pops

        if points is None and part is None:
            return coordinates, atomic_numbers, dict_energy
        elif part is None:
            return coordinates, atomic_numbers, dict_energy, dict_dens
        elif points is None:
            return coordinates, atomic_numbers, dict_energy, dict_pops
        # else:
        #     return coordinates, atomic_numbers, dict_energy, dict_dens, dict_pops


class GlobalConceptualDFT(BaseConceptualDFT):
    """
    Global conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    If only one molecule is provided, the frontier molecular orbital (FMO) approach is invoked,
    otherwise finite difference (FD) approach is taken.
    """

    def __init__(self, dict_values, model='quadratic', coordinates=None, numbers=None):
        """
        Initialize class.

        Parameters
        ----------
        dict_values : dict
            Dictionary of number_electron:energy
        model : str, default='quadratic'
            Energy model used to calculate global descriptive tools.
            The available models include: 'linear', 'quadratic', 'exponential', 'rational',
            and 'general'. Please see '' for more information.
        coordinates : np.ndarray
            Coordinates of atomic centers.
        numbers : np.ndarray
            Atomic number of atomic centers.
        """
        # available models for global tools
        dict_models = {'linear': LinearGlobalTool, 'quadratic': QuadraticGlobalTool,
                       'exponential': ExponentialGlobalTool, 'rational': RationalGlobalTool,
                       'general': GeneralGlobalTool}
        super(GlobalConceptualDFT, self).__init__(dict_values, dict_models, model,
                                                  coordinates, numbers)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith('_')]
        attrs.sort()
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith('_')]
        methods.sort()
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        for attr in attrs:
            value = getattr(self._tool, attr)
            if value is not None and not hasattr(value, '__iter__'):
                content += '\n%s   % .6f' % (attr.ljust(25), value)
            elif value is not None and hasattr(value, '__iter__'):
                content += '\n%s   [%s]' % (attr.ljust(25), ', '.join(['% .6f' % v for v in value]))
            else:
                content += '\n%s   %s' % (attr.ljust(25), ' ---')
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         # ('Parameters', self._tool.params),
                         ('Reference Energy', self._tool.energy_zero),
                         ('Reference #Electrons', self._tool.n0),
                         ('Energy Values', [getattr(self._tool, attr) for attr in
                                            ['_energy_minus', 'energy_zero', '_energy_plus']])])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model):
        """
        Initialize class from files.

        Parameters
        ----------
        filenames : str, list, tuple
            String specifying the path to molecule's file, or list/tuple of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate global properties.
        """
        molecules = cls.load_file(filenames)
        return cls.from_molecule(molecules, model)

    @classmethod
    def from_molecule(cls, molecules, model):
        """
        Initialize class from `BaseMolecule` objects.

        Parameters
        ----------
        molecules : list, tuple
            List/tuple containing `BaseMolecule` objects.
        model : str
            Energy model used to calculate global properties.
        """
        coordinates, numbers, dict_energy = cls.load_molecule(molecules)
        return cls(dict_energy, model, coordinates, numbers)


class LocalConceptualDFT(BaseConceptualDFT):
    """
    Local conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    If only one molecule is provided, the frontier molecular orbital (FMO) approach is invoked,
    otherwise finite difference (FD) approach is taken.

    Note: If FD approach is taken, the geometries of all molecules need to be the same.
    """

    def __init__(self, dict_values, model='quadratic', coordinates=None, numbers=None):
        """
        Initialize class.

        Parameters
        ----------
        dict_value : dict
            Dictionary of number_electron:density
        model : str, default='quadratic'
            Energy model used to calculate local descriptive tools.
            The available models include: 'linear', 'quadratic', and 'general'.
            Please see '' for more information.
        coordinates : np.ndarray
            Coordinates of atomic centers.
        numbers : np.ndarray
            Atomic number of atomic centers.
        """
        # available models for local tools
        dict_models = {'linear': LinearLocalTool, 'quadratic': QuadraticLocalTool}
        # check densities shape
        for index, value in enumerate(dict_values.values()):
            if index == 0:
                shape = value.shape
            elif value.shape != shape:
                raise ValueError(
                    'Densities should have the same shape! {0}!={1}'.format(shape, value.shape))
        super(LocalConceptualDFT, self).__init__(dict_values, dict_models, model,
                                                 coordinates, numbers)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool) + dir(self)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith('_')]
        attrs.sort()
        # remove n0, because it is both an attr of self._tool and self (duplicate)
        # attrs.remove('n0')
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith('_')]
        methods.sort()
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(attrs)
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference #Electrons', self._tool.n0)])
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
        molecules = cls.load_file(filenames)
        return cls.from_molecule(molecules, model, points)

    @classmethod
    def from_molecule(cls, molecules, model, points):
        """
        Initialize class from `BaseMolecule` objects.

        Parameters
        ----------
        filenames : str, list, tuple
            String specifying the path to molecule's files, or list/tuple of strings specifying
            path to molecule files.
        model : str
            Energy model used to calculate local properties.
            Available models are 'linear' and 'quadratic'.
        points : np.array
            Points on which the local properties are evaluated.
        """
        coords, numbers, _, dict_dens = cls.load_molecule(molecules, points=points)
        return cls(dict_dens, model, coords, numbers)


class CondensedConceptualDFT(BaseConceptualDFT):
    """
    Condensed conceptual density functional theory (DFT) analysis of quantum chemistry output files.

    If only one molecule is provided, the frontier molecular orbital (FMO) approach is invoked,
    otherwise finite difference (FD) approach is taken.

    Note: If FD approach is taken, the geometries of molecules can be different only for RMF
          approach.
    """

    def __init__(self, dict_values, model='quadratic', coordinates=None, numbers=None):
        """
        Initialize class.

        Parameters
        ----------
        dict_value : dict
            Dictionary of number_electron:atomic_populations
        model : str, default='quadratic'
            Energy model used to calculate local descriptive tools.
            The available models include: 'linear', 'quadratic', and 'general'.
            Please see '' for more information.
        coordinates : np.ndarray
            Coordinates of atomic centers.
        numbers : np.ndarray
            Atomic number of atomic centers.
        """
        # available models for local tools
        dict_models = {'linear': LinearLocalTool, 'quadratic': QuadraticLocalTool}
        super(CondensedConceptualDFT, self).__init__(dict_values, dict_models, model,
                                                     coordinates, numbers)

    def __repr__(self):
        """Print table of available class attributes and methods."""
        # get list of available attributes and methods
        avs = dir(self._tool)
        # get sorted list of public attributes
        attrs = [atr for atr in avs if not callable(getattr(self, atr)) and not atr.startswith('_')]
        attrs.sort()
        attrs.remove('n0')
        # get sorted list of public methods
        methods = [atr for atr in avs if callable(getattr(self, atr)) and not atr.startswith('_')]
        content = 'Available attributes in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(attrs)
        content += '\n\nAvailable methods in {0} global model:\n{1}\n'.format(self._model, '-' * 50)
        content += '\n'.join(methods) + '\n'
        return content

    def _log_init(self):
        """Print an initial informative message."""
        if log.do_medium:
            log('Initialize: Condensed Class')
            # log('Initialize: %s' % self.__class__)
            log.deflist([('Energy Model', self._model),
                         ('Reference #Electrons', self._tool.n0)])
            log.blank()

    @classmethod
    def from_file(cls, filenames, model, approach='FMR', scheme='h', grid=None, proatomdb=None):
        """
        Initialize class from files.

        Parameters
        ----------
        filenames : str, list, tuple
            String specifying the path to molecule's files, or list/tuple of strings specifying
            path to molecule files.
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
        molecules = cls.load_file(filenames)
        return cls.from_molecule(molecules, model, grid, approach, scheme, proatomdb)

    @classmethod
    def from_molecule(cls, molecules, model, grid=None, approach='FMR', scheme='mbis',
                      proatomdb=None, **kwargs):
        """
        Initialize class from `BaseMolecule` objects.

        Parameters
        ----------
        molecules : list,tuple
            List/tuple containing `BaseMolecule` instances.
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
        part = {'scheme': scheme, 'approach': approach, 'grid': grid, 'proatomdb': proatomdb,
                'kwargs': kwargs}
        coords, numbers, _, dict_pops = cls.load_molecule(molecules, part=part)
        return cls(dict_pops, model, coords, numbers)

    @staticmethod
    def condense_to_atoms(local_property, part):
        r"""
        Return condensed values of the local descriptor.

        Having the local descriptor
        :math:`p_{\text{local}}\left(\mathbf{r}\right)` into atomic contribution :math:`P_A`,
        the condensed value is defined as,

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
