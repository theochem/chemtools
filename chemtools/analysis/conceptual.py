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
from horton import IOData, BeckeMolGrid, ProAtomDB
from horton.scripts.wpart import wpart_schemes
from chemtools.tool.globaltool import LinearGlobalTool, QuadraticGlobalTool, ExponentialGlobalTool, RationalGlobalTool
from chemtools.tool.localtool import LinearLocalTool, QuadraticLocalTool
from chemtools.tool.condensedtool import CondensedTool
from chemtools.utils import CubeGen


class BaseConceptualDFT(object):
    '''
    Base class for conceptual density functional theory (DFT) analysis of quantum
    chemistry output file(s).
    '''
    def __init__(self, numbers, pseudo_numbers, coordinates, grid, model, energy_0, energy_p, energy_m, density_0, density_p, density_m, n_elec):
        '''
        Parameters
        ----------
        number : np.ndarray
            Atomic number array of atoms in the molecule.
        pseudo_numbers : np.ndarray
            Pseudo-number array of atoms in the molecule.
        coordinates : np.ndarray
            Cartesian coordinate array of atoms in the molecule.
        grid : instance of `BeckeMolGrid` or `CubeGen`
            Moleculr grid.
        model : str
            Interpolation model.
        energy_0 : float
            Energy of :math:`N_0`-electron system.
        energy_p : float
            Energy of :math:`(N_0+1)`-electron system.
        energy_m : float
            Energy of :math:`(N_0-1)`-electron system.
        density_0 : np.ndarray
            Electron density of :math:`N_0`-electron system.
        density_p : np.ndarray
            Electron density of :math:`(N_0+1)`-electron system.
        density_m : np.ndarray
            Electron density of :math:`(N_0-1)`-electron system.
        n_elec : float
            Number of electrons in reference system, i.e. :math:`N_0`.
        '''
        self._numbers = numbers
        self._pseudo_numbers = pseudo_numbers
        self._coordinates = coordinates
        self._grid = grid

        # Assign interpolation model
        if model not in ['linear', 'quadratic', 'exponential', 'rational']:
            raise ValueError('Argument model={0} is not supported.'.format(model))
        # if model is 'general' and energy_expr is None:
        #     raise ValueError('Argument energy_expr is required when model=\'general\'.')
        self._model = model

        dict_global = {'linear':LinearGlobalTool, 'quadratic':QuadraticGlobalTool,
                       'exponential':ExponentialGlobalTool, 'rational':RationalGlobalTool}
        dict_local = {'linear':LinearLocalTool, 'quadratic':QuadraticLocalTool}

        if self._model in dict_global.keys():
            self._globaltool = dict_global[self._model](energy_0, energy_p, energy_m, n_elec)
        else:
            self._globaltool = None

        if self._model in dict_local.keys():
            self._localtool = dict_local[self._model](density_0, density_p, density_m, n_elec)
        else:
            self._localtool = None

    @property
    def model(self):
        '''
        Energy model used to calculate descriptive tools.
        '''
        return self._model

    @property
    def numbers(self):
        '''
        Atomic number array of atoms in the molecule.
        '''
        return self._numbers

    @property
    def pseudo_numbers(self):
        '''
        Pseudo-number array of atoms in the molecule.
        '''
        return self._pseudo_numbers

    @property
    def coordinates(self):
        '''
        Cartesian coordinate array of atoms in the molecule.
        '''
        return self._coordinates

    @property
    def grid(self):
        '''
        Moleculr grid for computing local properties.
        '''
        return self._grid

    @property
    def globaltool(self):
        '''
        Instance of :mod:`chemtools.tool.globaltool`.
        '''
        return self._globaltool

    @property
    def localtool(self):
        '''
        Instance of :mod:`chemtools.tool.localtool`.
        '''
        return self._localtool

    @property
    def condensedtool(self):
        '''
        Instance of :mod:`chemtools.tool.localtool`.
        '''
        return self._condensedtool


class ConceptualDFT_1File(BaseConceptualDFT):
    '''
    Class for conceptual density functional theory (DFT) analysis of one quantum
    chemistry output file using the frontiner molecular orbital (FMO) approach.
    '''
    def __init__(self, molecule, model='quadratic', grid=None, part_scheme='b', proatoms=None):
        '''
        Parameters
        ----------
        molecule : str
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
        grid : instance of horton.BeckeMolGrid, default=None
            The BeckeMolGrid grid to calculate local tools. (under what conditions it can be a qubic grid as well.)
        part_scheme : str, default='b'
            Partitioning scheme used for condensing local descriptors.
        proatoms :
            Reference pro-atoms.
        '''
        # Load molecule
        if not isinstance(molecule, IOData):
            molecule = IOData.from_file(molecule)
        self._mol = molecule
        numbers = self._mol.numbers
        pseudo_numbers = self._mol.pseudo_numbers
        coordinates = self._mol.coordinates

        # Check or setup a default grid
        if grid is None:
            grid = BeckeMolGrid(coordinates, numbers, pseudo_numbers, agspec='fine',
                                random_rotate=False, mode='keep')
        else:
            assert np.all(abs(grid.numbers - numbers) < 1.e-6)
            assert np.all(abs(grid.pseudo_numbers - pseudo_numbers) < 1.e-6)
            assert np.all(abs(grid.coordinates - coordinates) < 1.e-6)

        #
        # Frontiner Molecular Orbital (FMO) Approach
        #

        # Get HOMO & LUMO indices, energies & densities
        homo_index = self._mol.exp_alpha.get_homo_index()
        lumo_index = self._mol.exp_alpha.get_lumo_index()
        homo_energy = self._mol.exp_alpha.homo_energy
        lumo_energy = self._mol.exp_alpha.lumo_energy
        homo_dens = self._mol.obasis.compute_grid_orbitals_exp(self._mol.exp_alpha,
                                                               grid.points,
                                                               np.array([homo_index]))**2
        lumo_dens = self._mol.obasis.compute_grid_orbitals_exp(self._mol.exp_alpha,
                                                               grid.points,
                                                               np.array([lumo_index]))**2
        n_elec = int(np.sum(self._mol.exp_alpha.occupations))
        # HACK: self._mol might not have the exp_beta attribute, making it crash
        beta = False
        if hasattr(self._mol, 'exp_beta'):
            if self._mol.exp_beta is not None:
                beta = True
        if beta:
            n_elec += int(np.sum(self._mol.exp_beta.occupations))
            homo_index_beta = self._mol.exp_beta.get_homo_index()
            lumo_index_beta = self._mol.exp_beta.get_lumo_index()
            if self._mol.exp_beta.homo_energy > homo_energy:
                homo_energy = self._mol.exp_beta.homo_energy
                homo_index = homo_index_beta
                homo_dens = self._mol.obasis.compute_grid_orbitals_exp(self._mol.exp_beta,
                                                                       grid.points,
                                                                       np.array([homo_index]))**2
            if self._mol.exp_beta.lumo_energy < lumo_energy:
                lumo_energy = self._mol.exp_beta.lumo_energy
                lumo_index = lumo_index_beta
                lumo_dens = self._mol.obasis.compute_grid_orbitals_exp(self._mol.exp_beta,
                                                                       grid.points,
                                                                       np.array([lumo_index]))**2
        homo_dens = homo_dens.flatten()
        lumo_dens = lumo_dens.flatten()

        # Compute E(N), E(N+1), & E(N-1)
        # Temporary check as HORTON does not store energy when reading WFN files.
        if hasattr(self._mol, 'energy'):
            energy = self._mol.energy
        else:
            raise ValueError('Argument molecule_filename does not contain energy value!')
        energy_plus = energy - lumo_energy
        energy_minus = energy - homo_energy

        # Compute rho(N), i.e. molecule's density
        dm_full = self._mol.get_dm_full()
        density = self._mol.obasis.compute_grid_density_dm(dm_full, grid.points)
        # Compute rho(N+1) & rho(N-1)
        density_plus = density + lumo_dens
        density_minus = density - homo_dens

        # Paritition electron density of N-electron system
        if isinstance(grid, BeckeMolGrid):
            if part_scheme not in wpart_schemes:
                raise ValueError('Argument part_scheme={0} should be one of wpart_schemes={1}.'.format(part_scheme, wpart_schemes))
            wpart = wpart_schemes[part_scheme]
            # TODO: add wpart.options to the kwargs
            kwargs = {}
            if isinstance(proatoms, ProAtomDB):
                kwargs['proatomdb'] = proatoms
            elif proatoms is not None:
                proatomdb = ProAtomDB.from_file(proatoms)
                proatomdb.normalize()
                kwargs['proatomdb'] = proatoms
            dens_part = wpart(coordinates, numbers, pseudo_numbers, grid, density, **kwargs)
            dens_part.do_all()

        elif isinstance(grid, CubeGen):
            raise NotImplementedError('Should implement partitioning with cubic grid!')

        else:
            raise ValueError('Argument grid={0} is not supported by partitioning class.'.format(grid))

        self._condensedtool = CondensedTool(dens_part)

        BaseConceptualDFT.__init__(self, numbers, pseudo_numbers, coordinates, grid, model,
                                   energy, energy_plus, energy_minus, density, density_plus, density_minus, n_elec)

    @property
    def condensedtool(self):
        '''
        Instance of :mod:`chemtools.tool.localtool`.
        '''
        return self._condensedtool
