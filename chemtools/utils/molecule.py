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
"""The Input-Output (IO) Module."""

import numpy as np
from horton import IOData

__all__ = ['Molecule']


class WaveFunction(object):
    """
    """
    def __init__(self, iodata):
        """
        """
        self._iodata = iodata

    def __getattr__(self, attr):
        """
        """
        value = getattr(self._iodata, attr, None)
        # if value is None:
        #     raise AttributeError('Attribute {0} does not exist!'.format(attr))
        return value

    @property
    def exp_alpha(self):
        """Alpha orbital expansion."""
        return self.obasis.exp_alpha

    @property
    def exp_beta(self):
        """Beta orbital expansion."""
        if hasattr(self._iodata, 'exp_beta') and self._iodata.exp_beta is not None:
            exp_beta = self._iodata.exp_beta
        else:
            exp_beta = self._iodata.exp_alpha
        return exp_beta

    @property
    def obasis(self):
        """Orbital basis instance."""
        return self._iodata.obasis

    @property
    def nbasis(self):
        """Number of basis functions."""
        return self.obasis.nbasis

    @property
    def nelectrons(self):
        """Number of electrons."""
        nelec = int(np.sum(self.exp_alpha.occupations))
        nelec += int(np.sum(self.exp_beta.occupations))
        return nelec

    @property
    def homo_index(self):
        """Index of alpha and beta HOMO orbital."""
        return self.exp_alpha.get_homo_index(), self.exp_beta.get_homo_index()

    @property
    def lumo_index(self):
        """Index of alpha and beta LUMO orbital."""
        return self.exp_alpha.get_lumo_index(), self.exp_beta.get_homo_index()

    @property
    def homo_energy(self):
        """Energy of alpha and beta HOMO orbital."""
        return self.exp_alpha.homo_energy, self.exp_beta.homo_energy

    @property
    def lumo_energy(self):
        """Energy of alpha and beta LUMO orbital."""
        return self.exp_alpha.lumo_energy, self.exp_beta.lumo_energy

    @property
    def orbital_occupation(self):
        """Orbital occupation of alpha and beta electrons."""
        return self.exp_alpha.occupations, self.exp_beta.occupations

    @property
    def orbital_energy(self):
        """Orbital energy of alpha and beta electrons."""
        return self.exp_alpha.energies, self.exp_beta.energies

    @property
    def orbital_coefficient(self):
        """Orbital coefficient of alpha and beta electrons."""
        return self.exp_alpha.coeffs, self.exp_beta.coeffs

    @property
    def dm_fulll(self):
        """ """
        return self._iodata.get_dm_full()

    @property
    def compute_density(self, points, spin=None):
        """Return electronic density evaluated on points."""
        if spin is None:
            dens = self.obasis.compute_grid_density_dm(self.dm_full, points)
        else:
            raise NotImplementedError('')
        return dens

    def compute_homo_density(self, points, spin='alpha', degeneracy=True):
        """
        """
        # index and expression of lumo orbital
        if spin == 'alpha':
            index = self.homo_index[0]
            exp = self.exp_alpha
        elif spin == 'beta':
            index = self.homo_index[1]
            exp = self.exp_beta
        else:
            raise ValueError('Spin={0} is not recognized.'.format(spin))
        # compute density
        dens = self.obasis.compute_grid_orbitals_exp(exp, points, np.array([index]))**2
        return dens.flatten()

    def compute_lumo_density(self, points, spin='alpha', degeneracy=True):
        """
        """
        spin_index = {'alpha': 0, 'beta': 1}
        if spin not in spin_index.keys():
            raise ValueError('Spin={0} is not recognized.'.format(spin))
        # index and expression of lumo orbital
        index = self.lumo_index[spin_index[spin]]
        exp = getattr(self.obasis, 'exp_' + spin)
        # compute density
        dens = self.obasis.compute_grid_orbitals_exp(exp, points, np.array([index]))**2
        return dens.flatten()

    def compute_orbital_density(self, points, index, spin='alpha', degeneracy=True):
        """
        """
        pass


class Molecule(object):
    """
    """
    def __init__(self, coordinates, numbers, iodata=None):
        """
        """
        # check whether iodata has wave-function infromation
        if iodata is not None and hasattr(iodata, 'obasis') and iodata.obasis is not None:
            self._wavefunction = WaveFunction(iodata)
        else:
            self._wavefunction = None
        # bare minimum to define a molecule
        self._coordinates = iodata.coordinates
        self._numbers = iodata.numbers

    def __getattr__(self, attr):
        """ """
        if self._wavefunction is None:
            return None
        value = getattr(self._wavefunction, attr, None)
        return value

    @classmethod
    def from_file(cls, filename):
        """
        Initialize class give a file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        """
        iodata = IOData.from_file(filename)
        return cls(iodata.coordinates, iodata.numbers, iodata)

    @property
    def coordinates(self):
        """ """
        return self._coordinates

    @property
    def numbers(self):
        """ """
        return self._numbers
