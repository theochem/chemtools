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
"""Wrapper Module."""

import numpy as np
from horton import IOData, DenseLinalgFactory
from chemtools.utils.molecule import BaseMolecule

__all__ = ['HortonMolecule']


class HortonMolecule(BaseMolecule):
    """Molecule class from HORTON package."""

    def __init__(self, iodata):
        """
        Initialize class.

        Parameters
        ----------
        iodata : horton.IOData
           An instance of horton.IOData object.
        """
        self._iodata = iodata
        # initialize base class
        super(HortonMolecule, self).__init__(self._iodata.coordinates, self._iodata.numbers)
        # assign alpha orbital expression
        self._exp_alpha = self._iodata.exp_alpha
        # assign beta orbital expression
        if hasattr(self._iodata, 'exp_beta') and self._iodata.exp_beta is not None:
            self._exp_beta = self._iodata.exp_beta
        else:
            self._exp_beta = self._iodata.exp_alpha

    @classmethod
    def from_file(cls, filename):
        """
        Initialize class given a file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        """
        # load molecule
        iodata = IOData.from_file(filename)
        return cls(iodata)

    def __getattr__(self, attr):
        """
        Return attribute.

        Parameters
        ----------
        attr : str
            The name of attribute to retrieve.
        """
        value = getattr(self._iodata, attr, None)
        return value

    @property
    def nbasis(self):
        """Number of basis functions."""
        return self._iodata.obasis.nbasis

    @property
    def nelectrons(self):
        """Number of alpha and beta electrons."""
        return np.sum(self._exp_alpha.occupations), np.sum(self._exp_beta.occupations)

    @property
    def homo_index(self):
        """Index of alpha and beta HOMO orbital."""
        # HORTON indexes the orbitals from 0, so 1 is added to get the intuitive index
        return self._exp_alpha.get_homo_index() + 1, self._exp_beta.get_homo_index() + 1

    @property
    def lumo_index(self):
        """Index of alpha and beta LUMO orbital."""
        # HORTON indexes the orbitals from 0, so 1 is added to get the intuitive index
        return self._exp_alpha.get_lumo_index() + 1, self._exp_beta.get_lumo_index() + 1

    @property
    def homo_energy(self):
        """Energy of alpha and beta HOMO orbital."""
        return self._exp_alpha.homo_energy, self._exp_beta.homo_energy

    @property
    def lumo_energy(self):
        """Energy of alpha and beta LUMO orbital."""
        return self._exp_alpha.lumo_energy, self._exp_beta.lumo_energy

    @property
    def orbital_occupation(self):
        """Orbital occupation of alpha and beta electrons."""
        return self._exp_alpha.occupations, self._exp_beta.occupations

    @property
    def orbital_energy(self):
        """Orbital energy of alpha and beta electrons."""
        return self._exp_alpha.energies, self._exp_beta.energies

    @property
    def orbital_coefficient(self):
        """Orbital coefficient of alpha and beta electrons.

        The alpha and beta orbital coefficients are each storied in a 2d-array in which
        the columns represent the basis coefficients of each molecular orbital.
        """
        return self._exp_alpha.coeffs, self._exp_beta.coeffs

    def compute_overlap(self):
        """Return the overlap matrix of molecular orbitals."""
        # make linear algebra factory
        lf = DenseLinalgFactory(self.nbasis)
        # compute overlap matrix
        arr = self._iodata.obasis.compute_overlap(lf)._array
        return arr

    def compute_density_matrix_array(self, spin='ab'):
        """
        Return the density matrix array for the specified spin orbitals.
        """
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        return dm._array

    def _get_density_matrix(self, spin):
        """
        Return HORTON density matrix object corresponding to the specified spin.

        Parameters
        ----------
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons
        """
        if spin not in ['a', 'b', 'alpha', 'beta', 'ab']:
            raise ValueError('Argument spin is not recognized!')

        if spin == 'ab':
            # get density matrix of alpha & beta electrons
            dm = self._iodata.get_dm_full()
        else:
            # get orbital expression of specified spin
            spin_type = {'a': 'alpha', 'alpha': 'alpha', 'b': 'beta', 'beta': 'beta'}
            exp = getattr(self._iodata.obasis, 'exp_' + spin_type[spin])
            # get density matrix of specified spin
            dm = exp.to_dm()
        return dm

    def compute_density(self, points, spin='ab', orbital_index=None, output=None):
        r"""
        Return electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons

        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), float)
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)

        # compute density
        if orbital_index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_density_dm(dm, points, output=output)
        else:
            # HORTON index the orbitals from 0
            orbs = np.copy(np.asarray(orbital_index)) - 1
            # include specified set of orbitals
            self._iodata.obasis.compute_grid_orbitals_exp(dm, points, orbs, output=output)**2
            output = output.flatten()
        return output

    def compute_gradient(self, points, spin='ab', orbital_index=None, output=None):
        r"""
        Return gradient of electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons

        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n, 3) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0], 3), float)
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)

        # compute gradient
        if orbital_index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_gradient_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_hessian(self, points, spin='ab', orbital_index=None, output=None):
        r"""
        Return hessian of electron density evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons

        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n, 6) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)

        # compute hessian
        if orbital_index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_hessian_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_esp(self, points, spin='ab', orbital_index=None, output=None, charges=None):
        r"""
        Return the molecular electrostatic potential on the given points for the specified spin.

        The molecular electrostatic potential at point :math:`\mathbf{r}` is caused by the
        electron density :math:`\rho` of the specified spin orbitals and set of point charges
        :math:`\{q_A\}_{A=1}^{N_\text{atoms}}` placed at the position of the nuclei. i.e,

        .. math::
           V \left(\mathbf{r}\right) =
             \sum_{A=1}^{N_\text{atoms}} \frac{q_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
             \int \frac{\rho \left(\mathbf{r}'\right)}{\rvert \mathbf{r}' - \mathbf{r} \lvert}
                  d\mathbf{r}'

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str, default='ab'
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons

        orbital_index : sequence, default=None
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray, default=None
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        charges : np.ndarray, default=None
           Array with shape (n,) representing the point charges at the position of the nuclei.
           When ``None``, the pseudo numbers are used.
        """
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), float)
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        # assign point charges
        if charges is None:
            charges = self.pseudo_numbers
        elif charges.shape != self.numbers.shape:
            raise ValueError('Argument charges should have {0} shape.'.format(self.numbers.shape))

        # compute esp
        if orbital_index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_esp_dm(dm, self.coordinates, charges, points,
                                                    output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_kinetic_energy_density(self, points, spin='ab', orbital_index=None, output=None):
        r"""
        Return positive definite kinetic energy density on the given points for the specified spin.

        Positive definite kinetic energy density is defined as,

        .. math::
           \tau \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - 'a' or 'alpha': consider alpha electrons
           - 'b' or 'beta': consider beta electrons
           - 'ab': consider alpha and beta electrons

        orbital_index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbtails are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), float)
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)

        # compute kinetic energy
        if orbital_index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_kinetic_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output
