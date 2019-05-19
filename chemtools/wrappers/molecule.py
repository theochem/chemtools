# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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


import logging
import numpy as np
from horton import IOData, DenseLinalgFactory
try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path


__all__ = ["Molecule"]


class Molecule(object):
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

        # if not (isinstance(coordinates, np.ndarray) and coordinates.ndim == 2):
        #     raise TypeError("Argument coordinates should be a 2d-array.")
        # if not (isinstance(numbers, np.ndarray) and numbers.ndim == 1):
        #     raise TypeError("Argument numbers should be a 1d-array.")
        # if coordinates.shape[0] != numbers.size:
        #     raise TypeError("Arguments coordinates and numbers should represent the same number "
        #                     "of atoms! {0} != {1}".format(coordinates.shape[0], numbers.size))
        # if coordinates.shape[1] != 3:
        #     raise TypeError("Argument coordinates should be a 2d-array with 3 columns.")

        self._coordinates = self._iodata.coordinates
        self._numbers = self._iodata.numbers
        if hasattr(self._iodata, 'exp_alpha'):
            # assign alpha orbital expression
            self._exp_alpha = self._iodata.exp_alpha
            # assign beta orbital expression
            if hasattr(self._iodata, 'exp_beta') and self._iodata.exp_beta is not None:
                self._exp_beta = self._iodata.exp_beta
            else:
                self._exp_beta = self._iodata.exp_alpha
        else:
            self._exp_alpha = None
            self._exp_beta = None

    @classmethod
    def from_file(cls, fname):
        """
        Initialize class given a file.

        Parameters
        ----------
        fname : str
            Path to molecule's files.
        """
        # load molecule
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        try:
            iodata = IOData.from_file(str(fname))
        except IOError as _:
            try:
                with path('chemtools.data.examples', str(fname)) as fname:
                    logging.info('Loading {0}'.format(str(fname)))
                    iodata = IOData.from_file(str(fname))
            except IOError as error:
                logging.info(error)
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
    def coordinates(self):
        """Cartesian coordinates of atomic centers."""
        return self._coordinates

    @property
    def numbers(self):
        """Atomic number of atomic centers."""
        return self._numbers

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

    def compute_orbital_overlap(self):
        """Return the overlap matrix of molecular orbitals."""
        # make linear algebra factory
        lf = DenseLinalgFactory(self.nbasis)
        # compute overlap matrix
        arr = self._iodata.obasis.compute_overlap(lf)._array
        return arr

    def compute_density_matrix(self, spin="ab"):
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

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons
        """
        # check orbital spin
        if spin not in ["a", "b", "alpha", "beta", "ab"]:
            raise ValueError("Argument spin={0} is not recognized!".format(spin))
        # compute density matrix
        if spin == "ab":
            # get density matrix of alpha & beta electrons
            dm = self._iodata.get_dm_full()
        else:
            # get orbital expression of specified spin
            spin_type = {"a": "alpha", "alpha": "alpha", "b": "beta", "beta": "beta"}
            exp = getattr(self, "_exp_" + spin_type[spin])
            # get density matrix of specified spin
            dm = exp.to_dm()
        return dm

    def compute_molecular_orbital(self, points, spin, index=None, output=None):
        """
        Return molecular orbitals evaluated on the given points for the spin orbitals.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals.

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons

        index : sequence, default=None
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray, default=None
           Array with shape (n, m) to store the output, where n in the number of points and m
           is the number of molecular orbitals. When ``None`` the array is allocated.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # assign orbital index (HORTON index the orbitals from 0)
        if index is None:
            # include all occupied orbitals of specified spin
            spin_index = {"a": 0, "alpha": 0, "b": 1, "beta": 1}
            index = np.arange(self.homo_index[spin_index[spin]])
        else:
            # include specified set of orbitals
            index = np.copy(np.asarray(index)) - 1
            if index.ndim == 0:
                index = np.array([index])

        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0], index.shape[0]), float)
        npoints, norbs = points.shape[0], index.shape[0]
        if output.shape != (npoints, norbs):
            raise ValueError("Argument output should be a {0} array.".format((npoints, norbs)))

        # get orbital expression of specified spin
        spin_type = {"a": "alpha", "alpha": "alpha", "b": "beta", "beta": "beta"}
        exp = getattr(self, "_exp_" + spin_type[spin])
        # compute mo expression
        self._iodata.obasis.compute_grid_orbitals_exp(exp, points, index, output=output)
        return output

    def compute_density(self, points, spin="ab", index=None, output=None):
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

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), float)
        if output.shape != (points.shape[0],):
            raise ValueError("Argument output should be a {0} array.".format((points.shape[0],)))

        # compute density
        if index is None:
            # get density matrix corresponding to the specified spin
            dm = self._get_density_matrix(spin)
            # include all orbitals
            self._iodata.obasis.compute_grid_density_dm(dm, points, output=output)
        else:
            # include subset of molecular orbitals
            if spin == "ab":
                # compute mo expression of alpha & beta orbitals
                mo_a = self.compute_molecular_orbital(points, "a", index)
                mo_b = self.compute_molecular_orbital(points, "b", index)
                # add density of alpha & beta molecular orbitals
                np.sum(mo_a**2, axis=1, out=output)
                output += np.sum(mo_b**2, axis=1)
            else:
                # compute mo expression of specified molecular orbitals
                mo = self.compute_molecular_orbital(points, spin, index)
                # add density of specified molecular orbitals
                np.sum(mo**2, axis=1, out=output)
        return output

    def compute_gradient(self, points, spin="ab", index=None, output=None):
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

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray
           Array with shape (n, 3) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0], 3), float)
        if output.shape != (points.shape[0], 3):
            raise ValueError("Argument output should be a {0} array.".format((points.shape[0], 3)))

        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        # compute gradient
        if index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_gradient_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_hessian(self, points, spin="ab", index=None, output=None):
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

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray
           Array with shape (n, 6) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0], 6), float)
        if output.shape != (points.shape[0], 6):
            raise ValueError("Argument output should be a {0} array.".format((points.shape[0], 6)))

        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        # compute hessian
        if index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_hessian_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_esp(self, points, spin="ab", index=None, output=None, charges=None):
        r"""
        Return the molecular electrostatic potential on the given points for the specified spin.

        The molecular electrostatic potential at point :math:`\mathbf{r}` is caused by the
        electron density :math:`\rho` of the specified spin orbitals and set of point charges
        :math:`\{q_A\}_{A=1}^{N_\text{atoms}}` placed at the position of the nuclei. i.e,

        .. math::
           V \left(\mathbf{r}\right) =
             \sum_{A=1}^{N_\text{atoms}} \frac{q_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
             \int \frac{\rho \left(\mathbf{r}"\right)}{\rvert \mathbf{r}" - \mathbf{r} \lvert}
                  d\mathbf{r}"

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str, default="ab"
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence, default=None
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray, default=None
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        charges : np.ndarray, default=None
           Array with shape (n,) representing the point charges at the position of the nuclei.
           When ``None``, the pseudo numbers are used.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), np.float)
        if output.shape != (points.shape[0],):
            raise ValueError("Argument output should be a {0} array.".format((points.shape[0],)))

        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        # assign point charges
        if charges is None:
            charges = self.pseudo_numbers
        elif not isinstance(charges, np.ndarray) or charges.shape != self.numbers.shape:
            raise ValueError("Argument charges should be a 1d-array "
                             "with {0} shape.".format(self.numbers.shape))
        # compute esp
        if index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_esp_dm(dm, self.coordinates, charges, points,
                                                    output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_ked(self, points, spin="ab", index=None, output=None):
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

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        output : np.ndarray
           Array with shape (n,) to store the output, where n in the number of points.
           When ``None`` the array is allocated.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")
        # allocate output array
        if output is None:
            output = np.zeros((points.shape[0],), float)
        if output.shape != (points.shape[0],):
            raise ValueError("Argument output should be a {0} array.".format((points.shape[0],)))
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)
        # compute kinetic energy
        if index is None:
            # include all orbitals
            self._iodata.obasis.compute_grid_kinetic_dm(dm, points, output=output)
        else:
            # include specified set of orbitals
            raise NotImplementedError()
        return output

    def compute_megga(self, points, spin='ab', index=None):
        """Return electron density, gradient, laplacian & kinetic energy density.

        Parameters
        ----------
        points : ndarray
           The 2d-array containing the cartesian coordinates of points on which density is
           evaluated. It has a shape (n, 3) where n is the number of points.
        spin : str
           The type of occupied spin orbitals. By default, the alpha and beta electrons (i.e.
           alpha and beta occupied spin orbitals) are used for computing the electron density.

           - "a" or "alpha": consider alpha electrons
           - "b" or "beta": consider beta electrons
           - "ab": consider alpha and beta electrons

        index : sequence
           Sequence of integers representing the index of spin orbitals. Alpha and beta spin
           orbitals are each indexed from 1 to :attr:`nbasis`.
           If ``None``, all occupied spin orbitals are included.
        """
        # check points
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2d-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2d-array of floats!")

        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin)

        # compute for the given set of orbitals
        if index is None:
            output = self._iodata.obasis.compute_grid_mgga_dm(dm, points)
        else:
            raise NotImplementedError()
        return output[:, 0], output[:, 1:4], output[:, 4], output[:, 5]
