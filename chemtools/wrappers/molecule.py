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

    def __init__(self, iodata, wavefunction=False):
        """
        Initialize class.

        Parameters
        ----------
        iodata : horton.IOData
           An instance of horton.IOData object.
        """
        self._iodata = iodata
        if hasattr(self._iodata, "obasis"):
            self._ao = AtomicOrbitals(self._iodata.obasis)

        self._coordinates = self._iodata.coordinates
        self._numbers = self._iodata.numbers

        if hasattr(self._iodata, "exp_alpha"):
            # assign alpha orbital expression
            self._exp_alpha = self._iodata.exp_alpha
            # assign beta orbital expression
            if hasattr(self._iodata, "exp_beta") and self._iodata.exp_beta is not None:
                self._exp_beta = self._iodata.exp_beta
            else:
                self._exp_beta = self._iodata.exp_alpha
        else:
            self._exp_alpha = None
            self._exp_beta = None
            if wavefunction:
                raise ValueError("There is no wave-function information!")

        # FIXME: following code is pretty hacky. it will be used for the orbital partitioning code
        # GOBasis object stores basis set to atom and angular momentum mapping
        # by shell and not by contraction. So it needs to be converted
        try:
            ind_shell_atom = np.array(iodata.obasis.shell_map)
            ind_shell_orbtype = np.array(iodata.obasis.shell_types)

            def num_contr(orbtype):
                """Return the number of contractions for the given orbital type.

                Parameters
                ----------
                orbtype : int
                    Horton's orbital type scheme.
                    Positive number corresponds to cartesian type and negative number corresponds to
                    spherical types.

                Returns
                -------
                num_contr : int
                    Number of contractions for the given orbital type.

                """
                if orbtype < 0:
                    return 1 + 2 * abs(orbtype)
                return 1 + orbtype + sum(i * (i + 1) / 2 for i in range(1, orbtype + 1))

            numcontr_per_shell = np.array([num_contr(i) for i in ind_shell_orbtype])
            # Duplicate each shell information by the number of contractions in shell
            self._ind_basis_center = np.repeat(ind_shell_atom, numcontr_per_shell)
            self._ind_basis_orbtype = np.repeat(ind_shell_orbtype, numcontr_per_shell)
        except AttributeError:
            pass

    @classmethod
    def from_file(cls, fname, wavefunction=False):
        """Initialize class given a file.

        Parameters
        ----------
        fname : str
            Path to molecule"s files.

        """
        # load molecule
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
        try:
            iodata = IOData.from_file(str(fname))
        except IOError as _:
            try:
                with path("chemtools.data.examples", str(fname)) as fname:
                    logging.info("Loading {0}".format(str(fname)))
                    iodata = IOData.from_file(str(fname))
            except IOError as error:
                logging.info(error)
        return cls(iodata, wavefunction)

    def __getattr__(self, attr):
        """Return attribute.

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
        return self._ao.nbasis

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

    def compute_density_matrix(self, spin="ab", index=None):
        """Compute the density matrix array for the specified spin orbitals.

        Parameters
        ----------
        spin : str, optional
           The type of occupied spin orbitals. Options are "a" (for alpha), "b" (for beta), and
           "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        # get density matrix corresponding to the specified spin
        dm = self._get_density_matrix(spin, index=index)
        return dm._array

    def _get_density_matrix(self, spin="ab", index=None):
        """Return HORTON density matrix object corresponding to the specified spin orbitals.

        Parameters
        ----------
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        # check orbital spin
        if spin not in ["a", "b", "ab"]:
            raise ValueError("Argument spin={0} is not recognized!".format(spin))
        # compute density matrix
        if spin == "ab":
            # get density matrix of alpha & beta electrons
            dm = self._iodata.get_dm_full()
        else:
            # get orbital expression of specified spin
            spin_type = {"a": "alpha", "b": "beta"}
            exp = getattr(self, "_exp_" + spin_type[spin])
            # get density matrix of specified spin
            dm = exp.to_dm()

        if index is not None:
            # convert to numpy array
            index = np.asarray(index)
            # check
            if index.ndim == 0:
                index = index.reshape(1)
            if index.ndim >= 2:
                raise ValueError("Indices should be given as a one-dimensional numpy array.")
            index -= 1
            if np.any(index < 0):
                raise ValueError(
                    "Indices cannot be less than 1. Note that indices start from 1."
                )
            # select the appropriate indices by hard-coding it into the HORTON's DenseTwoIndex
            # object.
            dm._array = dm._array[index[:, np.newaxis], index[np.newaxis, :]]
        return dm

    def compute_molecular_orbital(self, points, spin="ab", index=None):
        """Return molecular orbitals.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        self._check_argument(points)
        # assign orbital index (HORTON index the orbitals from 0)
        if index is None:
            # include all occupied orbitals of specified spin
            spin_index = {"a": 0, "b": 1}
            index = np.arange(self.homo_index[spin_index[spin]])
        else:
            # include specified set of orbitals
            index = np.copy(np.asarray(index)) - 1
            if index.ndim == 0:
                index = np.array([index])
            if np.any(index < 0):
                raise ValueError("Argument index={0} cannot be less than one!".format(index + 1))

        # get orbital expression of specified spin
        spin_type = {"a": "alpha", "b": "beta"}
        exp = getattr(self, "_exp_" + spin_type[spin])
        return self._ao.compute_orbitals(exp, points, index)

    def compute_density(self, points, spin="ab", index=None):
        r"""Return electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        self._check_argument(points)

        # allocate output array
        output = np.zeros((points.shape[0],), float)

        # compute density
        if index is None:
            # get density matrix corresponding to the specified spin
            dm = self._get_density_matrix(spin)
            # include all orbitals
            output = self._ao.compute_density(dm, points)
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

    def compute_gradient(self, points, spin="ab", index=None):
        r"""Return gradient of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        self._check_argument(points)
        dm = self._get_density_matrix(spin, index=index)
        return self._ao.compute_gradient(dm, points)

    def compute_hessian(self, points, spin="ab", index=None):
        r"""Return hessian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        self._check_argument(points)
        dm = self._get_density_matrix(spin, index=index)
        return self._ao.compute_hessian(dm, points)

    def compute_laplacian(self, points, spin="ab", index=None):
        r"""Return Laplacian of the electron density.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        hess = self.compute_hessian(points, spin, index)
        return np.trace(hess, axis1=1, axis2=2)

    def compute_esp(self, points, spin="ab", index=None, charges=None):
        r"""Return molecular electrostatic potential.

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
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.
        charges : np.ndarray, optional
           Array with shape (n,) representing the point charges at the position of the nuclei.
           When ``None``, the pseudo numbers are used.
        """
        self._check_argument(points)
        # assign point charges
        if charges is None:
            charges = self.pseudo_numbers
        elif not isinstance(charges, np.ndarray) or charges.shape != self.numbers.shape:
            raise ValueError("Argument charges should be a 1d-array "
                             "with {0} shape.".format(self.numbers.shape))
        dm = self._get_density_matrix(spin, index=index)
        return self._ao.compute_esp(dm, points, self.coordinates, charges)

    def compute_ked(self, points, spin="ab", index=None):
        r"""Return positive definite or Lagrangian kinetic energy density.

        .. math::
           \tau_\text{PD} \left(\mathbf{r}\right) =
           \tfrac{1}{2} \sum_i^N n_i \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        self._check_argument(points)
        dm = self._get_density_matrix(spin, index=index)
        return self._ao.compute_ked(dm, points)

    def _check_argument(self, points):
        """Check given arguments.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        if not isinstance(points, np.ndarray) or points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2D-array of floats!")


class AtomicOrbitals(object):
    """Gaussian Basis Set."""

    def __init__(self, basis):
        self._basis = basis

    @classmethod
    def from_molecule(cls, mol):
        """Initialize class given an instance of `Molecule`.

        Parameters
        ----------
        mol : Molecule
            An instance of `Molecule` class.

        """
        basis = mol._iodata.obasis
        return cls(basis)

    @classmethod
    def from_file(cls, fname):
        """Initialize class given a file.

        Parameters
        ----------
        fname : str
            Path to molecule"s files.

        """
        return cls.from_molecule(Molecule.from_file(fname))

    @property
    def nbasis(self):
        """int : number of basis functions."""
        return self._basis.nbasis

    def compute_orbitals(self, dm, points, index):
        """

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        index : sequence of int, optional
           Sequence of integers representing the occupied spin orbitals which are indexed
           from 1 to :attr:`nbasis`. If ``None``, all orbitals of the given spin(s) are included.

        """
        return self._basis.compute_grid_orbitals_exp(dm, points, index)

    def compute_density(self, dm, points):
        """Return electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return self._basis.compute_grid_density_dm(dm, points)

    def compute_gradient(self, dm, points):
        """Return gradient of the electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return self._basis.compute_grid_gradient_dm(dm, points)

    def compute_hessian(self, dm, points):
        """Return hessian of the electron density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        # compute upper triangular elements
        output = self._basis.compute_grid_hessian_dm(dm, points)
        # convert the (n, 6) shape to (n, 3, 3)
        hess = np.zeros((len(points), 9))
        # NOTE: hard coded in the indices of the upper triangular matrix in the flattened form
        # in C ordering. Maybe there is a numpy function that does this. This might fail if the
        # hessian is not in c-ordering
        hess[:, [0, 1, 2, 4, 5, 8]] = output
        hess = hess.reshape(len(points), 3, 3)
        hess += np.transpose(hess, axes=(0, 2, 1))
        for index in range(len(points)):
            hess[index][np.diag_indices(3)] /= 2.
        return hess

    def compute_esp(self, dm, points, coordinates, charges):
        """Return electrostatic potential evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return self._basis.compute_grid_esp_dm(dm, coordinates, charges, points)

    def compute_ked(self, dm, points):
        """Return positive definite kinetic energy density evaluated on the a set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        return self._basis.compute_grid_kinetic_dm(dm, points)
