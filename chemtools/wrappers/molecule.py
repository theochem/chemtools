# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2024 The ChemTools Development Team
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

# from horton.gbasis.cext import GOBasis
# from horton.meanfield.orbitals import Orbitals
from iodata import load_one
# from iodata.basis import HORTON2_CONVENTIONS, convert_conventions

from gbasis.wrappers import from_iodata
from gbasis.evals.eval import evaluate_basis
from gbasis.evals.density import evaluate_density, evaluate_density_gradient, \
    evaluate_density_laplacian, evaluate_density_hessian, evaluate_posdef_kinetic_energy_density
from gbasis.evals.electrostatic_potential import electrostatic_potential
from gbasis.integrals.overlap import overlap_integral


try:
    from importlib_resources import path
except ImportError:
    from importlib.resources import path

__all__ = ["Molecule"]


class Molecule(object):
    """Molecule class from IOData and HORTON a packages."""

    def __init__(self, iodata):
        """Initialize class.

        Parameters
        ----------
        iodata : iodata.IOData
            An instance of IOData object.
        """
        self._iodata = iodata
        self._coordinates = self._iodata.atcoords
        self._numbers = self._iodata.atnums
        self._pseudo_numbers = self._iodata.atcorenums

        #
        if hasattr(self._iodata, "obasis") and hasattr(self._iodata.obasis, "shells"):
            # self._ao = iodata.obasis
            self._ao = AtomicOrbitals(iodata.obasis)
        else:
            self._ao = None

        if hasattr(self._iodata, "mo") and hasattr(self._iodata.mo, "kind"):
            self._mo = MolecularOrbitals.from_molecule(self)


    @classmethod
    def from_file(cls, fname):
        """Initialize class given a file.

        Parameters
        ----------
        fname : str
            Wavefunction file path.

        """
        # load molecule
        logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
        try:
            iodata = load_one(str(fname))
        except IOError as _:
            try:
                with path("chemtools.data.examples", str(fname)) as fname:
                    logging.info("Loading {0}".format(str(fname)))
                    iodata = load_one(str(fname))
            except IOError as error:
                logging.info(error)
        return cls(iodata)

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
    def pseudo_numbers(self):
        """Pseudo-potential core charges."""
        return self._pseudo_numbers

    @property
    def ao(self):
        """Atomic orbital instance."""
        return self._ao

    @property
    def mo(self):
        """Molecular orbital instance."""
        return self._mo

    def compute_dm(self, spin="ab"):
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
        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            if self.mo:
                return self.mo.dm(spin)
            else:
                raise ValueError(f"Molecular Orbitals is {self.mo} Dm is not available without Molecular Orbitals information")

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
        #
        self._check_argument(points)
        # from IOData MolecularOrbital class
        # In case of restricted: shape = (nbasis, norba) = (nbasis, norbb).
        # In case of unrestricted: shape = (nbasis, norba + norbb).
        # In case of generalized: shape = (2 * nbasis, norb), where norb is the
        # total number of orbitals. (optional)

        # Load basis to Gbasis
        basis = from_iodata(self._iodata)
        # NOTE1: Gbasis shape == Horton.shape.T
        # NOTE2: molecular orbitals are return in order
        if spin == "a":
            basis_mo = evaluate_basis(basis, points, transform= self._iodata.mo.coeffsa.T)
        elif spin == "b":
            basis_mo = evaluate_basis(basis, points, transform=self._iodata.mo.coeffsb.T)
        # todo: check: this used to return the sum of a and b but I am not sure for, at least for
        #  restricted, if it makes sense. I think for the option ab what makes sense is to return
        #  both "a" and "b" spin orbitals matrices but not summed
        elif spin == "ab":
            basis_mo = evaluate_basis(basis, points, transform=self._iodata.mo.coeffsa.T)
        else:
            raise ValueError(f"Spin should be a,b or ab.Got{spin}")

        # Filter by index
        if isinstance(index, int):
            return basis_mo[index-1]
        elif isinstance(index, np.ndarray) and index.ndim == 0:
            basis_mo_index = basis_mo[index-1]
            return basis_mo_index
        elif isinstance(index, np.ndarray) and index.ndim == 1:
            basis_mo_index = [mo for i, mo in enumerate(basis_mo) if i + 1 in index]
            basis_mo_index = np.array(basis_mo_index)
            return basis_mo_index
        elif isinstance(index, list):
            basis_mo_index = [mo for i, mo in enumerate(basis_mo) if i + 1 in index]
            basis_mo_index = np.array(basis_mo_index)
            return basis_mo_index
        elif not index:
            return basis_mo
        else:
            raise ValueError(f"index should be integer, list, or numpy array")


    def compute_density(self, points, spin="ab"):
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

        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin=spin)
            density = evaluate_density(dm, basis, points)
            return density


    def compute_gradient(self, points, spin="ab"):
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

        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin=spin)
            gradient = evaluate_density_gradient(dm, basis, points, deriv_type="direct")
            return gradient

    def compute_hessian(self, points, spin="ab"):
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

        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin=spin)
            hessian = evaluate_density_hessian(dm, basis, points, deriv_type="direct")
            return hessian

    def compute_laplacian(self, points, spin="ab"):
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
        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin=spin)
            laplacian = evaluate_density_laplacian(dm, basis, points, deriv_type="direct")
            return laplacian

    # Todo: Include gbasis API options to compute_esp
    def compute_esp(self, points, spin="ab", charges=None, charges_coords=None):
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
            charges_coords = self.coordinates
        elif not isinstance(charges, np.ndarray) or charges.shape != self.numbers.shape:
            raise ValueError(
                f"Argument charges should be a 1d-array with {self.numbers.shape} shape.")
        elif not isinstance(charges_coords, np.ndarray) or charges_coords.shape != self.coordinates.shape:
            raise ValueError(
                f"Argument charges should be a 2d-array with {self.coordinates.shape} shape.")

        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin)
            esp = electrostatic_potential(basis, dm, points, charges_coords, charges)
            return esp

    def compute_ked(self, points, spin="ab"):
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
        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if spin not in ["a", "b", "ab"]:
            raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            dm = self.compute_dm(spin)
            ked = evaluate_posdef_kinetic_energy_density(dm, basis, points)
            return ked

    def compute_overlap(self, type_ovlp="atomic", spin="ab"):
        r"""Return overlap matrix :math:`\mathbf{S}` of molecular or atomic orbitals.

                .. math:: [\mathbf{S}]_{ij} = \int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}
                """
        #
        # Load basis to Gbasis
        basis = from_iodata(self._iodata)

        if type_ovlp == "atomic":
            transform = None
        elif type_ovlp == "molecular":
            if spin == "a":
                transform = self.coeffs_a.T
            elif spin == "b":
                transform = self.coeffs_b.T
            elif spin == "ab":
                transform = self.coeffs_a.T
            else:
                raise ValueError(f"Specify spin a, b or ab. Got {spin}")
        else:
            raise ValueError(f"Type of overlap must be either atomic or molecular. Got {type_ovlp}")


        # todo: same situation as compute_molecular_orbital. For now it returns overlap using
        #  coeffs a as transformation matrix because for restricted coffsa == coeffsb
        ovlp = overlap_integral(basis, transform=transform)

        return ovlp

    def _check_argument(self, points):
        """Check given arguments.

        Parameters
        ----------
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.

        """
        if not isinstance(points, np.ndarray):
            raise ValueError(f"Argument points should be of type array. Got {type(points)}")
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D-array with 3 columns.")
        if not np.issubdtype(points.dtype, np.float64):
            raise ValueError("Argument points should be a 2D-array of floats!")
        if self._ao is None:
            raise AttributeError("Atomic Orbitals information is needed!")
        if self._mo is None:
            raise AttributeError("Molecular Orbitals information is needed!")


class MolecularOrbitals(object):
    """Molecular orbital class."""

    def __init__(self, occs_a, occs_b, energy_a, energy_b, coeffs_a, coeffs_b):
        self._occs_a, self._occs_b = occs_a, occs_b
        self._energy_a, self._energy_b = energy_a, energy_b
        self._coeffs_a, self._coeffs_b = coeffs_a, coeffs_b

    @classmethod
    def from_molecule(cls, mol):
        """Initialize class given an instance of `Molecule`.

        Parameters
        ----------
        mol : Molecule
            An instance of `Molecule` class.

        """
        # from IOData MolecularOrbital class
        # In case of restricted: shape = (nbasis, norba) = (nbasis, norbb).
        # In case of unrestricted: shape = (nbasis, norba + norbb).
        # In case of generalized: shape = (2 * nbasis, norb), where norb is the
        # total number of orbitals. (optional)
        if mol._iodata.mo.kind == "generalized":
            raise NotImplementedError(
                f"The generalized MOs are not supported, got {mol._iodata.mo.kind}"
            )

        occs_a = mol._iodata.mo.occsa
        occs_b = mol._iodata.mo.occsb
        energy_a = mol._iodata.mo.energiesa
        energy_b = mol._iodata.mo.energiesb
        coeffs_a = mol._iodata.mo.coeffsa
        coeffs_b = mol._iodata.mo.coeffsb

        return cls(occs_a, occs_b, energy_a, energy_b, coeffs_a, coeffs_b)

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
    def homo_index(self):
        """Index of alpha and beta HOMO orbital."""
        # HORTON indexes the orbitals from 0, so 1 is added to get the intuitive index
        index_a = np.argwhere(self._occs_a == 0.0)[0, 0]
        index_b = np.argwhere(self._occs_b == 0.0)[0, 0]
        return index_a, index_b

    @property
    def lumo_index(self):
        """Index of alpha and beta LUMO orbital."""
        # HORTON indexes the orbitals from 0, so 1 is added to get the intuitive index
        return self.homo_index[0] + 1, self.homo_index[1] + 1

    @property
    def homo_energy(self):
        """Energy of alpha and beta HOMO orbital."""
        return self._energy_a[self.homo_index[0] - 1], self._energy_b[self.homo_index[1] - 1]

    @property
    def lumo_energy(self):
        """Energy of alpha and beta LUMO orbital."""
        return self._energy_a[self.lumo_index[0] - 1], self._energy_b[self.lumo_index[1] - 1]

    @property
    def occupation(self):
        """Orbital occupation of alpha and beta electrons."""
        return self._occs_a, self._occs_b

    @property
    def energy(self):
        """Orbital energy of alpha and beta electrons."""
        return self._energy_a, self._energy_b

    @property
    def coefficient(self):
        """Orbital coefficient of alpha and beta electrons.

        The alpha and beta orbital coefficients are each storied in a 2d-array in which
        the columns represent the basis coefficients of each molecular orbital.
        """
        return self._coeffs_a, self._coeffs_b

    @property
    def nelectrons(self):
        """Number of alpha and beta electrons."""
        return np.sum(self._occs_a), np.sum(self._occs_b)

    def compute_overlap(self):
        return NotImplementedError

    def dm(self, spin="ab"):
        """Returns density matrix object corresponding to the specified spin orbitals.

        Parameters
        ----------
        spin : str, optional
           Type of occupied spin orbitals which can be either "a" (for alpha), "b" (for
           beta), and "ab" (for alpha + beta).

        """

        if spin == "ab":
            return self.dm("a") + self.dm("b")
        elif spin == "a":
            arr = np.dot(self._coeffs_a * self._occs_a, self._coeffs_a.T)
        elif spin == "b":
            arr = np.dot(self._coeffs_b * self._occs_b, self._coeffs_b.T)
        else:
            raise ValueError("Argument spin={0} is not recognized!".format(spin))

        return arr


class AtomicOrbitals(object):
    """Gaussian Basis Functions or Atomic Orbitals."""

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

    @property
    def center_index(self):
        # print('puta bidaaaa')
        numcontr_per_shell = [int(shell.nbasis) for shell in self._basis.shells]
        ind_shell_atom = [int(shell.icenter) for shell in self._basis.shells]
        # Duplicate each shell information by the number of contractions in shell
        ind_basis_center = np.repeat(ind_shell_atom, numcontr_per_shell)
        return np.array(ind_basis_center)

    def compute_overlap(self):
        r"""Return overlap matrix :math:`\mathbf{S}` of atomic orbitals.

        .. math:: [\mathbf{S}]_{ij} = \int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}

        """
        # make linear algebra factory
        # compute overlap matrix
        arr = self._basis.compute_overlap()
        return arr

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
        # todo: this functionality is already coded in the compute_gradient in Molecule class
        #  _basis == ao
        """Return gradient of the electron density evaluated on a set of points.

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
            hess[index][np.diag_indices(3)] /= 2.0
        return hess

    def compute_esp(self, dm, points, coordinates, charges):
        """Return electrostatic potential evaluated on the set of points.

        Parameters
        ----------
        dm : ndarray
           First order reduced density matrix of B basis sets given as a 2D array of (B, B) shape.
        points : ndarray
           Cartesian coordinates of N points given as a 2D-array with (N, 3) shape.
        coordinates: ndarray
            Grid to compute numerical calculation of electrostatic potential
        charges: ndarray
            N point charges to use for computing electrostatic potential with electron density

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
