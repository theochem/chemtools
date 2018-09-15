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
# pragma pylint: disable=invalid-name
"""Orbital-Based Local Tools."""


import numpy as np

from scipy.optimize import bisect

from chemtools.wrappers.molecule import Molecule
from chemtools.denstools.densitybased import DensityLocalTool


__all__ = ["OrbitalLocalTool"]


class OrbitalLocalTool(DensityLocalTool):
    """Class of orbital-based descriptive tools."""

    def __init__(self, molecule, points):
        r"""Initialize class using instance of `Molecule` and grid points.

        Parameters
        ----------
        molecule : Molecule
            An instance of `Molecule` class
        points : np.ndarray
            Grid points, given as a 2D array with 3 columns, used for calculating local properties.
        """
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("Argument points should be a 2D array with 3 columns.")

        self._molecule = molecule
        self._points = points
        # compute density, gradient & hessian on grid
        dens = self._molecule.compute_density(self._points)
        grad = self._molecule.compute_gradient(self._points)
        # hess = self._molecule.compute_hessian(self._points)
        super(OrbitalLocalTool, self).__init__(dens, grad, None)

    @classmethod
    def from_file(cls, filename, points):
        """Initialize class from file.

        Parameters
        ----------
        filename : str
            Path to molecule's files.
        points : np.ndarray
            Grid points, given as a 2D array with 3 columns, used for calculating local properties.
        """
        molecule = Molecule.from_file(filename)
        return cls(molecule, points)

    @property
    def kinetic_energy_density(self):
        r"""Positive definite kinetic energy density.

        .. math::
           \tau \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2
        """
        return self._molecule.compute_kinetic_energy_density(self._points)

    @property
    def electron_localization_function(self):
        r"""Electron Localization Function introduced by Becke and Edgecombe.

        .. math::
           ELF (\mathbf{r}) =
                \frac{1}{\left( 1 + \left(\frac{D_{\sigma}(\mathbf{r})}
                {D_{\sigma}^0 (\mathbf{r})} \right)^2\right)}

        with XXX, XXX, and positive definite kinetic energy density defined as, respectively,

        .. math::
            D_{\sigma} (\mathbf{r}) &= \tau_{\sigma} (\mathbf{r}) -
               \frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}}

           D_{\sigma}^0 (\mathbf{r}) &=
              \frac{3}{5} (6 \pi^2)^{2/3} \rho_{\sigma}^{5/3} (\mathbf{r})

           \tau_{\sigma} (\mathbf{r}) =
                 \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2
        """
        elfd = self.kinetic_energy_density - self.weizsacker_kinetic_energy_density
        tf = np.ma.masked_less(self.thomas_fermi_kinetic_energy_density, 1.0e-30)
        tf.filled(1.0e-30)
        elf = 1.0 / (1.0 + (elfd / tf)**2.0)
        return elf

    @property
    def electrostatic_potential(self):
        r"""Molecular Electrostatic Potential.

        .. math::
           V \left(\mathbf{r}\right) = \sum_A \frac{Z_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
             \int \frac{\rho \left(\mathbf{r}"\right)}{\rvert \mathbf{r}" -
             \mathbf{r} \lvert} d\mathbf{r}"
        """
        return self._molecule.compute_esp(self._points)

    @property
    def local_ionization_potential(self):
        r"""Local Ionization Potential.

        .. math::
           IP \left(\mathbf{r}\right) = \frac{\sum_{i \in \mathrm{MOs}} n_i \epsilon_i
           \phi_i(\mathbf{r}) \phi_i^*(\mathbf{r})}{\rho(\mathbf{r})}
        """
        iorbs = np.arange(1, self._molecule.nbasis + 1)

        orbitals_alpha = self.compute_orbital_expression(iorbs, spin="alpha") ** 2
        orbitals_beta = self.compute_orbital_expression(iorbs, spin="beta") ** 2
        result = np.dot(self._molecule.orbital_occupation[0] * self._molecule.orbital_energy[0],
                        orbitals_alpha.T)
        result += np.dot(self._molecule.orbital_occupation[1] * self._molecule.orbital_energy[1],
                         orbitals_beta.T)
        result /= self.density
        return result

    def compute_orbital_expression(self, index, spin="alpha"):
        r"""Compute molecular orbital expression of the specified orbitals on the grid.

        Parameters
        ----------
        index : np.ndarray, int, list, tuple
            The indexes of the orbitals to be computed.
            As is common in chemistry we start the orbital numbering at 1 and
            are not using the python numbering.
        spin : str
            the spin of the orbitals to be calculated.
        """
        index = np.copy(np.asarray(index))
        return self._molecule.compute_molecular_orbital(self._points, spin, index=index)

    def compute_spin_chemical_potential(self, temperature, maxiter=500):
        r"""Compute alpha and beta spin chemical potentials at the given temperature.

        The spin chemical potential, :math:`\mu_{\sigma}`, is obtained by solving the following
        one-dimensional nonlinear equations:

        .. math::
           N_{\sigma} = \sum_{i = 1}^{N_{basis}} \frac{1}
           {1 + e^{\beta(\epsilon_{i \sigma} - \mu_{sigma})}}

        where :math:`\beta = \frac{1}{k_B T}` is the so-called thermodynamic beta,
        :math:`\epsilon_{i \sigma}` is the molecular orbital energy of orbital :math:`i \sigma` and
        :math:`N_{\sigma}` is the number of electrons with spin :math:`\sigma = \{ \alpha, \beta\}`.

        Parameters
        ----------
        temperature : float
            Temperature at which to evaluate the spin chemical potential (in Kelvin).
        maxiter : int, optional
            Maximum number of iterations used in solving the equation with scipy.optimize.bisect.

        Returns
        -------
        spin_mu_a : float
            Alpha spin chemical potential.
        spin_mu_b : float
            Beta spin chemical potential.
        """
        kb = 3.1668144e-6  # Boltzman constant in Hartree/Kelvin
        bt = np.divide(1.0, (kb * temperature))

        e_alpha = np.array(self._molecule.orbital_energy[0], dtype=np.float128, copy=True)
        e_beta = np.array(self._molecule.orbital_energy[1], dtype=np.float128, copy=True)
        n_alpha = np.sum(self._molecule.orbital_occupation[0])
        n_beta = np.sum(self._molecule.orbital_occupation[1])

        spin_pot_a = bisect(lambda x: (np.sum(1. / (1. + np.exp(bt * (e_alpha - x)))) - n_alpha),
                            e_alpha[0], e_alpha[-1], maxiter=maxiter)

        spin_pot_b = bisect(lambda x: (np.sum(1. / (1. + np.exp(bt * (e_beta - x)))) - n_beta),
                            e_beta[0], e_beta[-1], maxiter=maxiter)

        return spin_pot_a, spin_pot_b

    def compute_temperature_dependent_density(self, temperature):
        r"""Compute temperature-dependent density of alpha and beta electrons on the grid.

        Parameters
        ----------
        temperature : float
            Temperature at which to evaluate the spin chemical potential (in Kelvin).

        Returns
        -------
        dens_temp : np.array
            Temperature-dependent density evaluated on the grid points.
        """
        kb = 3.1668144e-6  # Boltzman constant in Hartree/Kelvin
        bt = np.divide(1.0, (kb * temperature))
        nbf = self._molecule.nbasis
        tempdens = np.zeros(self._points.shape[0])
        iorbs = np.arange(1, nbf + 1)

        spin_chemical_potential = self.compute_spin_chemical_potential(temperature)

        orbitals_alpha = self.compute_orbital_expression(iorbs, spin="alpha")**2
        orbitals_beta = self.compute_orbital_expression(iorbs, spin="beta")**2

        for i in range(0, nbf):
            denom = np.exp(bt * (self._molecule.orbital_energy[0][i] - spin_chemical_potential[0]))
            tempdens[:] += orbitals_alpha[:, i] / (1. + denom)
            denom = np.exp(bt * (self._molecule.orbital_energy[1][i] - spin_chemical_potential[1]))
            tempdens[:] += orbitals_beta[:, i] / (1. + denom)

        return tempdens
