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
        # boltzmann constant in hartree/kelvin
        self._kb = 3.1668144e-6
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

    def compute_spin_chemical_potential(self, temperature, maxiter=500, tolerance=1.e-12):
        r"""Compute temperature-dependent spin of alpha and beta electrons. on the grid.

        Spin chemical potential, :math:`\mu_{\sigma, T}`, at temperature :math:`T` is found by
        solving the following one-dimensional nonlinear equation:

        .. math::
            N_{\sigma} = \sum_{i = 1}^{N_{\text{basis}}} \left( \frac{1}
            {1 + e^{\frac{(\epsilon_{i \sigma} - \mu_{\sigma, T})}{k_{\text{B}} T}}} \right)

        where the :math:`\sigma \in \{ \alpha, \beta\}` denotes the electron spin, the
        :math:`N_{\sigma}` is the number of :math:`\sigma`-spin electrons, the
        :math:`\epsilon_{i \sigma}` specifies the energy of :math:`i^{\text{th}}`
        :math:`\sigma`-molecular orbital, and the :math:`\mu_{\sigma, T}` represents the chemical
        potential of :math:`\sigma`-electrons at temperature :math:`T`.
        The :math:`k_{\text{B}}` is the Boltzmann constant.

        This equation is solved using ``scipy.optimize.bisect`` routine for finding root of a
        function within :math:`[a, b]` interval. The first and last :math:`\sigma`-molecular
        orbital energies have been used as bracketing interval to find :math:`\mu_{\sigma, T}`
        at the given temperature :math:`T`.

        Parameters
        ----------
        temperature : float
            Temperature at which to evaluate the spin chemical potential (in Kelvin).
        maxiter : int, optional
            Maximum number of iterations of ``scipy.optimize.bisect`` routine.
        tolerance : float, optional
            Convergence tolerance of ``scipy.optimize.bisect`` routine.

        Returns
        -------
        spin_mu_a : float
            Alpha spin chemical potential.
        spin_mu_b : float
            Beta spin chemical potential.
        """
        bt = 1.0 / (self._kb * temperature)
        # get number and energy of alpha and beta electrons
        n_a, n_b = self._molecule.nelectrons
        energy_a, energy_b = self._molecule.orbital_energy
        # find spin chemical potential of alpha electrons
        spin_pot_a = bisect(lambda x: np.sum(1. / (1. + np.exp(bt * (energy_a - x)))) - n_a,
                            energy_a[0], energy_a[-1], maxiter=maxiter, xtol=tolerance)
        # find spin chemical potential of beta electrons
        spin_pot_b = bisect(lambda x: (np.sum(1. / (1. + np.exp(bt * (energy_b - x)))) - n_b),
                            energy_b[0], energy_b[-1], maxiter=maxiter, xtol=tolerance)
        return spin_pot_a, spin_pot_b

    def compute_temperature_dependent_density(self, temperature):
        r"""Compute temperature-dependent density of alpha and beta electrons on the grid.

        The temperature-dependent :math:`\sigma`-density at temperature :math:`T` is defined as,

        .. math::
            \rho_{\sigma, T} \left(\mathbf{r}\right) = \sum_{i = 1}^{N_{\text{basis}}}
            \left( \frac{1}{1 + e^{\frac{(\epsilon_{i\sigma} - \mu_{\sigma, T})}{k_{\text{B}} T}}}
            \right) |\phi_i (\mathbf{r})|^2

        where the :math:`\sigma \in \{ \alpha, \beta\}` denotes the electron spin, the
        :math:`\epsilon_{i \sigma}` specifies the energy of :math:`i^{\text{th}}`
        :math:`\sigma`-molecular orbital, the :math:`\mu_{\sigma, T}` represents the temperature-
        dependent spin chemical potential of :math:`\sigma`-electrons, and the
        :math:`\phi_{i\sigma}(\mathbf{r})` denotes the :math:`i^{\text{th}}` :math:`\sigma`-
        molecular orbital. The :math:`k_{\text{B}}` is the Boltzmann constant.

        Parameters
        ----------
        temperature : float
            Temperature at which to evaluate the spin chemical potential (in Kelvin).

        Returns
        -------
        dens_temp : np.array
            Temperature-dependent density evaluated on the grid points.
        """
        bt = 1.0 / (self._kb * temperature)
        # compute spin chemical potential & energies of alpha and beta orbitals
        spin_mu_a, spin_mu_b = self.compute_spin_chemical_potential(temperature)
        energy_a, energy_b = self._molecule.orbital_energy
        # compute density of each alpha and beta orbital on grid points
        index = np.arange(1, self._molecule.nbasis + 1)
        dens_a = self.compute_orbital_expression(index, spin="alpha")**2
        dens_b = self.compute_orbital_expression(index, spin="beta")**2
        # compute temperature-dependent density of alpha and beta orbitals
        dens_a /= (1. + np.exp(bt * (energy_a - spin_mu_a)))
        dens_b /= (1. + np.exp(bt * (energy_b - spin_mu_b)))
        # sum temperature-dependent density of all alpha and beta orbitals
        return np.sum(dens_a + dens_b, axis=1)
