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
"""Orbital-Based Local Tools."""


import numpy as np
from scipy.optimize import bisect
from chemtools.denstools.densitybased import DensityLocalTool

__all__ = ['OrbitalLocalTool']


class OrbitalLocalTool(DensityLocalTool):
    """Class of orbital-based descriptive tools."""

    def __init__(self, molecule, points):
        r"""
        Initialize OrbitalLocalTool class using gridpoints, basisset and orbital expansion.

        Parameters
        ----------
        molecule :
            An instance of `Molecule` class
        points : np.ndarray
            Gridpoints used to calculate the properties.
        """
        self._molecule = molecule
        self._points = points

        # Compute density & gradient on grid
        dens = self._molecule.compute_density(self._points)
        grad = self._molecule.compute_gradient(self._points)
        super(OrbitalLocalTool, self).__init__(dens, grad, hessian=None)

    @property
    def kinetic_energy_density(self):
        r"""
        Positive definite kinetic energy density.

        Positive definite kinetic energy density defined as,

        .. math::
           \tau \left(\mathbf{r}\right) =
           \sum_i^N n_i \frac{1}{2} \rvert \nabla \phi_i \left(\mathbf{r}\right) \lvert^2
        """
        return self._molecule.compute_kinetic_energy_density(self._points)

    @property
    def elf(self):
        r"""
        Electron Localization Function.

        The Electron Localization Function introduced by Becke and Edgecombe,

        .. math::
           ELF (\mathbf{r}) =
                \frac{1}{\left( 1 + \left(\frac{D_{\sigma}(\mathbf{r})}
                {D_{\sigma}^0 (\mathbf{r})} \right)^2\right)}

        with

        .. math::
            D_{\sigma} (\mathbf{r}) &= \tau_{\sigma} (\mathbf{r}) -
               \frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}}

           D_{\sigma}^0 (\mathbf{r}) &=
              \frac{3}{5} (6 \pi^2)^{2/3} \rho_{\sigma}^{5/3} (\mathbf{r})

        where :math:`\tau_{\sigma} (\mathbf{r})` is the positive definite kinetic energy density,

        .. math::
           \tau_{\sigma} (\mathbf{r}) =
                 \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2
        """
        elfd = self.kinetic_energy_density - self.weizsacker_kinetic_energy_density
        tf = np.ma.masked_less(self.thomas_fermi_kinetic_energy_density, 1.0e-30)
        tf.filled(1.0e-30)
        elf = 1.0 / (1.0 + (elfd / tf)**2.0)
        return elf

    def mep(self, coordinates, pseudo_numbers):
        r"""
        Molecular Electrostatic Potential.

        Molecular Electrostatic Potential defined as,

        .. math::
           V \left(\mathbf{r}\right) = \sum_A \frac{Z_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
             \int \frac{\rho \left(\mathbf{r}'\right)}{\rvert \mathbf{r}' -
             \mathbf{r} \lvert} d\mathbf{r}'
        """
        return self._molecule.compute_esp(self._points)

    def orbitals_exp(self, iorbs, spin='alpha'):
        r"""
        Compute the orbital expectation value on the grid.

        Parameters
        ----------
        iorbs : np.ndarray, int, list, tuple
            The indexes of the orbitals to be computed.
            As is common in chemistry we start the orbital numbering at 1 and
            are not using the python numbering.
        spin : str
            the spin of the orbitals to be calculated.
        """
        iorbs = np.copy(np.asarray(iorbs))
        return self._molecule.compute_molecular_orbital(self._points, spin, orbital_index=iorbs)

    @property
    def local_ip(self):
        r"""
        Local Ionization Potential defined as,

        .. math::
           IP \left(\mathbf{r}\right) = \frac{\sum_{i \in \mathrm{MOs}} n_i \epsilon_i
           \phi_i(\mathbf{r}) \phi_i^*(\mathbf{r})}{\rho(\mathbf{r})}
        """
        iorbs = np.arange(1, self._molecule.nbasis+1)

        orbitals_alpha = self.orbitals_exp(iorbs, spin='alpha')**2
        orbitals_beta = self.orbitals_exp(iorbs, spin='beta')**2
        result = np.dot(self._molecule.orbital_occupation[0]*self._molecule.orbital_energy[0],
                        orbitals_alpha.T)
        result += np.dot(self._molecule.orbital_occupation[1]*self._molecule.orbital_energy[1],
                         orbitals_beta.T)
        result /= self.density
        return result

    def spin_chemical_potential(self, temperature, maxiter=500):
        r"""
        Spin Chemical Potential.

        Spin Chemical Potential, :math:`\mu_{\sigma}` , found by solving the following
        one-dimensional nonlinear equations:

        .. math::
           N_{\sigma} = \sum_{i = 1}^{N_{basis}} \frac{1}
           {1 + e^{\beta(\epsilon_{i \sigma} - \mu_{sigma})}}

        with :math:`\beta = \frac{1}{k_B T}`, the so-called thermodynamic beta,
        :math:`\epsilon_{i \sigma}` the molecular orbital energies and
        :math:`N_{\sigma}` the number of electrons with spin
        :math:`\sigma = \{ \alpha, \beta\}`

        Parameters
        ----------
        temperature : float
            The temperature at which to evaluate the spin chemical potential (in Kelvin).
        maxiter : int, default=500
            The maximum number of iterations.

        Returns
        -------
        np.array, shape=(2,)
            the spin chemical potential as `[alpha_chemical_potential, beta_chemical_potential]`

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

        return np.array([spin_pot_a, spin_pot_b])

    def temperature_dependent_density(self, temperature, spin_chemical_potential=None):
        r"""
        Temperature-Dependent Density.

        Parameters
        ----------
        temperature : float
            The temperatire at which to evaluate the spin chemical potential (in Kelvin).
        spin_chemical_potential : np.array, shape=(2,), default=None
            The spin chemical potential, when not provided it is calculated.

        Returns
        -------
        np.array
            The Temperature-Dependent Density at the gridpoints.
        """
        kb = 3.1668144e-6  # Boltzman constant in Hartree/Kelvin
        bt = np.divide(1.0, (kb * temperature))
        nbf = self._molecule.nbasis
        tempdens = np.zeros(self._points.shape[0])
        iorbs = np.arange(1, nbf + 1)

        if spin_chemical_potential is None:
            spin_chemical_potential = self.spin_chemical_potential(temperature)

        orbitals_alpha = self.orbitals_exp(iorbs, spin='alpha')**2
        orbitals_beta = self.orbitals_exp(iorbs, spin='beta')**2

        for i in range(0, nbf):
            denom = (1. + np.exp(bt * (self._molecule.orbital_energy[0][i]
                                       - spin_chemical_potential[0])))
            tempdens[:] += orbitals_alpha[:, i] / denom

            denom = (1. + np.exp(bt * (self._molecule.orbital_energy[1][i]
                                       - spin_chemical_potential[1])))
            tempdens[:] += orbitals_beta[:, i] / denom

        return tempdens
