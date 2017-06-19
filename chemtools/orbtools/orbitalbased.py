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

    def __init__(self, points, obasis, exp_alpha, exp_beta=None):
        r"""
        Initialize OrbitalLocalTool class using gridpoints, basisset and orbital expansion.

        Parameters
        ----------
        points : np.ndarray
            Gridpoints used to calculate the properties.
        obasis :
            An instance of `GOBasis` class from `HORTON` library representing the orbital basis.
        exp_alpha :
            An expansion of the alpha orbitals in a basis set,
            with orbital energies and occupation numbers.
        exp_beta : default=None
            An expansion of the beta orbitals in a basis set,
            with orbital energies and occupation numbers.
        """
        self._obasis = obasis
        self._exp_alpha = exp_alpha

        if exp_beta is None:
            self._dm = self._exp_alpha.to_dm(factor=2.0)
        else:
            self._exp_beta = exp_beta
            self._dm = self._exp_alpha.to_dm(factor=1.0)
            self._exp_beta.to_dm(self._dm, factor=1.0, clear=False)

        self._points = points

        # Compute density & gradient on grid
        dens = self._obasis.compute_grid_density_dm(self._dm, self._points)
        grad = self._obasis.compute_grid_gradient_dm(self._dm, self._points)
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
        return self._obasis.compute_grid_kinetic_dm(self._dm, self._points)

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
        # compute mep with HORTON
        return self._obasis.compute_grid_esp_dm(self._dm, coordinates, pseudo_numbers, self._points)

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
        # Our orbital numbering starts at 1, but HORTON starts at 0.
        iorbs -= 1
        if iorbs.ndim == 0:
            iorbs = np.array([iorbs])
        if spin == 'alpha':
            return self._obasis.compute_grid_orbitals_exp(self._exp_alpha, self._points, iorbs)
        elif spin == 'beta':
            return self._obasis.compute_grid_orbitals_exp(self._exp_beta, self._points, iorbs)
        else:
            raise ValueError('Argument spin={0} is not known.'.format(spin) +
                             'Choose either "alpha" or "beta".')

    @property
    def local_ip(self):
        r"""
        Local Ionization Potential defined as,

        .. math::
           IP \left(\mathbf{r}\right) = \frac{\sum_{i \in \mathrm{MOs}} n_i \epsilon_i
           \phi_i(\mathbf{r}) \phi_i^*(\mathbf{r})}{\rho(\mathbf{r})}
        """
        iorbs = np.arange(1, self._obasis.nbasis+1)
        if hasattr(self, '_exp_beta'):
            orbitals_alpha = self.orbitals_exp(iorbs, spin='alpha')**2
            orbitals_beta = self.orbitals_exp(iorbs, spin='beta')**2
            result = np.dot(self._exp_alpha.occupations*self._exp_alpha.energies, orbitals_alpha.T)
            result += np.dot(self._exp_beta.occupations*self._exp_beta.energies, orbitals_beta.T)
        else:
            orbitals = self.orbitals_exp(iorbs)**2
            result = np.dot(2.0*self._exp_alpha.occupations*self._exp_alpha.energies, orbitals.T)
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

        e_alpha = np.array(self._exp_alpha.energies, dtype=np.float128, copy=True)
        n_alpha = np.sum(self._exp_alpha.occupations)

        spin_pot_a = bisect(lambda x: _spin_chemical_function(x, n_alpha, bt, e_alpha),
                            e_alpha[0], e_alpha[-1], maxiter=maxiter)

        if hasattr(self, '_exp_beta'):
            e_beta = np.array(self._exp_beta.energies, dtype=np.float128)
            n_beta = np.sum(self._exp_beta.occupations)
            spin_pot_b = bisect(lambda x: _spin_chemical_function(x, n_beta, bt, e_beta),
                                e_beta[0], e_beta[-1], maxiter=maxiter)

            return np.array([spin_pot_a, spin_pot_b])

        else:
            return np.array([spin_pot_a, spin_pot_a])

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
        nbf = self._exp_alpha.energies.shape[0]
        tempdens = np.zeros(self._points.shape[0])
        iorbs = np.arange(1, self._obasis.nbasis+1)

        if spin_chemical_potential is None:
            spin_chemical_potential = self.spin_chemical_potential(temperature)

        orbdens = self.orbitals_exp(iorbs, spin='alpha')**2.

        for i in range(0, nbf):
            denom = (1. + np.exp(bt * (self._exp_alpha.energies[i] - spin_chemical_potential[0])))
            tempdens[:] += orbdens[:, i] / denom

        if hasattr(self, '_exp_beta'):
            orbdens = self.orbitals_exp(iorbs, spin='beta')**2.

            for i in range(0, nbf):
                denom = (1. + np.exp(bt * (self._exp_beta.energies[i] -
                                           spin_chemical_potential[1])))
                tempdens[:] += orbdens[:, i] / denom
            return tempdens

        else:
            return 2.0*tempdens


def _spin_chemical_function(x, n, bt, en):
    r"""
    Function for calculating the spin chemical potential.

    The function:
    .. math::
        f(\mu_{sigma}, \beta, N_{\sigma}) =  \sum_{i = 1}^{N_{basis}} \frac{1}
        {1 + e^{\beta(\epsilon_{i \sigma} - \mu_{sigma})}} - N_{\sigma}

    Parameters
    ----------

    x : float
        \mu_{sigma}, the Spin Chemical Potential
    bt : float
        :math:`\beta = \frac{1}{k_B T}`
    n : float
        the number of electrons, :math:`N_{\sigma}`, with spin :math:`\sigma = \{ \alpha, \beta\}`
    en : np.array
        :math:`\epsilon_{i \sigma}` the molecular orbital energies

    Returns
    -------
    float
        The value of :math:`f(\mu_{sigma}, N)`.
    """
    return np.sum(1. / (1. + np.exp(bt * (en - x)))) - n
