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
# pragma pylint: disable=invalid-name
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Mixed Energy Models."""


from chemtools.conceptual import LinearGlobalTool, LinearLocalTool, LinearCondensedTool
from chemtools.conceptual import QuadraticGlobalTool, QuadraticLocalTool, QuadraticCondensedTool


__all__ = ["MixedGlobalTool", "MixedLocalTool", "MixedCondensedTool"]


class MixedGlobalTool(object):
    """Class of global conceptual DFT reactivity descriptors based on mixed energy models."""

    def __init__(self, dict_energy):
        r"""Initialize mixed energy model to compute global reactivity descriptors.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): E(N_0 - 1), N_0: E(N_0),
            (N_0 + 1): E(N_0 + 1)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        """
        model = LinearGlobalTool(dict_energy)
        self.ip, self.ea = model.ip, model.ea

    @property
    def chemical_potential_gcv(self):
        r"""Global chemical potential definition of Gazquez, Cedillo & Vela.

         Equation [11] of J. Phys. Chem. A (2007) 111, 1966-1970:

         .. math::
            \mu^+_{\text{GCV}} = -\frac{I + 3A}{4} \\
            \mu^-_{\text{GCV}} = -\frac{3I + A}{4}

        Returns
        -------
        mu_p : float
            Chemical potential from above, :math:`\mu^+_{\text{GCV}}`.
        mu_m : float
            Chemical potential from below, :math:`\mu^-_{\text{GCV}}`.
        """
        mu_p = - (self.ip + 3 * self.ea) / 4.
        mu_m = - (3 * self.ip + self.ea) / 4.
        return mu_p, mu_m

    @property
    def electron_transfer_power_gcv(self):
        r"""Global Electron accepting & electron donating powers of Gazquez, Cedillo & Vela.

         Equation [10] & [9] of J. Phys. Chem. A (2007) 111, 1966-1970:

         .. math::
            \omega^+_{\text{GCV}} = -\frac{(I + 3A)^2}{32(I - A)} \\
            \omega^-_{\text{GCV}} = -\frac{(3I + A)^2}{32(I - A)}

        Returns
        -------
        omega_p : float
            Electron accepting power, :math:`\omega^+_{\text{GCV}}`.
        omega_m : float
            Electron donating power, :math:`\omega^-_{\text{GCV}}`.
        """
        omega_p = (self.ip + 3 * self.ea)**2 / (32. * (self.ip - self.ea))
        omega_m = (3 * self.ip + self.ea)**2 / (32. * (self.ip - self.ea))
        return omega_p, omega_m

    def chemical_potential_ma(self, alpha):
        r"""Adjusted Chemical potential of Miranda-Quintana and Ayers.

        Equation [65] & [64] of Phys. Chem. Chem. Phys. (2016) 18, 15070-15080:

        .. math::
           \mu_{\text{acid}}(\alpha) &= -\frac{\alpha I + A}{1 + \alpha} \\
           \mu_{\text{base}}(\alpha) &= -\frac{I + \alpha A}{1 + \alpha}

        Note: For :math:`\alpha = 3`, this equals to the :attr:`chemical_potential_gcv`.

        Parameters
        ----------
        alpha : float
            The alpha parameter, :math:`\alpha`.

        Returns
        -------
        mu_a : float
            Chemical potential of acid, :math:`\mu_{\text{acid}}(\alpha)`.
        mu_b : float
            Chemical potential of base, :math:`\mu_{\text{base}}(\alpha)`.
        """
        # check alpha

        # compute chemical potential
        mu_a = - (alpha * self.ip + self.ea) / (1. + alpha)
        mu_b = - (self.ip + alpha * self.ea) / (1. + alpha)
        return mu_a, mu_b


class MixedLocalTool(object):
    """Class of local conceptual DFT reactivity descriptors based on mixed energy models."""

    def __init__(self, dict_energy, dict_density):
        r"""Initialize to compute mixed local reactivity descriptors.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): E(N_0 - 1), N_0: E(N_0),
            (N_0 + 1): E(N_0 + 1)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        dict_density : dict
            Dictionary of number of electrons (keys) and corresponding density array (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): \rho_{N_0 - 1}\left(\mathbf{
            r}\right), N_0: \rho_{N_0}\left(\mathbf{r}\right), (N_0 + 1): \rho_{N_0 + 1}\left(
            \mathbf{r}\right)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        """
        # check matching number of electrons
        if sorted(dict_energy.keys()) != sorted(dict_density.keys()):
            nums_e, nums_d = sorted(dict_energy.keys()), sorted(dict_density.keys())
            raise ValueError("The number of electrons (keys) in dict_energy and dict_density "
                             "arguments should match! {0} != {1}".format(nums_e, nums_d))
        # make quadratic global and local classes
        self.quad_g = QuadraticGlobalTool(dict_energy)
        n_max, softness = self.quad_g.n_max, self.quad_g.softness
        self.quad_l = QuadraticLocalTool(dict_density, n_max, softness)
        # make linear global and local classes
        self.lin_g = LinearGlobalTool(dict_energy)
        n_max, softness = self.lin_g.n_max, self.lin_g.softness
        self.lin_l = LinearLocalTool(dict_density, n_max, softness)

    @property
    def softness_yp(self):
        r"""Local softness of Yang and Parr.

        Equation [18] of Proc. Natl. Acad. Sci. USA (1985) 82, 6723-6726:

        .. math::
           s^+(\mathbf{r}) &= S f^+(\mathbf{r}) \\
           s^0(\mathbf{r}) &= S f^0(\mathbf{r}) \\
           s^-(\mathbf{r}) &= S f^-(\mathbf{r})

        where :math:`f^{+,0,-}(\mathbf{r})` is Fukui function from the linear energy model,
        and :math:`S={}^1/_{\eta}` is global chemical softness (inverse of global chemical
        hardness) from the quadratic energy model.

        Returns
        -------
        softness_p : ndarray
            Local softness from above measuring nucleophilic attack, :math:`s^+(\mathbf{r})`.
        softness_0 : ndarray
            Local softness (centered) measuring radical attack, :math:`s^0(\mathbf{r})`.
        softness_m : ndarray
            Local softness from below measuring electrophilic attack, :math:`s^-(\mathbf{r})`.
        """
        softness_p = self.lin_l.ff_plus / self.quad_g.chemical_hardness
        softness_0 = self.lin_l.ff_zero / self.quad_g.chemical_hardness
        softness_m = self.lin_l.ff_minus / self.quad_g.chemical_hardness
        return softness_p, softness_0, softness_m

    @property
    def philicity_mgvgc(self):
        r"""Local philicity measure of Morell, Gazquez, Vela, Guegana & Chermette.

        Equation [46], [15] & [47] of Phys. Chem. Chem. Phys. (2014) 16, 26832-26842:

        .. math::
           \omega^+(\mathbf{r}) &= -(\frac{\mu^+}{\eta}) f^+(\mathbf{r}) +
                            \frac{1}{2} \left(\frac{\mu^+}{\eta}\right)^2 f^{(2)}(\mathbf{r}) \\
           \omega^0(\mathbf{r}) &= -(\frac{\mu^0}{\eta}) f^0(\mathbf{r}) +
                            \frac{1}{2} \left(\frac{\mu^0}{\eta}\right)^2 f^{(2)}(\mathbf{r}) \\
           \omega^-(\mathbf{r}) &= +(\frac{\mu^-}{\eta}) f^-(\mathbf{r}) +
                            \frac{1}{2} \left(\frac{\mu^-}{\eta}\right)^2 f^{(2)}(\mathbf{r})

        where :math:`\mu^{+,0,-}` is global chemical potential from the linear energy model,
        :math:`\eta` is global chemical hardness from the quadratic energy model,
        :math:`f^{+,0,-}(\mathbf{r})` is Fukui function from the linear density model, and
        :math:`f^{(2)}(\mathbf{r})` is dual descriptor from the quadratic density model.

        Returns
        -------
        omega_p : ndarray
            Local philicity index from above measuring nucleophilic attack,
            :math:`\omega^+(\mathbf{r})`.
        omega_0 : ndarray
            Local philicity index (centered) measuring radical attack,
            :math:`\omega^0(\mathbf{r})`.
        omega_m : ndarray
            Local philicity index from below measuring electrophilic attack,
            :math:`\omega^-(\mathbf{r})`.
        """
        coeff_p = -1. * self.lin_g.mu_plus / self.quad_g.eta
        omega_p = coeff_p * self.lin_l.ff_plus + 0.5 * coeff_p**2 * self.quad_l.dual_descriptor
        coeff_0 = -1. * self.lin_g.mu_zero / self.quad_g.eta
        omega_0 = coeff_0 * self.lin_l.ff_zero + 0.5 * coeff_0**2 * self.quad_l.dual_descriptor
        coeff_m = self.lin_g.mu_minus / self.quad_g.eta
        omega_m = coeff_m * self.lin_l.ff_minus + 0.5 * coeff_m**2 * self.quad_l.dual_descriptor
        return omega_p, omega_0, omega_m

    @property
    def philicity_cms(self):
        r"""Local philicity index of Chattaraj, Maiti & Sarkar.

        Equation [12] of J. Phys. Chem. A (2003) 107, 4973–4975:

        .. math::
           \omega^+(\mathbf{r}) &= \omega \text{ } f^+(\mathbf{r}) \\
           \omega^0(\mathbf{r}) &= \omega \text{ } f^0(\mathbf{r}) \\
           \omega^-(\mathbf{r}) &= \omega \text{ } f^-(\mathbf{r})

        where :math:`\omega` is global electrophilicity from quadratic energy model, and
        :math:`f^{+,0,-}(\mathbf{r})` is Fukui function from linear energy model.

        Returns
        -------
        omega_p : ndarray
            Local philicity index from above measuring nucleophilic attack,
            :math:`\omega^+(\mathbf{r})`.
        omega_0 : ndarray
            Local philicity index (centered) measuring radical attack,
            :math:`\omega^0(\mathbf{r})`.
        omega_m : ndarray
            Local philicity index from below measuring electrophilic attack,
            :math:`\omega^-(\mathbf{r})`.
        """
        omega_p = self.quad_g.electrophilicity * self.lin_l.ff_plus
        omega_0 = self.quad_g.electrophilicity * self.lin_l.ff_zero
        omega_m = self.quad_g.electrophilicity * self.lin_l.ff_minus
        return omega_p, omega_0, omega_m


class MixedCondensedTool(object):
    """Class of condensed conceptual DFT reactivity descriptors based on mixed energy models."""

    def __init__(self, dict_energy, dict_population):
        r"""Initialize to compute mixed condensed reactivity descriptors.

        Parameters
        ----------
        dict_energy : dict
            Dictionary of number of electrons (keys) and corresponding energy (values).
            This model expects three energy values corresponding to three consecutive number of
            electrons differing by one, i.e. :math:`\{(N_0 - 1): E(N_0 - 1), N_0: E(N_0),
            (N_0 + 1): E(N_0 + 1)\}`. The :math:`N_0` value is considered as the reference number
            of electrons.
        dict_population : dict
            Dictionary of number of electrons (keys) and corresponding atomic population array
            (values). This model expects three atomic populations corresponding to three consecutive
            number of electrons differing by one, i.e. :math:`\{(N_0 - 1): \{p^{(N_0 - 1)}_A\},
            N_0: \{p^{(N_0)}_A\}, N_0 + 1: \{p^{(N_0 + 1)}_A\}\}`.
            The :math:`N_0` value is considered as the reference number of electrons.
        """
        # check matching number of electrons
        if sorted(dict_energy.keys()) != sorted(dict_population.keys()):
            nums_e, nums_d = sorted(dict_energy.keys()), sorted(dict_population.keys())
            raise ValueError("The number of electrons (keys) in dict_energy and dict_population "
                             "arguments should match! {0} != {1}".format(nums_e, nums_d))
        # make quadratic global and condensed classes
        self.quad_g = QuadraticGlobalTool(dict_energy)
        n_max, softness = self.quad_g.n_max, self.quad_g.softness
        self.quad_c = QuadraticCondensedTool(dict_population, n_max, softness)
        # make linear global and condensed classes
        self.lin_g = LinearGlobalTool(dict_energy)
        n_max, softness = self.lin_g.n_max, self.lin_g.softness
        self.lin_c = LinearCondensedTool(dict_population, n_max, softness)

    @property
    def softness_yp(self):
        r"""Condensed softness of Yang and Parr.

        Atom-condensed implementation of local :attr:`MixedLocalTool.softness_yp`:

        .. math::
           s^+_A &= S \text{ } f^+_A \\
           s^0_A &= S \text{ } f^0_A \\
           s^-_A &= S \text{ } f^-_A

        where :math:`f^{+,0,-}_A` is condensed Fukui function from the linear energy model,
        and :math:`S={}^1/_{\eta}` is global chemical softness (inverse of global chemical
        hardness) from the quadratic energy model.

        Returns
        -------
        softness_p : ndarray
            Condensed softness from above measuring nucleophilic attack,
            :math:`\{s^+_A\}_{A=1}^{N_{\text{atoms}}}`.
        softness_0 : ndarray
            Condensed softness (centered) measuring radical attack,
            :math:`\{s^0_A\}_{A=1}^{N_{\text{atoms}}}`.
        softness_m : ndarray
            Condensed softness from below measuring electrophilic attack,
            :math:`\{s^-_A\}_{A=1}^{N_{\text{atoms}}}`.
        """
        softness_p = self.lin_c.ff_plus / self.quad_g.chemical_hardness
        softness_0 = self.lin_c.ff_zero / self.quad_g.chemical_hardness
        softness_m = self.lin_c.ff_minus / self.quad_g.chemical_hardness
        return softness_p, softness_0, softness_m

    @property
    def philicity_mgvgc(self):
        r"""Condensed philicity measure of Morell, Gazquez, Vela, Guegana & Chermette.

        Atom-condensed implementation of local :attr:`MixedLocalTool.philicity_mgvgc`:

        .. math::
           \omega^+_A &= -\left(\frac{\mu^+}{\eta}\right) f^+_A +
                            \frac{1}{2} \left(\frac{\mu^+}{\eta}\right)^2 f^{(2)}_A \\
           \omega^0_A &= -\left(\frac{\mu^0}{\eta}\right) f^0_A +
                            \frac{1}{2} \left(\frac{\mu^0}{\eta}\right)^2 f^{(2)}_A \\
           \omega^-_A &= +\left(\frac{\mu^-}{\eta}\right) f^-_A +
                            \frac{1}{2} \left(\frac{\mu^-}{\eta}\right)^2 f^{(2)}_A

        where :math:`\mu^{+,0,-}` is global chemical potential from the linear energy model,
        :math:`\eta` is global chemical hardness from the quadratic energy model,
        :math:`f^{+,0,-}_A` is condensed Fukui function from the linear energy model, and
        :math:`f^{(2)}_A` is condensed dual descriptor from the quadratic energy model.

        Returns
        -------
        omega_p : ndarray
            Local philicity index from above measuring nucleophilic attack,
            :math:`\{\omega^+_A\}_{A=1}^{N_{\text{atoms}}}`.
        omega_0 : ndarray
            Local philicity index (centered) measuring radical attack,
            :math:`\{\omega^0_A\}_{A=1}^{N_{\text{atoms}}}`.
        omega_m : ndarray
            Local philicity index from below measuring electrophilic attack,
            :math:`\{\omega^-_A\}_{A=1}^{N_{\text{atoms}}}`.
        """
        coeff_p = -1. * self.lin_g.mu_plus / self.quad_g.eta
        omega_p = coeff_p * self.lin_c.ff_plus + 0.5 * coeff_p**2 * self.quad_c.dual_descriptor
        coeff_0 = -1. * self.lin_g.mu_zero / self.quad_g.eta
        omega_0 = coeff_0 * self.lin_c.ff_zero + 0.5 * coeff_0**2 * self.quad_c.dual_descriptor
        coeff_m = self.lin_g.mu_minus / self.quad_g.eta
        omega_m = coeff_m * self.lin_c.ff_minus + 0.5 * coeff_m**2 * self.quad_c.dual_descriptor
        return omega_p, omega_0, omega_m

    @property
    def philicity_cms(self):
        r"""Condensed philicity index of Chattaraj, Maiti & Sarkar.

        Atom-condensed implementation of local :attr:`MixedLocalTool.philicity_cms`:

        .. math::
           \omega^+_A &= \omega \text{ } f^+_A \\
           \omega^0_A &= \omega \text{ } f^0_A \\
           \omega^-_A &= \omega \text{ } f^-_A

        where :math:`f^{+,0,-}_A` is condensed Fukui function from linear energy model, and
        :math:`\omega` is global electrophilicity from quadratic energy model.

        Returns
        -------
        omega_p : ndarray
            Condensed philicity index from above measuring nucleophilic attack,
            :math:`\{\omega^+_A\}_{A=1}^{N_{\text{atoms}}}`.
        omega_0 : ndarray
            Condensed philicity index (centered) measuring radical attack,
            :math:`\{\omega^0_A\}_{A=1}^{N_{\text{atoms}}}`.
        omega_m : ndarray
            Condensed philicity index from below measuring electrophilic attack,
            :math:`\{\omega^-_A\}_{A=1}^{N_{\text{atoms}}}`.
        """
        omega_p = self.quad_g.electrophilicity * self.lin_c.ff_plus
        omega_0 = self.quad_g.electrophilicity * self.lin_c.ff_zero
        omega_m = self.quad_g.electrophilicity * self.lin_c.ff_minus
        return omega_p, omega_0, omega_m

    @property
    def philicity_rkgp(self):
        r"""Relative electrophilicity & nucleophilicity of Roy, Krishnamurti, Geerlings & Pal.

        Based on J. Phys. Chem. A (1998), 102, 3746–3755:

        .. math::
           \epsilon_{\text{electrophilicity}, A} &= \frac{s^+_A}{s^-_A} = \frac{f^+_A}{f^-_A} \\
           \epsilon_{\text{nucleophilicity}, A} &= \frac{s^-_A}{s^+_A} = \frac{f^-_A}{f^+_A}

        Returns
        -------
        epsilon_e : ndarray
            Condensed relative electrophilicity, :math:`\epsilon_{\text{electrophilicity},A}`.
        epsilon_n : ndarray
            Condensed relative nucleophilicity, :math:`\epsilon_{\text{nucleophilicity},A}`.
        """
        softness_plus, _, softness_minus = self.softness_yp
        epsilon_e = softness_plus / softness_minus
        epsilon_n = softness_minus / softness_plus
        return epsilon_e, epsilon_n
