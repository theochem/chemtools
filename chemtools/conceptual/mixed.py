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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on Mixed Energy Models."""


from chemtools.conceptual import LinearGlobalTool


__all__ = ["MixedGlobalTool"]


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
            \mu^+_{\text{GCV}} = -\frac{I + 3A}{4}
            \mu^-_{\text{GCV}} = -\frac{3I + A}{4}

        Returns
        -------
        mu_p : float
            Chemical potential from above, :math:`\mu^+_{\text{GCV}}`
        mu_m : float
            Chemical potential from below, :math:`\mu^-_{\text{GCV}}`
        """
        mu_p = - (self.ip + 3 * self.ea) / 4.
        mu_m = - (3 * self.ip + self.ea) / 4.
        return mu_p, mu_m

    @property
    def electron_transfer_power_gcv(self):
        r"""Global Electron accepting & electron donating powers of Gazquez, Cedillo & Vela.

         Equation [10] & [9] of J. Phys. Chem. A (2007) 111, 1966-1970:

         .. math::
            \omega^+_{\text{GCV}} = -\frac{(I + 3A)^2}{32(I - A)}
            \omega^-_{\text{GCV}} = -\frac{(3I + A)^2}{32(I - A)} \\

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
           \mu_{\text{base}}(\alpha) &= -\frac{I + \alpha A}{1 + \alpha} \\

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
