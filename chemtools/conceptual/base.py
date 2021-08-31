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
"""Global Conceptual Density Functional Theory (DFT) Reactivity Tools.

This module contains various global tool classes corresponding to
linear, quadratic, exponential, general energy models.
"""


import logging
import numpy as np
import sympy as sp

from scipy.optimize import newton


__all__ = ["BaseGlobalTool", "BaseLocalTool", "BaseCondensedTool"]


class BaseGlobalTool(object):
    """Base class of global conceptual DFT reactivity descriptors."""

    def __init__(self, n0, n_max):
        r"""Initialize class.

        Parameters
        ----------
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        n_max : float
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
        """
        if n0 <= 0:
            raise ValueError("Argument n0 should be positive! Given n0={0}".format(n0))
        if n_max is not None and n_max < 0:
            raise ValueError("Argument n_max cannot be negative! Given n_max={0}".format(n_max))
        self._n0 = n0
        self._n_max = n_max
        # calculate ionization potential & electron affinity
        self._ip = self.energy(self._n0 - 1) - self.energy(self._n0)
        self._ea = self.energy(self._n0) - self.energy(self._n0 + 1)

    @property
    def n0(self):
        """Reference number of electrons, i.e. :math:`N_0`."""
        return self._n0

    @property
    def n_max(self):
        r"""
        Maximum number of electrons that the system can accept.

        .. math:: N_{\text{max}} = \underbrace {\min }_N E(N)
        """
        return self._n_max

    @property
    def ionization_potential(self):
        r"""
        Ionization potential (IP) of the :math:`N_0`-electron system.

        .. math:: IP = E\left(N_0 - 1\right) - E\left(N_0\right)
        """
        return self._ip

    @property
    def ip(self):
        """The same as :attr:`ionization_potential`."""
        return self.ionization_potential

    @property
    def electron_affinity(self):
        r"""
        Electron affinity (EA) of the :math:`N_0`-electron system.

        .. math:: EA = E\left(N_0\right) - E\left(N_0 + 1\right)
        """
        return self._ea

    @property
    def ea(self):
        """The same as :attr:`electron_affinity`."""
        return self.electron_affinity

    @property
    def electronegativity(self):
        r"""
        Mulliken electronegativity defined as negative :attr:`chemical_potential`.

        .. math:: \chi_{\text{Mulliken}} = - \mu
        """
        if self.chemical_potential is None:
            return None
        value = -1 * self.chemical_potential
        return value

    @property
    def electrophilicity(self):
        r"""
        Electrophilicity of the :math:`N_0`-electron system.

        .. math::
           \omega_{\text{electrophilicity}} = \text{sgn}\left(N_{\text{max}} - N_0\right)
                                              \times \left(E(N_0) - E(N_{\text{max}})\right)
        """
        if self._n_max is None:
            return None
        sign = np.sign(self._n_max - self._n0)
        value = sign * (self.energy(self._n0) - self.energy(self._n_max))
        return value

    @property
    def nucleofugality(self):
        r"""
        Nucleofugality of the :math:`N_0`-electron system.

        .. math::
           \nu_{\text{nucleofugality}} = \text{sgn}\left(N_0 + 1 - N_{\text{max}}\right)
                                         \times \left(E(N_0 + 1) - E(N_{\text{max}})\right)
        """
        if self._n_max is None:
            return None
        sign = np.sign(self._n0 + 1 - self._n_max)
        value = sign * (self.energy(self._n0 + 1) - self.energy(self._n_max))
        return value

    @property
    def electrofugality(self):
        r"""
        Electrofugalityof the :math:`N_0`-electron system.

        .. math::
           \nu_{\text{electrofugality}} = \text{sgn}\left(N_{\text{max}} - N_0 + 1\right)
                                          \times \left(E(N_0 - 1) - E(N_{\text{max}})\right)
        """
        if self._n_max is None:
            return None
        sign = np.sign(self._n_max - self._n0 + 1)
        value = sign * (self.energy(self._n0 - 1) - self.energy(self._n_max))
        return value

    @property
    def chemical_potential(self):
        r"""
        Chemical potential of the :math:`N_0`-electron system.

        The chemical potential is defined as the first derivative of the energy model w.r.t.
        the number of electrons, at fixed external potential, evaluated at :math:`N_0`,

        .. math::
           \mu = \left. \left(\frac{\partial E}{\partial N}
                        \right)_{v(\mathbf{r})} \right|_{N = N_0}
        """
        value = self.energy_derivative(self._n0, order=1)
        return value

    @property
    def mu(self):
        """The same as :attr:`chemical_potential`."""
        return self.chemical_potential

    @property
    def chemical_hardness(self):
        r"""
        Chemical hardness of the :math:`N_0`-electron system.

        This chemical hardness is defined as the second derivative of the energy model
        w.r.t. the number of electrons, at fixed external potential, evaluated at :math:`N_0`.

        .. math::
           \eta  = \left. \left(\frac{\partial^2 E}{\partial N^2}
                          \right)_{v(\mathbf{r})} \right|_{N = N_0}
        """
        value = self.energy_derivative(self._n0, order=2)
        return value

    @property
    def eta(self):
        """The same as :attr:`chemical_hardness`."""
        return self.chemical_hardness

    @property
    def softness(self):
        r"""
        Chemical softness of the :math:`N_0`-electron system.

        The chemical softness is defined as the second derivative of the grand potential model
        w.r.t the number of electrons, at fixed external potential, evaluated at :math:`N_0`.
        This is equal to the inverse chemical hardness.

        .. math::
           S = - \left. \left(\frac{\partial^2 \Omega}{\partial \mu^2}
                        \right)_{v(\mathbf{r})}\right|_{N = N_0} = \frac{1}{\eta}
        """
        # compute 2nd-order derivative of grand potential w.r.t. mu at N0
        value = self.grand_potential_derivative(self._n0, 2)
        if value is not None:
            value *= -1.0
        return value

    def hyper_hardness(self, order=2):
        r"""
        Return the :math:`n^{\text{th}}`-order hyper-hardness of the :math:`N_0`-electron system.

        The :math:`n^{\text{th}}`-order hyper-hardness is defined as the :math:`(n+1)^{\text{th}}`
        -order derivative, where :math:`n \geq 2`, of the energy model w.r.t the number of
        electrons, at fixed external potential, evaluated at :math:`N_0`.

        .. math::
           \eta^{(n)} = \left. \left(\frac{\partial^{n+1} E}{\partial N^{n+1}}
                        \right)_{v(\mathbf{r})} \right|_{N = N_0} \quad \text{for} \quad n \geq 2

        Parameters
        ----------
        order : int, default=2
            The order of hyper-hardness denoted by :math:`n \geq 2` in the formula.
        """
        if not (isinstance(order, int) and order >= 2):
            raise ValueError('Argument order should be an integer greater than or equal to 2.')
        value = self.energy_derivative(self._n0, order + 1)
        return value

    def hyper_softness(self, order):
        r"""
        Return the :math:`n^{\text{th}}`-order hyper-softness of the :math:`N_0`-electron system.

        The :math:`n^{\text{th}}`-order hyper softness is defined as the :math:`(n+1)^{\text{th}}`
        -order derivative, where :math:`n \geq 2`, of the grand potential model w.r.t the number of
        electrons at fixed external potential evaluated at :math:`N_0`.

        .. math::
           S^{(n)} = - \left. \left(\frac{\partial^{n+1} \Omega}{\partial \mu^{n+1}}
                       \right)_{v(\mathbf{r})} \right|_{N = N_0} \quad \text{for} \quad n \geq 2

        Parameters
        ----------
        order : int, default=2
            The order of hyper-hardness denoted by :math:`n \geq 2` in the formula.
        """
        if not (isinstance(order, int) and order >= 2):
            raise ValueError('Argument order should be an integer greater than or equal to 2.')
        # compute derivative of grand potential w.r.t. mu at N0
        value = self.grand_potential_derivative(self._n0, order + 1)
        if value is not None:
            value *= -1.0
        return value

    def energy(self, n_elec):
        r"""
        Return the energy model :math:`E(N)` evaluated for the specified number of electrons.

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        """
        raise NotImplementedError

    def energy_derivative(self, n_elec, order=1):
        r"""
        Return the :math:`n^{\text{th}}`-order derivative of energy w.r.t. the number of electrons.

        This returns the :math:`n^{\text{th}}`-order derivative of energy model :math:`E(N)` w.r.t.
        to the number of electrons, at fixed chemical potential, evaluated for the specified number
        of electrons.

        .. math::
           \left. \left(\frac{\partial^n E}{\partial N^n}
                  \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{\text{elec}} = N_0` the first, second and higher order derivatives are equal
        to the :attr:`BaseGlobalTool.chemical_potential`, :attr:`BaseGlobalTool.chemical_hardness`
        and :attr:`BaseGlobalTool.hyper_hardness`, respectively.
        """
        raise NotImplementedError

    def grand_potential(self, n_elec):
        r"""
        Return the grand potential model evaluated for the specified number of electrons.

        .. math::
           \Omega [\mu(N_{\text{elec}}); v(\mathbf{r})] &=
            E(N_{\text{elec}}) - \mu(N_{\text{elec}}) \times N_{\text{elec}} \\
            &= E(N_{\text{elec}}) - \left.\left(\frac{\partial E(N)}{\partial N}
                    \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}} \times N_{\text{elec}}

        Parameters
        ----------
        n_elec : float
            Number of electrons, :math:`N_{\text{elec}}`.
        """
        if n_elec is None or self.energy_derivative(n_elec, 1) is None:
            return None
        # compute grand potential as a function of N
        value = self.energy(n_elec) - self.energy_derivative(n_elec, 1) * n_elec
        return value

    def grand_potential_derivative(self, n_elec, order=1):
        r"""
        Evaluate the :math:`n^{\text{th}}`-order derivative of grand potential at the given n_elec.

        This returns the :math:`n^{\text{th}}`-order derivative of grand potential model w.r.t.
        to the chemical potential, at fixed external potential, evaluated for the specified
        number of electrons :math:`N_{\text{elec}}`.

        That is,

        .. math::
             \left. \left(\frac{\partial^n \Omega}{\partial \mu^n}
                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}} =
             \left. \left(\frac{\partial^{n-1}}{\partial \mu^{n-1}}
                          \frac{\partial \Omega}{\partial \mu}
                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}} =
           - \left. \left(\frac{\partial^{n-1} N}{\partial \mu^{n-1}}
                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}}
             \quad  n = 1, 2, \dots

        These derivatives can be computed using the derivative of energy model w.r.t. number of
        electrons, at fixed external potential, evaluated at :math:`N_{\text{elec}}`.
        More specifically,

        .. math::
             \left. \left(\frac{\partial \Omega}{\partial \mu}
                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}} &= - N_{\text{elec}} \\
             \left. \left(\frac{\partial^2 \Omega}{\partial \mu^2}
                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}} &= -\frac{1}{\eta^{(1)}}

        where :math:`\eta^{(n)}` denotes the :math:`(n+1)^{\text{th}}`-order derivative of energy
        w.r.t. number of electrons evaluated at :math:`N_{\text{elec}}`, i.e.

        .. math::
           \eta^{(n)} = \left. \left(\frac{\partial^{n+1} E}{\partial N^{n+1}}
                               \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}}
                        \qquad n = 1, 2, \dots

        To compute higher-order derivatives, Faa di Bruno formula which generalizes the chain rule
        to higher derivatives can be used. i.e. for :math:`n \geq 2`,

        .. math::
           \left. \left(\frac{\partial^n \Omega}{\partial \mu^n}
                  \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}} =
           \frac{-\displaystyle\sum_{k=1}^{n-2} \left.\left(\frac{\partial^k \Omega}{\partial \mu^k}
                                    \right)_{v(\mathbf{r})} \right|_{N = N_{\text{elec}}}
                 \cdot B_{n-1,k} \left(\eta^{(1)}, \eta^{(2)}, \dots, \eta^{(n-k)} \right)}
                {B_{n-1,n-1} \left(\eta^{(1)}\right)}

        where :math:`B_{n-1,k} \left(x_1, x_2, \dots, x_{n-k}\right)` denotes the Bell polynomials.

        Parameters
        ----------
        n_elec : float
            Number of electrons, :math:`N_{\text{elec}}`.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.
        """
        if n_elec is not None and n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! #elec={0}'.format(n_elec))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')

        if n_elec is None:
            deriv = None
        elif order == 1:
            # 1st order derivative is minus number of electrons
            deriv = - n_elec
        elif order == 2:
            # 2nd order derivative is inverse hardness
            hardness = self.energy_derivative(n_elec, order=2)
            if hardness is not None and hardness != 0.0:
                deriv = -1.0 / hardness
            else:
                deriv = None
        else:
            # higher-order derivatives are compute with Faa Di Bruno formula
            # list of hyper-hardneses (derivatives of energy w.r.t. N)
            e_deriv = [self.energy_derivative(n_elec, i + 1) for i in range(1, order)]
            g_deriv = [self.grand_potential_derivative(n_elec, k + 1) for k in range(1, order - 1)]
            if any([item is None for item in e_deriv]) or any([item is None for item in g_deriv]):
                deriv = None
            else:
                deriv = 0
                for k in range(1, order - 1):
                    deriv -= g_deriv[k - 1] * sp.bell(order - 1, k, e_deriv[:order - k])
                deriv /= sp.bell(order - 1, order - 1, [e_deriv[0]])
        return deriv

    def grand_potential_mu(self, mu):
        r"""
        Evaluate the grand potential model for the specified chemical potential :math:`\mu`.

        To evaluate grand potential model, first the number of electrons corresponding to the
        specified :math:`\mu` is found, i.e. :math:`N(\mu)=\mu^{-1}(N)`, then the grand potential
        in computed by,

        .. math:: \Omega [\mu(N); v(\mathbf{r})] = E(N(\mu)) - \mu \times N(\mu)

        Parameters
        ----------
        mu : float
            Chemical potential :math:`\mu`.
        """
        # find N corresponding to the given mu
        n_elec = self.convert_mu_to_n(mu)
        # evaluate grand potential as a function of N
        value = self.grand_potential(n_elec)
        return value

    def grand_potential_mu_derivative(self, mu, order=1):
        r"""
        Evaluate the :math:`n^{\text{th}}`-order derivative of grand potential at the given mu.

        This returns the :math:`n^{\text{th}}`-order derivative of grand potential model w.r.t.
        chemical potential, at fixed external potential, evaluated for the specified chemical
        potential :math:`\mu`.

        That is,

        .. math::
             \left. \left(\frac{\partial^n \Omega}{\partial \mu^n}
                    \right)_{v(\mathbf{r})} \right|_{N = N\left(\mu\right)} =
             \left. \left(\frac{\partial^{n-1}}{\partial \mu^{n-1}}
                          \frac{\partial \Omega}{\partial \mu}
                    \right)_{v(\mathbf{r})} \right|_{N = N\left(\mu\right)} =
           - \left. \left(\frac{\partial^{n-1} N}{\partial \mu^{n-1}}
                    \right)_{v(\mathbf{r})} \right|_{N = N\left(\mu\right)}
             \quad  n = 1, 2, \dots

        To evaluate this expression, the number of electrons corresponding to the specified
        :math:`\mu` should is found, i.e. :math:`N(\mu)=\mu^{-1}(N)`.

        Parameters
        ----------
        mu : float
            Chemical potential, :math:`\mu`.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.
        """
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # find N corresponding to the given mu
        n_elec = self.convert_mu_to_n(mu)
        # evaluate grand potential derivative w.r.t. mu
        value = self.grand_potential_derivative(n_elec, order)
        return value

    def convert_mu_to_n(self, mu, guess=None):
        r"""
        Return the number of electrons :math:`N` matching the given chemical potential :math:`\mu`.

        Chemical potential is a function of the number of electrons, :math:`\mu(N)`, as it is the
        first derivative of energy model :math:`E(N)` with respect to the number of electrons at
        fixed external potential,

        .. math:: \mu(N) = \left(\frac{\partial E(N)}{\partial N}\right)_{v(\mathbf{r})}

        Here we solve for :math:`N` which results in the specified :math:`\mu` according to the
        equation above, i.e. :math:`N(\mu) = \mu^{-1}(N)`, using ``scipy.optimize.newton``.

        Parameters
        ----------
        mu : float
            Chemical potential, :math:`\mu`.
        guess : float, default=None
            Initial guess used for solving for :math:`N`.
            If ``None``, the reference number of electrons :math:`N_0` is used as an initial guess.
        """
        # assign an initial guess for N
        if guess is None:
            guess = self._n0
        # solve for N corresponding to the given mu using scipy.optimize.newton
        try:
            n_elec = newton(lambda n: self.energy_derivative(n, 1) - mu,
                            guess,
                            fprime=lambda n: self.energy_derivative(n, 2),
                            fprime2=lambda n: self.energy_derivative(n, 3))
        except ValueError:
            raise ValueError(
                'Number of electrons corresponding to mu={0} could not be found!'.format(mu))
        # according to scipy.optimize.newton notes: the stopping criterion used here is the step
        # size and there is no guarantee that a zero has been found. Consequently the result
        # should be verified.
        if not abs(self.energy_derivative(n_elec, 1) - mu) < 1.e-4:
            logging.warning("Solved number of electrons {0} corresponding to {1} gives, "
                            "mu(N={0})={2}".format(n_elec, mu, self.energy_derivative(n_elec, 1)))
        return n_elec


class BaseLocalTool(object):
    """Base class of local conceptual DFT reactivity descriptors."""

    def __init__(self, n0, n_max=None, global_softness=None):
        r"""
        Initialize class.

        Parameters
        ----------
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.
        """
        if n0 <= 0:
            raise ValueError("Argument n0 should be positive! Given n0={0}".format(n0))
        if n_max is not None and n_max < 0:
            raise ValueError("Argument n_max cannot be negative! Given n_max={0}".format(n_max))
        self._n0 = n0
        self._n_max = n_max
        self._global_softness = global_softness

    @property
    def n0(self):
        r"""Reference number of electrons, i.e. :math:`N_0`."""
        return self._n0

    @property
    def n_max(self):
        r"""Maximum number of electrons that the system accepts, i.e. :math:`N_{\text{max}}`."""
        return self._n_max

    @property
    def global_softness(self):
        r"""Global softness."""
        return self._global_softness

    def density(self, n_elec):
        r"""
        Evaluate density model :math:`\rho_N(\mathbf{r})` at the :math:`N_{\text{elec}}`.

        The functional derivative of energy model :math:`E(N)` w.r.t. external potential
        at fixed number of electrons, evaluated at the given number of electrons
        :math:`N_{\text{elec}}`, i.e.

        .. math::
           \left.\rho_N(\mathbf{r}) = {\left(\frac{\delta E(N)}{\delta v(\mathbf{r})}\right)_N}
           \right|_{N = N_{\text{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        """
        raise NotImplementedError

    def density_derivative(self, n_elec, order=1):
        r"""
        Evaluate n-th derivative of density w.r.t. number of electrons at :math:`N_{\text{elec}}`.

        The n-th order derivative of density model :math:`\rho_N(\mathbf{r})` w.r.t. the number
        of electrons, at fixed external potential, evaluated at the given number of electrons
        :math:`N_{\text{elec}}` is:

        .. math::
           \left. \left(\frac{\partial^n \rho_N(\mathbf{r})}{\partial N^n}
                  \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        order : int, optional
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{\text{elec}} = N_0` the first, second and higher order density
        derivatives correspond to the :attr:`fukui function <BaseLocalTool.fukui_function>`,
        :attr:`dual descriptor <BaseLocalTool.dual_descriptor>` and
        :attr:`hyper fukui function <BaseLocalTool.hyper_fukui_function>`, respectively.
        """
        raise NotImplementedError

    @property
    def fukui_function(self):
        r"""
        Fukui function of :math:`N_0`-electron system.

        This is defined as the 1st derivative of density model :math:`\rho_N(\mathbf{r})` w.r.t.
        the number of electrons, at fixed external potential, evaluated at :math:`N_0`, or the
        functional derivative of chemical potential w.r.t. external potential, at fixed number of
        electrons, i.e.

        .. math::
           f_{N_0}(\mathbf{r}) = {\left( \frac{\delta \mu}{\delta v(\mathbf{r})} \right)_N}
                   = \left. \left(\frac{\partial \rho_N(\mathbf{r})}{\partial N}
                            \right)_{v(\mathbf{r})}\right|_{N = N_0}

        where :math:`\mu` is the :attr:`chemical potential <BaseGlobalTool.chemical_potential>`.
        """
        return self.density_derivative(self._n0, order=1)

    @property
    def dual_descriptor(self):
        r"""
        Dual descriptor of :math:`N_0`-electron system.

        This is defined as the 2nd derivative of density model :math:`\rho_N(\mathbf{r})` w.r.t. the
        number of electrons, at fixed external potential, evaluated at :math:`N_0`, or the
        functional derivative of chemical hardness w.r.t. external potential, at fixed number of
        electrons, i.e.

        .. math::
           \Delta f_{N_0}(\mathbf{r}) = {\left( \frac{\delta \eta}{\delta v(\mathbf{r})} \right)_N}
                   = \left. \left(\frac{\partial^2 \rho_N(\mathbf{r})}{\partial N^2}
                            \right)_{v(\mathbf{r})}\right|_{N = N_0}

        where :math:`\eta` is the :attr:`chemical hardness <BaseGlobalTool.chemical_hardness>`.
        """
        return self.density_derivative(self._n0, order=2)

    # def hyper_fukui_function(self, order=3):
    #     return NotImplementedError

    @property
    def softness(self):
        r"""
        Chemical softness of :math:`N_0`-electron system.

        .. math::
           s_{N_0}\left(\mathbf{r}\right) = S \cdot f_{N_0}\left(\mathbf{r}\right)

        where :math:`S` is the :attr:`global softness <BaseGlobalTool.global_softness>`,
        and :math:`f_{N_0}\left(\mathbf{r}\right)` is
        :attr:`fukui function <BaseLocalTool.fukui_function>`.
        """
        if self._global_softness is None or self.fukui_function is None:
            return None
        return self._global_softness * self.fukui_function

    @property
    def hyper_softness(self):
        r"""
        Chemical hyper-softness of :math:`N_0`-electron system.

        .. math::
           s_{N_0}^{(2)}\left(\mathbf{r}\right) = S^2 \cdot \Delta f_{N_0}\left(\mathbf{r}\right)

        where :math:`S` is the :attr:`global softness <BaseGlobalTool.global_softness>`,
        and :math:`\Delta f_{N_0}\left(\mathbf{r}\right)` is the
        :attr:`dual descriptor <BaseLocalTool.dual_descriptor>`.
        """
        if self._global_softness is None or self.dual_descriptor is None:
            return None
        return self._global_softness**2 * self.dual_descriptor


class BaseCondensedTool(object):
    """Base class of condensed conceptual DFT reactivity descriptors."""

    def __init__(self, n0, n_max=None, global_softness=None):
        r"""
        Initialize class.

        Parameters
        ----------
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        n_max : float, optional
            Maximum number of electrons that system can accept, i.e. :math:`N_{\text{max}}`.
            See :attr:`BaseGlobalTool.n_max`.
        global_softness : float, optional
            Global softness. See :attr:`BaseGlobalTool.softness`.
        """
        if n0 <= 0:
            raise ValueError("Argument n0 should be positive! Given n0={0}".format(n0))
        if n_max is not None and n_max < 0:
            raise ValueError("Argument n_max cannot be negative! Given n_max={0}".format(n_max))
        self._n0 = n0
        self._n_max = n_max
        self._global_softness = global_softness

    @property
    def n_ref(self):
        r"""Reference number of electrons, i.e. :math:`N_0`."""
        return self._n0

    @property
    def n_max(self):
        r"""Maximum number of electrons that the system accepts, i.e. :math:`N_{\text{max}}`."""
        return self._n_max

    @property
    def global_softness(self):
        r"""Global softness."""
        return self._global_softness

    def population(self, n_elec):
        r"""
        Evaluate atomic populations at the given number of electrons :math:`N_{\text{elec}}`.

        ..math::
          N_A = Z_A - \int w_A(\mathbf{r}) \rho_N(\mathbf{r}) d\mathbf{r}

        where :math:`w_A(\mathbf{r})` is the atomic weight of atom :math:`A` at point \mathbf{r}.

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        """
        raise NotImplementedError

    def population_derivative(self, n_elec, order=1):
        r"""
        Evaluate n-th derivative of atomic populations w.r.t. number of electrons.

        The n-th order derivative of atomic populations :math:`\rho_N(\mathbf{r})` w.r.t.
        the number of electrons, at fixed chemical potential, evaluated at the given number
        of electrons :math:`N_{\text{elec}}` is:

        .. math::
           \left. \left(\frac{\partial^n \rho_N(\mathbf{r})}{\partial N^n}
                  \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons, :math:`N_{\text{elec}}`.
        order : int, optional
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{\text{elec}} = N_0` the first, second and higher order density
        derivatives correspond to the condensed :attr:`fukui function
        <BaseCondensedTool.fukui_function>`,
        :attr:`dual descriptor <BaseCondensedTool.dual_descriptor>` and
        :attr:`hyper fukui function <BaseCondensedTool.hyper_fukui_function>`, respectively.
        """
        raise NotImplementedError

    @property
    def fukui_function(self):
        r"""
        Atomic Fukui function of :math:`N_0`-electron system.

        This is defined as the 1st derivative of density model :math:`\rho_N(\mathbf{r})` w.r.t.
        the number of electrons, at fixed external potential, evaluated at :math:`N_0`, or the
        functional derivative of chemical potential w.r.t. external potential, at fixed number of
        electrons, i.e.

        .. math::
           f_{N_0}(\mathbf{r}) = {\left( \frac{\delta \mu}{\delta v(\mathbf{r})} \right)_N}
                   = \left. \left(\frac{\partial \rho_N(\mathbf{r})}{\partial N}
                            \right)_{v(\mathbf{r})}\right|_{N = N_0}

        where :math:`\mu` is the :attr:`chemical potential <BaseGlobalTool.chemical_potential>`.
        """
        return self.population_derivative(self._n0, order=1)

    @property
    def dual_descriptor(self):
        r"""
        Atomic dual descriptor of :math:`N_0`-electron system.

        This is defined as the 2nd derivative of density model :math:`\rho_N(\mathbf{r})` w.r.t. the
        number of electrons, at fixed external potential, evaluated at :math:`N_0`, or the
        functional derivative of chemical hardness w.r.t. external potential, at fixed number of
        electrons, i.e.

        .. math::
           \Delta f_{N_0}(\mathbf{r}) = {\left( \frac{\delta \eta}{\delta v(\mathbf{r})} \right)_N}
                   = \left. \left(\frac{\partial^2 \rho_N(\mathbf{r})}{\partial N^2}
                            \right)_{v(\mathbf{r})}\right|_{N = N_0}

        where :math:`\eta` is the :attr:`chemical hardness <BaseGlobalTool.chemical_hardness>`.
        """
        return self.population_derivative(self._n0, order=2)

    # def hyper_fukui_function(self, order=3):
    #     return NotImplementedError

    @property
    def softness(self):
        r"""
        Atomic chemical softness of :math:`N_0`-electron system.

        .. math::
           s_A\left(N\right) = S \cdot f_A\left(\mathbf{r}\right)

        where :math:`S` is the :attr:`global softness <BaseGlobalTool.global_softness>`.
        """
        if self._global_softness is None or self.fukui_function is None:
            return None
        return self._global_softness * self.fukui_function

    @property
    def hyper_softness(self):
        r"""
        Chemical hyper-softness of :math:`N_0`-electron system, :math:`s_N^{(2)}(\mathbf{r})`.

        .. math::
           s_N^{(2)}\left(\mathbf{r}\right) = S^2 \cdot \Delta f_N\left(\mathbf{r}\right)

        where :math:`S` is the :attr:`global softness <BaseGlobalTool.global_softness>` and
        :math:`\Delta f_N` is the :attr:`dual descriptor <BaseLocalTool.dual_descriptor>`.
        """
        if self._global_softness is None or self.dual_descriptor is None:
            return None
        return self._global_softness**2 * self.dual_descriptor
