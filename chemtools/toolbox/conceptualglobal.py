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
"""Global Conceptual Density Functional Theory (DFT) Reactivity Tools.

   This module contains various global tool classes corresponding to
   linear, quadratic, exponential, general energy models.
"""

import math
import warnings
from abc import ABCMeta, abstractmethod
import numpy as np
import sympy as sp
from scipy.optimize import root, newton, least_squares
from chemtools.utils.utils import doc_inherit

__all__ = ['BaseGlobalTool', 'LinearGlobalTool', 'QuadraticGlobalTool',
           'ExponentialGlobalTool', 'RationalGlobalTool', 'GeneralGlobalTool']


class BaseGlobalTool(object):
    """
    Base class of global conceptual DFT reactivity descriptors.
    """
    __metaclass__ = ABCMeta

    def __init__(self, energy_zero, energy_plus, energy_minus, n0, n_max):
        """
        Parameters
        ----------
        energy_zero : float
            Energy of the :math:`N_0` -electron system, i.e. :math:`E(N_0)`.
        energy_plus : float
            Energy of the :math:`(N_0 + 1)` -electron system, i.e. :math:`E(N_0 + 1)`.
        energy_minus : float
            Energy of the :math:`(N_0 - 1)` -electron system, i.e. :math:`E(N_0 - 1)`.
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        """
        if n0 <= 0:
            raise ValueError('Argument n0 should be positive! Given n0={0}'.format(n0))
        if n_max is not None and n_max < 0:
            raise ValueError('Argument n_max cannot be negative! Given n0={0}'.format(n_max))
        self._n0 = n0
        self._n_max = n_max
        self._energy_zero = energy_zero
        self._energy_plus = energy_plus
        self._energy_minus = energy_minus
        # calculate ionization potential and electron affinity
        self._ip = energy_minus - energy_zero
        self._ea = energy_zero - energy_plus

    @property
    def n0(self):
        """
        Reference number of electrons, i.e. :math:`N_0`.
        """
        return self._n0

    @property
    def n_max(self):
        r"""
        Maximum number of electrons that the system can accept, defined as,

        .. math:: N_{\text{max}} = \underbrace {\min }_N E(N)
        """
        return self._n_max

    @property
    def energy_zero(self):
        """
        Energy of the system with :math:`N_0` electrons, i.e. :math:`E(N_0)`.
        """
        return self._energy_zero

    @property
    def ionization_potential(self):
        """
        Ionization potential (IP) of the :math:`N_0` -electron system defined as,

        .. math:: IP = E(N_0 - 1) - E(N_0)
        """
        return self._ip

    @property
    def ip(self):
        """
        The same as :attr:`ionization_potential`.
        """
        return self.ionization_potential

    @property
    def electron_affinity(self):
        """
        Electron affinity (EA) of the :math:`N_0` -electron system defined as,

        .. math:: EA = E(N_0) - E(N_0 + 1)
        """
        return self._ea

    @property
    def ea(self):
        """
        The same as :attr:`electron_affinity`.
        """
        return self.electron_affinity

    @property
    def electronegativity(self):
        r"""
        Mulliken electronegativity defined as negative :attr:`chemical_potential`,

        .. math:: \chi_{\text{Mulliken}} = - \mu
        """
        if self.chemical_potential is None:
            return None
        value = -1 * self.chemical_potential
        return value

    @property
    def electrophilicity(self):
        r"""
        Electrophilicity defined as,

        .. math::
           \omega_{\text{electrophilicity}} = \text{sgn}\left(N_{\text{max}} - N_0\right)
                                              \times \left(E(N_0) - E(N_{\text{max}})\right)
        """
        if self._n_max is None:
            return None
        sign = np.sign(self._n_max - self._n0)
        value = sign * (self._energy_zero - self.energy(self._n_max))
        return value

    @property
    def nucleofugality(self):
        r"""
        Nucleofugality defined as,

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
        Electrofugality defined as,

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
        Chemical potential defined as the first derivative of the energy model w.r.t.
        the number of electrons at fixed external potential evaluated at :math:`N_0`,

        .. math::
           \mu = \left. \left(\frac{\partial E}{\partial N}
                        \right)_{v(\mathbf{r})} \right|_{N = N_0}
        """
        value = self.energy_derivative(self._n0, order=1)
        return value

    @property
    def mu(self):
        """
        The same as :attr:`chemical_potential`.
        """
        return self.chemical_potential

    @property
    def chemical_hardness(self):
        r"""
        Chemical hardness defined as the second derivative of the energy model w.r.t.
        the number of electrons at fixed external potential evaluated at :math:`N_0`,

        .. math::
           \eta  = \left. \left(\frac{\partial^2 E}{\partial N^2}
                          \right)_{v(\mathbf{r})} \right|_{N = N_0}
        """
        value = self.energy_derivative(self._n0, order=2)
        return value

    @property
    def eta(self):
        """
        The same as :attr:`chemical_hardness`.
        """
        return self.chemical_hardness

    @property
    def softness(self):
        r"""
        Chemical softness defined as the second derivative of the grand potential model w.r.t
        the number of electrons at fixed external potential evaluated at :math:`N_0`. This is
        equal to the inverse chemical hardness.

        .. math::
           S = - \left. \left(\frac{\partial^2 \Omega}{\partial \mu^2}
                        \right)_{v(\mathbf{r})}\right|_{N = N_0} = \frac{1}{\eta}
        """
        # compute 2nd-order derivative of grand potential w.r.t. mu at N0
        value = - self.grand_potential_derivative(self._n0, 2)
        return value

    def hyper_hardness(self, order=2):
        r"""
        :math:`n^{\text{th}}`-order hyper hardness defined as the :math:`(n+1)^{\text{th}}`-order
        derivative, where :math:`n \geq 2`, of the energy model w.r.t the number of electrons at
        fixed external potential evaluated at :math:`N_0`.

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
        :math:`n^{\text{th}}`-order hyper softness defined as the :math:`(n+1)^{\text{th}}`-order
        derivative, where :math:`n \geq 2`, of the grand potential model w.r.t the number of
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
        value = - self.grand_potential_derivative(self._n0, order + 1)
        return value

    @abstractmethod
    def energy(self, n_elec):
        r"""
        Return energy model :math:`E(N)` evaluated for the specified number of electrons,
        i.e. :math:`E(N_{\text{elec}})`.

        Parameters
        ----------
        n_elec: float
            Number of electrons.
        """
        pass

    @abstractmethod
    def energy_derivative(self, n_elec, order=1):
        r"""
        Return the :math:`n^{\text{th}}`-order derivative of energy model :math:`E(N)` w.r.t.
        to the number of electrons at fixed chemical potential evaluated for the specified number
        of electrons.

        .. math::
           \left. \left(\frac{\partial^n E}{\partial N^n}
                  \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{\text{elec}} = N_0` the first, second and higher order derivatives are equal
        to the :attr:`BaseGlobalTool.chemical_potential`, :attr:`BaseGlobalTool.chemical_hardness`
        and :attr:`BaseGlobalTool.hyper_hardness`, respectively.
        """
        pass

    def grand_potential(self, n_elec):
        r"""
        Return the grand potential model evaluated for the specified number of electrons
        :math:`N_{\text{elec}}`.

        That is,

        .. math::
           \Omega [\mu(N_{\text{elec}}); v(\mathbf{r})] &=
            E(N_{\text{elec}}) - \mu(N_{\text{elec}}) \times N_{\text{elec}} \\
            &= E(N_{\text{elec}}) - \left.\left(\frac{\partial E(N)}{\partial N}
                    \right)_{v(\mathbf{r})}\right|_{N = N_{\text{elec}}} \times N_{\text{elec}}

        Parameters
        ----------
        n_elec : float
            Number of electrons.
        """
        if n_elec is None:
            return None
        # compute grand potential as a function of N
        value = self.energy(n_elec) - self.energy_derivative(n_elec, 1) * n_elec
        return value

    def grand_potential_derivative(self, n_elec, order=1):
        r"""
        Return the :math:`n^{\text{th}}`-order derivative of grand potential model w.r.t.
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
            Number of electrons.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.
        """
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! #elec={0}'.format(n_elec))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')

        if order == 1:
            # 1st order derivative is minus number of electrons
            deriv = - n_elec
        elif order == 2:
            # 2nd order derivative is inverse hardness
            hardness = self.energy_derivative(n_elec, order=2)
            if hardness is None:
                deriv = None
            elif hardness != 0.0:
                deriv = - 1.0 / hardness
            else:
                deriv = None
        else:
            # higher-order derivatives are compute with Faa Di Bruno formula
            # list of hyper-hardneses (derivatives of energy w.r.t. N)
            energy_deriv = [self.energy_derivative(n_elec, i + 1) for i in range(1, order)]
            deriv = 0
            for k in xrange(1, order - 1):
                grand_potential = self.grand_potential_derivative(n_elec, k + 1)
                deriv -= grand_potential * sp.bell(order - 1, k, energy_deriv[:order - k])
            deriv /= sp.bell(order - 1, order - 1, [energy_deriv[0]])
        return deriv

    def grand_potential_mu(self, mu):
        r"""
        Return the grand potential model evaluated for the specified chemical potential :math:`\mu`.

        To evaluate grand potential model,first the number of electrons corresponding to the
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
        Return the :math:`n^{\text{th}}`-order derivative of grand potential model w.r.t.
        to the chemical potential, at fixed external potential, evaluated for the specified
        chemical potential :math:`\mu`.

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
            Chemical potential.
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
        Return the number of electrons, :math:`N`, corresponding to the specified chemical potential
        :math:`\mu`.

        Chemical potential is a function of the number of electrons, :math:`\mu(N)`, as it is the
        first derivative of energy model :math:`E(N)` with respect to the number of electrons at
        fixed external potential,

        .. math:: \mu(N) = \left(\frac{\partial E(N)}{\partial N}\right)_{v(\mathbf{r})}

        Here we solve for :math:`N` which results in the specified :math:`\mu` according to the
        equation above, i.e. :math:`N(\mu) = \mu^{-1}(N)`, using ``scipy.optimize.newton``.

        Parameters
        ----------
        mu : float
            Chemical potential.
        guess : float, default=None
            Initial guess used for solving for :math:`N`.
            If ``None``, the reference number of electrons :math:`N_0` is used as an initial guess.
        """
        # assign an initial guess for N
        if guess is None:
            guess = self._n0
        # solve for N corresponding to the given mu using scipy.optimize.newton
        func = lambda n: self.energy_derivative(n, 1) - mu
        fprime = lambda n: self.energy_derivative(n, 2)
        fprime2 = lambda n: self.energy_derivative(n, 3)
        try:
            n_elec = newton(func, guess, fprime=fprime, fprime2=fprime2)
        except:
            raise ValueError(
                'Number of electrons corresponding to mu={0} could not be found!'.format(mu))
        # according to scipy.optimize.newton notes: the stopping criterion used here is the step
        # size and there is no guarantee that a zero has been found. Consequently the result
        # should be verified.
        if not abs(self.energy_derivative(n_elec, 1) - mu) < 1.e-4:
            warnings.warn('Solved number of electrons {0} corresponding to {1}'.format(n_elec, mu) +
                          ' gives, mu(N={0})={1}'.format(n_elec, self.energy_derivative(n_elec, 1)),
                          RuntimeWarning)
        return n_elec


class LinearGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the linear energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values.

    The energy is approximated as a piece-wise linear function of the number of electrons:

    .. math::
       E(N) = \begin{cases}
               (N_0 - N) \times E(N_0 - 1) + (N - N_0 + 1) \times E(N_0) & N \leqslant N_0 \\
               (N_0 + 1 - N) \times E(N_0) + (N - N_0) \times E(N_0 + 1) & N \geqslant N_0 \\
              \end{cases}

    Because the energy model is not differentiable at integer number of electrons, the first
    derivative of the energy w.r.t. number if electrons is calculated from above, from below
    and averaged:

    .. math::
       \mu^{-} &= -IP \\\
       \mu^{0} &= \frac{\mu^{+} + \mu^{-}}{2} \\
       \mu^{+} &= -EA \\
    """
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        if energy_zero < energy_plus:
            n_max = n0
        else:
            n_max = None
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def mu_minus(self):
        r"""
        Chemical potential from below, i.e. :math:`N_0^{-}`, given by,

        .. math:: \mu^{-} = E\left(N_0\right) - E\left(N_0 - 1\right)  = -IP
        """
        return -1 * self._ip

    @property
    def mu_plus(self):
        r"""
        Chemical potential from above, i.e. :math:`N_0^{+}`, given by,

        .. math:: \mu^{+} = E\left(N_0 + 1\right) - E\left(N_0\right)  = -EA
        """
        return -1 * self._ea

    @property
    def mu_zero(self):
        r"""
        Chemical potential averaged, given by,

        .. math::
           \mu^{0} = \frac{\mu^{+} + \mu^{-}}{2}
                   = \frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2} = - \frac{IP + EA}{2}
        """
        return -0.5 * (self._ea + self._ip)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec is None:
            return None
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        # evaluate energy
        value = self._energy_zero
        if n_elec < self._n0:
            value += (self._n0 - n_elec) * self._ip
        elif n_elec > self._n0:
            value += (self._n0 - n_elec) * self._ea
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec is None:
            return None
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy derivative evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if n_elec == self._n0:
            deriv = None
        elif order >= 2:
            deriv = 0.0
        elif n_elec < self._n0:
            deriv = - self._ip
        elif n_elec > self._n0:
            deriv = - self._ea
        return deriv


class QuadraticGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the quadratic energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy.

    The energy is approximated as a quadratic function of the number of electrons,
    and the three unknown parameters are obtained by fitting to the given values of energy.

    .. math:: E(N) = a + b N + c {N^2}

    First, second and higher order derivatives of the quadratic energy model with respect to
    the number of electrons at fixed external potential are given by:

    .. math::
       \left(\frac{\partial E}{\partial N}\right)_{v(\mathbf{r})} &= b + 2 c N \\
       \left(\frac{\partial^2 E}{\partial N^2}\right)_{v(\mathbf{r})} &= 2 c \\
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} &= 0
             \quad \text{for} \quad n \geq 2
    """
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # calculate parameters a, b, c of quadratic energy model
        c = 0.5 * (energy_minus - 2 * energy_zero + energy_plus)
        b = 0.5 * (energy_plus - energy_minus) - 2 * c * n0
        a = energy_zero - b * n0 - c * (n0**2)
        self._params = [a, b, c]
        # calculate Nmax (number of electrons for which energy is minimum)
        n_max = - b / (2 * c)
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        """
        Parameters :math:`a`, :math:`b` and :math:`c` of energy model.
        """
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        # evaluate energy
        value = self._params[0] + self._params[1] * n_elec + self._params[2] * n_elec**2
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy derivative evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if order == 1:
            deriv = self._params[1] + 2 * n_elec * self._params[2]
        elif order == 2:
            deriv = 2 * self._params[2]
        elif order >= 2:
            deriv = 0.0
        return deriv


class ExponentialGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the exponential energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy.

    The energy is approximated as a exponential function of the number of electrons,
    and the three unknown parameters are obtained by interpolating to the given values of energy.

    .. math:: E(N) = A \exp(-\gamma(N-N_0)) + B

    The :math:`n^{\text{th}}`-order derivative of the rational energy model with respect to
    the number of electrons at fixed external potential is given by:

    .. math::
       \left(\frac{\partial^n E}{\partial N^n}\right)_{v(\mathbf{r})} =
              A (-\gamma)^n \exp(-\gamma (N - N_0))
    """
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # check energy values are monotonic, i.e. E(N-1) > E(N) > E(N+1)
        if not (energy_minus > energy_zero > energy_plus):
            energies = [energy_minus, energy_zero, energy_plus]
            n_values = [n0 - 1, n0, n0 + 1]
            raise ValueError('To interpolate exponential energy model, E values vs. N should be ' +
                             'monotonic! Given E={0} for N={1}.'.format(energies, n_values))

        # calculate the A, B, gamma parameters of the model and N_max
        param_a = (energy_minus - energy_zero) * (energy_zero - energy_plus)
        param_a /= (energy_minus - 2 * energy_zero + energy_plus)
        param_b = energy_zero - param_a
        gamma = (energy_minus - 2 * energy_zero + energy_plus) / (energy_plus - energy_zero)
        gamma = math.log(1. - gamma)
        self._params = [param_a, gamma, param_b]
        # calulate N_max
        n_max = float('inf')
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        r"""
        Parameters :math:`A`, :math:`\gamma` and :math:`B` of energy model.
        """
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        # evaluate energy
        if np.isinf(n_elec):
            # limit of E(N) as N goes to infinity equals B
            value = self._params[2]
        else:
            dn = n_elec - self._n0
            value = self._params[0] * math.exp(- self._params[1] * dn) + self._params[2]
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy derivative evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if np.isinf(n_elec):
            # limit of E(N) derivatives as N goes to infinity equals zero
            deriv = 0.0
        else:
            dn = n_elec - self._n0
            deriv = self._params[0] * (- self._params[1])**order * math.exp(- self._params[1] * dn)
        return deriv


class RationalGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the 3-parameter
    rational energy model, given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)`
    known values of energy.

    The energy is approximated as a rational function of the number of electrons,
    and the three unknown parameters are obtained by interpolating to the given values of energy.

    .. math:: E(N) = \frac{a_0 + a_1 N}{1 + b_1 N}

    The :math:`n^{\text{th}}`-order derivatives of the rational energy model with respect to
    the number of electrons at fixed external potential is given by:

    .. math::
       \left(\frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})} =
             \frac{b_1^{n - 1} (a_1 - a_0 b_1) n!}{(1 + b_1 N)^{2n}}
    """
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # check energy values are monotonic, i.e. E(N-1) > E(N) > E(N+1)
        if not energy_minus > energy_zero > energy_plus:
            energies = [energy_minus, energy_zero, energy_plus]
            n_values = [n0 - 1, n0, n0 + 1]
            raise ValueError('To interpolate rational energy model, E values vs. N should be ' +
                             'monotonic! Given E={0} for N={1}.'.format(energies, n_values))

        # calculate parameters a0, a1 and b1 of rational energy model
        b1 = - (energy_plus - 2 * energy_zero + energy_minus)
        b1 /= ((n0 + 1) * energy_plus - 2 * n0 * energy_zero + (n0 - 1) * energy_minus)
        a1 = (1 + b1 * n0) * (energy_plus - energy_zero) + (b1 * energy_plus)
        a0 = - a1 * n0 + energy_zero * (1 + b1 * n0)
        self._params = [a0, a1, b1]
        # calculate Nmax
        n_max = float('inf')
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        """
        Parameters :math:`a_0`, :math:`a_1` and :math:`b_1` of energy model.
        """
        return self._params

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        # evaluate energy
        if np.isinf(n_elec):
            # limit of E(N) as N goes to infinity equals a1/b1
            value = self._params[1] / self._params[2]
        else:
            value = (self._params[0] + self._params[1] * n_elec) / (1 + self._params[2] * n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy derivative evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # evaluate derivative
        if np.isinf(n_elec):
            # limit of E(N) derivatives as N goes to infinity equals zero
            deriv = 0.0
        else:
            deriv = (-self._params[2])**(order - 1)
            deriv *= (self._params[1] - self._params[0] * self._params[2]) * math.factorial(order)
            deriv /= (1 + self._params[2] * n_elec)**(order + 1)
        return deriv


class GeneralGlobalTool(BaseGlobalTool):
    """
    Class of global conceptual DFT reactivity descriptors based on the user-specified
    symbolic energy model and given known values of energy.

    The energy is approximated as a symbolic function of the number of electrons,
    and the unknown parameters are obtained by interpolating to the given values of energy; i.e.
    the number of parameters in the model should equal the number of given energy values and
    the corresponding number of electrons.

    The :math:`n^{\text{th}}`-order derivative of the symbolic energy model with respect to the
    number of electrons at fixed external potential is calculated symbolically.
    """
    def __init__(self, expr, n0, n_energies, n_symbol=None, n0_symbol=None, guess=None, opts=None):
        """
        Parameters
        ----------
        expr : sp.Exp
            The energy expression representing the dependence of energy on the number of electrons.
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        n_energies : dict
            The energy values of `expr` at different electron-numbers.  The dict has int
            (electron-number) keys, float (energy) values.
        n_symbol : sp.Symbol, default=sp.symbols('N')
            The symbol in `expr` that represents the number of electrons.
        n0_symbol: sp.Symbol, optional
            The symbol in `expr` that represents the electron-number at which to evaluate
            `expr`.  If not specified, assume that it is already expressed numerically in
            `expr`.
        guess : dict, optional
            Guesses at the values of the parameters of `expr`.  The dict has sp.Symbol
            keys, float values.
        opts : dict, optional
            Optional keyword arguments to pass to the :py:meth:`scipy.optimize.root` solver
            that is used to solve for the parameters in the model.
        """
        # make sure that the energy expression depends on number of electrons
        if n_symbol is None:
            n_symbol = sp.symbols('N')
        if n_symbol not in expr.atoms(sp.Symbol):
            raise ValueError(
                'The expr={0} does not contain {1} symbol representing '.format(expr, n_symbol) +
                'the number of electrons.')
        self._n_symb = n_symbol
        # store minimum and maximum number of electrons used for interpolation
        self._n_min, self._n_max = np.min(n_energies.keys()), np.max(n_energies.keys())

        # substitute N0 in energy expression
        if n0_symbol:
            expr = expr.subs(n0_symbol, n0)

        # list of energy model parameters
        params = expr.atoms(sp.Symbol)
        params.remove(self._n_symb)
        # assign initial values for parameters of energy model
        if guess is None:
            guess = {}
        guess.update({param: 1. for param in params if param not in guess})
        # solve for the parameters of energy model
        self._params = self._solve_parameters(expr, n_energies, guess, opts)

        # substitute values of parameters in energy expression
        self._expr = expr.subs(self._params.items())

        # solve for N_max (number of electrons for which the 1st derivative of energy is zero)
        n_max = self._solve_nmax(n0)

        # calculate E(N0 - 1), E(N0) and E(N0 + 1) values
        energy_zero = self.energy(n0)
        energy_plus = self.energy(n0 + 1)
        energy_minus = self.energy(n0 - 1)
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0, n_max)

    @property
    def params(self):
        """
        Parameters dictionary of energy model.
        """
        return self._params

    @property
    def n_symbol(self):
        """
        Symbol used to denote the number of electrons.
        """
        return self._n_symb

    @property
    def expression(self):
        """
        Energy expression as a function of number of electrons, :math:`E(N)`.
        """
        return self._expr

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n_min <= n_elec <= self._n_max:
            warnings.warn('Energy evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n_min, self._n_max))
        # evaluate energy
        value = self._expr.subs(self._n_symb, n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            warnings.warn('Energy derivative evaluated for n_elec={0} outside of '.format(n_elec) +
                          'interpolation region [{0}, {1}].'.format(self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')

        # obtain derivative expression
        deriv = self._expr.diff(self._n_symb, order)
        # evaluate derivative expression at n_elec
        deriv = deriv.subs(self._n_symb, n_elec)
        return deriv

    def _solve_parameters(self, expr, n_energies, guess, opts=None):
        r"""
        Solve for the unknown parameters of the energy model using the given
        :math:`\left(N, E(N)\right)` pair. This requires solving a system of
        (non-)linear equations

        Parameters
        ----------
        See __init__().

        Returns
        -------
        parameters : dict
            A dictionary of sympy.Symbol keys corresponding to the
            value of the expression's solved parameters.
        """
        # obtain set of parameters in the energy expression
        params = guess.keys()
        if len(params) == 0:
            raise ValueError(
                'There is no parameters in the energy_expression={0} to solve for.'.format(expr))
        if not all([param in expr.atoms(sp.Symbol) for param in params]):
            raise ValueError('The expr={0} does not contain parameters given in guess={1}'.format(expr, guess))
        if len(params) > len(n_energies):
            raise ValueError('Underdetermined system of equations: Number of unknowns parameters '
                             'in the energy model is more than number of given known energies.')

        # initial guess for the parameters in the energy model
        guess = np.array([guess[param] for param in params])

        # construct system of equations to solve
        system_eqns = []
        d_system_eqns = []
        for n, energy in n_energies.iteritems():
            eqn = sp.lambdify((params,), expr.subs(self._n_symb, n) - energy, 'numpy')
            system_eqns.append(eqn)
            d_eqn_row = []
            for p in params:
                d_eqn = sp.lambdify((params,), expr.diff(p).subs(self._n_symb, n), 'numpy')
                d_eqn_row.append(d_eqn)
            d_system_eqns.append(d_eqn_row)

        def objective(args):
            """
            Evaluate the system of equations for the given values of parameters.

            Parameters
            ----------
            args : array representing the value of parameters.
                The expression for the property.
            """
            return np.array([eqn(args) for eqn in system_eqns])

        def jacobian(args):
            """
            Evaluate the Jacobian of the above objective function for the given
            values of parameters.

            Parameters
            ----------
            See objective().
            """
            jac = []
            for row in d_system_eqns:
                jac.append([eqn(args) for eqn in row])
            return np.array(jac)

        # solve for the parameters in the energy model
        if opts is None:
            opts = {}
        result = least_squares(objective, guess, jac=jacobian, **opts)
        if not result.success:
            raise ValueError(
                'The system of equations for parameters could not be solved. '
                'message:{0}'.format(result.message))
        # make dictionary of parameter values
        parameters = dict([(param, result.x[i]) for i, param in enumerate(params)])
        return parameters

    def _solve_nmax(self, guess):
        """
        Solve for the :math:`N_{\text{max}}` of the energy model.
        """
        d_expr = self._expr.diff(self._n_symb)
        n_max_eqn = sp.lambdify(self._n_symb, d_expr, 'numpy')
        result = root(n_max_eqn, guess)
        print result
        if result.success:
            n_max = np.asscalar(result.x)
            # n_ceil = math.ceil(n_max)
            # n_floor = math.floor(n_max)
            # e_ceil = self._expr.subs(self._n_symb, math.ceil(n_max))
            # e_floor = self._expr.subs(self._n_symb, math.floor(n_max))
            # n_max = n_floor if e_floor < e_ceil else n_ceil
        else:
            for sign in (+1, -1):
                n_inf = n_max_eqn(sign * np.inf)
                if np.isfinite(n_inf):
                    n_max = sign * np.inf
                    break
            else:
                n_max = None
                warnings.warn('The system of equations for Nmax could not be solved; Nmax=`None`.'
                              ' message:{0}'.format(result.message))
        return n_max
