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
'''Global Conceptual Density Functional Theory (DFT) Reactivity Tools.

   This module contains various global tool classes corresponding to
   linear, quadratic, exponential, general energy models.
'''

import math
import numpy as np
import sympy as sp
from scipy.optimize import root
from chemtools.utils import doc_inherit
from abc import ABCMeta, abstractmethod


class BaseGlobalTool(object):
    __metaclass__ = ABCMeta
    '''
    Base class of global conceptual DFT reactivity descriptors.
    '''
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        '''
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
        '''
        if n0 <= 0:
            raise ValueError('Argument n0 should be positive.')
        self._n0 = n0
        self._energy_zero = energy_zero
        self._energy_plus = energy_plus
        self._energy_minus = energy_minus
        # Calculate ionization potential and electron affinity
        self._ip = energy_minus - energy_zero
        self._ea = energy_zero - energy_plus

    @property
    def n0(self):
        '''
        Reference number of electrons, i.e. :math:`N_0`.
        '''
        return self._n0

    @property
    def n_max(self):
        r'''
        Maximum number of electrons that the system can accept, defined as:

        .. math:: N_{max} = \text{argmin} E(N)
        '''
        return self._n_max

    @property
    def energy_zero(self):
        '''
        Energy of the system with :math:`N_0` electrons, i.e. :math:`E(N_0)`.
        '''
        return self._energy_zero

    @property
    def ionization_potential(self):
        '''
        Ionization potential (IP) of the :math:`N_0` -electron system defined as:

        .. math:: I = E(N_0 - 1) - E(N_0)
        '''
        return self._ip

    @property
    def ip(self):
        '''
        The same as :attr:`ionization_potential`.
        '''
        return self.ionization_potential

    @property
    def electron_affinity(self):
        '''
        Electron affinity (EA) of the :math:`N_0` -electron system defined as:

        .. math:: A = E(N_0) - E(N_0 + 1)
        '''
        return self._ea

    @property
    def ea(self):
        '''
        The same as :attr:`electron_affinity`.
        '''
        return self.electron_affinity

    @property
    def electronegativity(self):
        r'''
        Mulliken electronegativity defined as negative :attr:`chemical_potential`.

        .. math:: \chi_{Mulliken} = - \mu
        '''
        value = -1 * self.chemical_potential
        return value

    @property
    def electrophilicity(self):
        r'''
        Electrophilicity defined as:

        ..math:: \omega_{\text {electrophile}} &= \text {sgn}(N_0 - N_{max}) (E(N_0) - E(N_{max}))
        '''
        sign = math.copysign(1, self._n0 - self._n_max)
        value = sign * (self._energy_zero - self.energy(self._n_max))
        #value = math.pow(self.mu, 2) / (2 * self.eta)
        return value

    @property
    def nucleofugality(self):
        r'''
        Nucleofugality defined as:

        .. math:: \nu_{\text {nucleofuge}} &= \text {sgn}(N_0 + 1 - N_{max}) (E(N_0 + 1) - E(N_{max}))
        '''
        sign = math.copysign(1, self._n0 + 1 - self._n_max)
        value = sign * (self.energy(self._n0 + 1) - self.energy(self._n_max))
        # value = math.pow(self._ip - 3 * self._ea, 2)
        # value /= (8 * (self._ip - self._ea))
        return value

    @property
    def electrofugality(self):
        r'''
        Electrofugality defined as:

        .. math:: \nu_{\text {electrofuge}} &= \text {sgn}(N_0 - 1 - N_{max}) (E(N_0 - 1) - E(N_{max}))
        '''
        sign = math.copysign(1, self._n0 - 1 - self._n_max)
        value = sign * (self.energy(self._n0 - 1) - self.energy(self._n_max))
        # value = math.pow(3 * self._ip - self._ea, 2)
        # value /= (8 * (self._ip - self._ea))
        return value

    @property
    def chemical_potential(self):
        r'''
        Chemical potential defined as the first derivative of the energy model w.r.t.
        the number of electrons at fixed external potential evaluated at :math:`N_0`:

        .. math::

           \mu  = {\left. {\left( {\frac{{\partial E}}{{\partial N}}} \right)_{v(r)}} \right|_{N = N_0}}
        '''
        value = self.energy_derivative(self._n0, order=1)
        return value

    @property
    def mu(self):
        '''
        The same as :attr:`chemical_potential`.
        '''
        return self.chemical_potential

    @property
    def chemical_hardness(self):
        r'''
        Chemical hardness defined as the second derivative of the energy model w.r.t.
        the number of electrons at fixed external potential evaluated at :math:`N_0`:

        .. math::

           \eta  = {\left. {\left( {\frac{{{\partial ^2}E}}{{\partial {N^2}}}} \right)_{v(r)}} \right|_{N = N_0}}
        '''
        value = self.energy_derivative(self._n0, order=2)
        return value

    @property
    def eta(self):
        '''
        The same as :attr:`chemical_hardness`.
        '''
        return self.chemical_hardness

    @property
    def softness(self):
        r'''
        Chemical softness defined as the second derivative of the grand potential model w.r.t
        the number of electrons at fixed external potential evaluated at :math:`N_0`. This is
        equal to the inverse chemical hardness.

        .. math::

           S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(r)} \right|_{N = N_0} = \frac{1}{\eta}
        '''
        hardness = self.chemical_hardness
        if hardness != 0.0:
            value = 1.0 / self.chemical_hardness
        else:
            value = None
        #TODO: when grand potential derivative is ready, it is better to calculate softness through that.
        #value = self.grand_potential_derivative(self._n0, 1)
        return value

    def hyper_hardness(self, order=2):
        r'''
        :math:`n^{th}`-order Hyper hardness defined as the :math:`(n+1)^{th}` -order derivative, where :math:`n \geq 2`,
        of the energy model w.r.t the number of electrons at fixed external potential evaluated
        at :math:`N_0`.

        .. math::

           \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(r)} \right|_{N = N_0} \text{for } n \geq 2

        Parameters
        ----------
        order : int, default=2
            The order of hyper-hardness denoted by :math:`n \geq 2` in the formula.
        '''
        if not (isinstance(order, int) and order >= 2):
            raise ValueError('Argument order should be an integer greater than or equal to 2.')
        value = self.energy_derivative(self._n0, order + 1)
        return value

    def hyper_softness(self, order):
        r'''
        :math:`n^{th}`-order Hyper softness defined as the :math:`(n+1)^{th}` -order derivative, where :math:`n \geq 2`,
        of the grand potential model w.r.t the number of electrons at fixed external potential
        evaluated at :math:`N_0`.

        .. math::

           S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)} \right|_{N = N_0} \text{for } n \geq 2

        Parameters
        ----------
        order : int, default=2
            The order of hyper-hardness denoted by :math:`n \geq 2` in the formula.
        '''
        if not (isinstance(order, int) and order >= 2):
            raise ValueError('Argument order should be an integer greater than or equal to 2.')
        value = self.grand_potential_derivative(self._n0, order + 1)
        return value

    @abstractmethod
    def energy(self, n_elec):
        '''
        Return the value of energy model for the specified number of electrons,
        i.e. :math:`E(N_{elec})`.

        Parameters
        ----------
        n_elec: float
            Number of electrons.
        '''
        pass

    @abstractmethod
    def energy_derivative(self, n_elec, order=1):
        r'''
        Return the :math:`n^{the}` -order derivative of energy model w.r.t. to the number of electrons
        at fixed chemical potential evaluated for the specified number of electrons.

        .. math::

           {\left. {\left( {\frac{{\partial^{n} E}}{{\partial N^n}}} \right)_{v(r)}} \right|_{N = N_{elec}}}

        Parameters
        ----------
        n_elec: float
            Number of electrons.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{elec} = N_0` the first, second and higher order derivatives are equal
        to the :attr:`BaseGlobalTool.chemical_potential`, :attr:`BaseGlobalTool.chemical_hardness`
        and :attr:`BaseGlobalTool.hyper_hardness`, respectively.
        '''
        pass

    def grand_potential(self, n_elec):
        r'''
        Return the grand potential model defined as:

        .. math:: \Omega = E(\left\langle N \right\rangle) - \mu \cdot \left\langle N \right\rangle

        Parameters
        ----------
        n_elec : float
            Number of electrons.
        '''
        value = self.energy(n_elec) - self.chemical_potential * n_elec
        return value

    def grand_potential_derivative(self, n_elec, order):
        r'''
        Return the :math:`n^{the}` -order derivative of grand potential model w.r.t. to the chemical potential
        at fixed chemical potential evaluated for the specified number of electrons.

        .. math::

           - \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)}
                  = - \left( \frac{\partial^n}{\partial\mu^n} \frac{\partial\Omega}{\partial\mu} \right)_{v(r)}
                  = \left( \frac{\partial^n N}{\partial \mu^n} \right)_{v(r)}

        Parameters
        ----------
        n_elec : float
            Number of electrons.
        order : int, default=1
            The order of derivative denoted by :math:`n` in the formula.

        Note
        ----
        For :math:`N_{elec} = N_0` the second and higher order derivatives are equal
        to the :attr:`chemical_softness`, :attr:`hyper_softness`, respectively.
        '''
        # implementing the Faa di Bruno identity
        # decide whether we use sympy to evaluate Bell polynomilas or code our own Bell polynomial
        pass


class LinearGlobalTool(BaseGlobalTool):
    r'''
    Class of global conceptual DFT reactivity descriptors based on the linear energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` values.

    The energy is approximated as a piece-wise linear function of the number of electrons:

    .. math::

       E(N) = \begin{cases}
                    (N_0 - N) E(N_0 - 1) + (N - (N_0 - 1)) E(N_0) & N < N_0 \\
                    (N_0 + 1 + N) E(N_0 - 1) + (N - N_0) E(N_0 + 1) & N \geqslant N_0 \\
              \end{cases}

    Because the energy model is not differentiable at integer number of electrons, the first derivative of the
    energy w.r.t. number if electrons is calculated from above, from below and averaged:

    .. math::

       \mu^{-} &= -I \\\
       \mu^{0} &= \frac{\mu^{+} + \mu^{-}}{2} \\
       \mu^{+} &= -A \\
    '''
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0)
        self._n_max = self._n0

    @property
    def mu_minus(self):
        r'''
        Chemical potential from below, i.e. :math:`N_0^{-}`, given by:

        .. math:: \mu^{-} = E(N_0) - E(N_0 - 1)  = -I
        '''
        return -1 * self._ip

    @property
    def mu_plus(self):
        r'''
        Chemical potential from above, i.e. :math:`N_0^{+}`, given by:

        .. math:: \mu^{+} = E(N_0 + 1) - E(N_0)  = -A
        '''
        return -1 * self._ea

    @property
    def mu_zero(self):
        r'''
        Chemical potential averaged, given by:

        .. math::

            \mu^{0} = \frac{\mu^{+} + \mu^{-}}{2} = \frac{E(N_0 + 1) - E(N_0 - 1)}{2} = - \frac{I + A}{2}
        '''
        return -0.5 * (self._ea + self._ip)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < self._n0:
            return (self._n0 - n_elec) * self._ip + self._energy_zero
        elif n_elec > self._n0:
            return (self._n0 - n_elec) * self._ea + self._energy_zero
        else:
            return self._energy_zero

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        assert isinstance(order, int) and order > 0
        if order >= 2:
            return None
        if n_elec < self._n0:
            return - self._ip
        elif n_elec > self._n0:
            return - self._ea
        elif n_elec == self._n0:
            return None
        else:
            raise ValueError('')


class QuadraticGlobalTool(BaseGlobalTool):
    r'''
    Class of global conceptual DFT reactivity descriptors based on the quadratic energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy.

    The energy is approximated as a quadratic function of the number of electrons,
    and the three unknown parameters are obtained by fitting to the given values of energy.

    .. math:: E(N) = a + b N + c {N^2}

    First, second and higher order derivatives of the quadratic energy model with respect to
    the number of electrons at fixed external potential are given by:

    .. math::

       {\left( {\frac{{\partial E}}{{\partial N}}} \right)_{v(r)}} &= b + 2 c N \\
       {\left( {\frac{{{\partial ^2}E}}{{\partial {N^2}}}} \right)_{v(r)}} &= 2 c \\
       {\left( {\frac{{{\partial ^n}E}}{{\partial {N^n}}}} \right)_{v(r)}} &= 0 \text{ for } n \geq 2
    '''
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # calculate the a, b, c parameters of the model and N_max
        self._c = 0.5 * (energy_minus - 2 * energy_zero + energy_plus)
        self._b = 0.5 * (energy_plus - energy_minus) - 2 * self._c * n0
        self._a = energy_zero - self._b * n0 - self._c * (n0**2)
        self._n_max = - self._b / (2 * self._c)
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        value = self._a + self._b * n_elec + self._c * n_elec**2
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if not(isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        if order == 1:
            deriv = self._b + 2 * n_elec * self._c
        elif order == 2:
            deriv = 2 * self._c
        elif order >= 2:
            deriv = 0.0
        return deriv


class ExponentialGlobalTool(BaseGlobalTool):
    r'''
    Class of global conceptual DFT reactivity descriptors based on the exponential energy model,
    given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)` known values of energy.

    The energy is approximated as a exponential function of the number of electrons,
    and the three unknown parameters are obtained by fitting to the given values of energy.

    .. math:: E(N) = A \exp(-\gamma(N-N_0)) + B

    The :math:`n^{th}` -order derivative of the rational energy model with respect to the number
    of electrons at fixed external potential is given by:

    .. math::

       {\left( {\frac{{{\partial ^n}E}}{{\partial {N^n}}}} \right)_{v(r)}} =
                A {(-\gamma) ^n} \exp(-\gamma (N - N_0))
    '''
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        # calculate the A, B, gamma parameters of the model and N_max
        self._A = (energy_minus - energy_zero) * (energy_zero - energy_plus)
        self._A /= (energy_minus - 2 * energy_zero + energy_plus)
        self._B = energy_zero - self._A
        self._gamma = math.log(1. - (energy_minus - 2 * energy_zero + energy_plus) / (energy_plus - energy_zero))
        # Calulate N_max
        self._n_max = 0.0
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        dN = n_elec - self._n0
        value = self._A * math.exp(- self._gamma * dN) + self._B
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if not(isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        dN = n_elec - self._n0
        deriv = self._A * math.pow(- self._gamma, order) * math.exp(- self._gamma * dN)
        return deriv


class RationalGlobalTool(BaseGlobalTool):
    r'''
    Class of global conceptual DFT reactivity descriptors based on the 3-parameter
    rational energy model, given :math:`E(N_0 - 1)`, :math:`E(N_0)` and :math:`E(N_0 + 1)`
    known values of energy.

    The energy is approximated as a rational function of the number of electrons,
    and the three unknown parameters are obtained by fitting to the given values of energy.

    .. math:: E(N) = \frac{a_0 + a_1 N}{1 + b_1 N}

    The :math:`n^{th}` -order derivative of the rational energy model with respect to the number
    of electrons at fixed external potential is given by:

    .. math::

       {\left( {\frac{{{\partial ^n}E}}{{\partial {N^n}}}} \right)_{v(r)}} =
            \frac{b_1^{n - 1} (a_1 - a_0 b_1) n!}{(1 + b_1 N)^{2n}}
    '''
    @doc_inherit(BaseGlobalTool)
    def __init__(self, energy_zero, energy_plus, energy_minus, n0):
        #
        # calculate the a0, a1, b1 parameters of the model and N_max
        #
        self._b1 = -1 / (n0 + 1)
        self._a1 = (2 * n0 * energy_plus + (2 * n0 - 1) * energy_zero) / (n0 + 1)
        self._a0 = energy_zero * (1 + self._b1 * n0) - (self._a1 * n0)
        self._n_max = 0
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        value = (self._a0 + self._a1 * n_elec) / (1 + self._b1 * n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if not(isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        deriv = math.pow(self._b1, order - 1) * (self._a1 - self._a0 * self._b1) * math.factorial(order)
        deriv /= math.pow(1 + self._b1 * n_elec, 2 * order)
        return deriv


class GeneralGlobalTool(BaseGlobalTool):
    '''
    Class of global conceptual DFT reactivity descriptors based on the user-specified
    symbolic energy model given known values of energy.

    The energy is approximated as a symbolic function of the number of electrons,
    and the unknown parameters are obtained by fitting to the given values of energy; i.e.
    the number of parameters in the molde should equal the number of given energy values and
    the corresponding number of electrons.

    The :math:`n^{th}` -order derivative of the symbolic energy model with respect to the number
    of electrons at fixed external potential is calculated symbolically.
    '''
    def __init__(self, expr, n0, n_energies, n_symbol, n0_symbol=None, guess={}, opts={}):
        '''
        Parameters
        ----------
        expr : sp.Exp
            The energy expression representing the dependence of energy on the number of electrons.
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`.
        n_energies : dict
            The energy values of `expr` at different electron-numbers.  The dict has int
            (electron-number) keys, float (energy) values.
        n_symbol : sp.Symbol
            The symbol in `expr` that represents the number of electrons.
        n0_symbol: sp.Symbol, optional
            The symbol in `expr` that represents the electron-number at which to evaluate
            `expr`.  If not specified, assume that it is already expressed numerically in
            `expr`.
        guess : dict, optional
            Guesses at the values of the parameters of `expr`.  The dict has sp.Symbol
            keys, float values.
        opts : dict, optional
            Optional keyword arguments to pass to the :py:meth:`scipy.optimimze.root` solver
            that is used to solve for the prameters in the model.
        '''
        # Make sure that the energy expression depends on number of electrons
        if n_symbol not in expr.atoms(sp.Symbol):
            raise ValueError(
                'The expr={0} does not contain {1} symbol representing the number of electrons.'.format(expr, n_symbol))
        self._n_symb = n_symbol

        # Replace the value of N0 in the energy expression, if N0 exists in the expression.
        if n0_symbol:
            expr = expr.subs(n0_symbol, n0)

        # Solve for the parameters of energy expression
        self._params = self._solve(expr, n_energies, guess, opts)
        # Substitute the parameters in the energy expression
        self._expr = expr.subs(self._params.items())

        # solve for N_max
        # the N for which the 1st derivative is equal to zero

        # Calculate the E(N0 - 1), E(N0) and E(N0 + 1)
        energy_zero = self.energy(n0)
        energy_plus = self.energy(n0 + 1)
        energy_minus = self.energy(n0 - 1)
        # Initilaze the base class based on the known energy values
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus, n0)

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        value = self._expr.subs(self._n_symb, n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if not(isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')
        # obtain derivative expression
        deriv = self._expr.diff(self._n_symb, order)
        # evaluate derivative expression at n_elec
        deriv = deriv.subs(self._n_symb, n_elec)
        return deriv

    def _solve(self, expr, n_energies, guess={}, opts={}):
        '''
        Solve for the parameters of the tool's property expression.

        Parameters
        ----------
        See __init__().

        Returns
        -------
        parameters : dict
            A dictionary of sympy.Symbol keys corresponding to the
            value of the expression's solved parameters.
        '''
        # Obtain set of parameters in the energy expression
        params = expr.atoms(sp.Symbol)
        params.remove(self._n_symb)
        if len(params) == 0:
            raise ValueError(
                'There is no parameters in the energy_expression={0} to solve for.'.format(expr))
        if len(params) > len(n_energies):
            print n_energies
            print params
            raise ValueError('Underdetermined system of equations')

        # Set the initial guess of the parameters not specified to 0.
        guess.update({param: 0. for param in params if param not in guess})
        # intial guess for the parameters in the energy model
        guess = np.array([guess[param] for param in params])

        # Construct system of equations to solve
        system_eqns = []
        for n, energy in n_energies.iteritems():
            eqn = sp.lambdify((params,), expr.subs(self._n_symb, n) - energy, 'numpy')
            system_eqns.append(eqn)

        def objective(args):
            '''
            Evaluate the system of equations for the given values of parameters.

            Parameters
            ----------
            args : array representing the value of parameters.
                The expression for the property.
            '''
            return np.array([eqn(args) for eqn in system_eqns])

        #TODO: Add jacobina

        # solve for the parameters in the energy model
        result = root(objective, guess, **opts)
        if not result.success:
            raise ValueError('The system of equations could not be solved. message:{0}'.format(result.message))
        # make the dictionary of parameter values
        parameters = dict([(param, result.x[i]) for i, param in enumerate(params)])

        return parameters
