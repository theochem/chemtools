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
"""Conceptual Density Functional Theory (DFT) Reactivity Tools Based on General Energy Model.

This module contains the global and local tool classes corresponding to user-specified energy
models.
"""

import numpy as np
import sympy as sp
from horton import log
from scipy.optimize import root, least_squares
from chemtools.utils.utils import doc_inherit
from chemtools.conceptual.base import BaseGlobalTool

__all__ = ['GeneralGlobalTool']


class GeneralGlobalTool(BaseGlobalTool):
    r"""
    Class of global conceptual DFT reactivity descriptors based on the user-specified energy model.

    The energy model is approximated as a symbolic function of the number of electrons,
    and the unknown parameters are obtained by interpolating to the given values of energy; i.e.
    the number of parameters in the model should equal the number of given energy values and
    the corresponding number of electrons.

    The :math:`n^{\text{th}}`-order derivative of the symbolic energy model with respect to the
    number of electrons at fixed external potential is calculated symbolically.
    """

    def __init__(self, expr, n0, n_energies, n_symbol=None, n0_symbol=None, guess=None, opts=None):
        """
        Initialize class.

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
        """Parameter dictionary of energy model."""
        return self._params

    @property
    def n_symbol(self):
        """Symbol used to denote the number of electrons."""
        return self._n_symb

    @property
    def expression(self):
        """Energy expression as a function of number of electrons, :math:`E(N)`."""
        return self._expr

    @doc_inherit(BaseGlobalTool)
    def energy(self, n_elec):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n_min <= n_elec <= self._n_max:
            log.warn('Energy evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n_min, self._n_max))
        # evaluate energy
        value = self._expr.subs(self._n_symb, n_elec)
        return value

    @doc_inherit(BaseGlobalTool)
    def energy_derivative(self, n_elec, order=1):
        if n_elec < 0.0:
            raise ValueError('Number of electrons cannot be negativ! n_elec={0}'.format(n_elec))
        if not self._n0 - 1 <= n_elec <= self._n0 + 1:
            log.warn('Energy derivative evaluated for n_elec={0} outside of interpolation '
                     'region [{1}, {2}].'.format(n_elec, self._n0 - 1, self._n0 + 1))
        if not (isinstance(order, int) and order > 0):
            raise ValueError('Argument order should be an integer greater than or equal to 1.')

        # obtain derivative expression
        deriv = self._expr.diff(self._n_symb, order)
        # evaluate derivative expression at n_elec
        deriv = deriv.subs(self._n_symb, n_elec)
        return deriv

    def _solve_parameters(self, expr, n_energies, guess, opts=None):
        r"""
        Solve for the unknown parameters of the energy model.

        The parameters of energy model is found using the given :math:`\left(N, E(N)\right)` pairs.
        This requires solving a system of (non-)linear equations

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
            raise ValueError('The expr={0} does not contain parameters given in guess={1}'
                             ''.format(expr, guess))
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
            Evaluate jacobian of the objective function for the given values of parameters.

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
        r"""Solve for the :math:`N_{\text{max}}` of the energy model."""
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
                log.warn('The system of equations for Nmax could not be solved; Nmax=`None`. '
                         'message:{0}'.format(result.message))
        return n_max
