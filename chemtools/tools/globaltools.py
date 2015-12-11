#!/usr/bin/env python
'''Global Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


import math
import numpy as np
import sympy as sp
from scipy.optimize import root


class BaseGlobalTool(object):
    '''
    Base Class of Global Conceptual DFT Reactivity Descriptors.
    '''
    def __init__(self, energy_zero, energy_plus, energy_minus):
        '''
        Parameters
        ----------
        energy_zero : float
            The energy of inital state
        energy_plus : float
            The energy of initial state with one more electron
        energy_minus : float
            The energy of inital state with one less electron
        '''
        self._energy_zero = energy_zero
        self._energy_minus = energy_minus
        self._energy_plus = energy_plus
        self._ip = energy_minus - energy_zero
        self._ea = energy_zero - energy_plus

    @property
    def energy_zero(self):
        '''Energy of inital state'''
        return self._energy_zero

    @property
    def ip(self):
        '''Ionization Potential (IP).'''
        return self._ip

    @property
    def ionization_potential(self):
        '''Ionization Potential (IP).'''
        return self._ip

    @property
    def electron_affinity(self):
        '''Electron Affinity (EA).'''
        return self._ea

    @property
    def ea(self):
        '''Electron Affinity (EA).'''
        return self._ea


class QuadraticGlobalTool(BaseGlobalTool):
    '''
    Class of Global Conceptual DFT Reactivity Descriptors based on the Quadratic Energy Model.
    '''
    def __init__(self, energy_zero, energy_plus, energy_minus):
        '''
        Parameters
        ----------
        energy_zero : float
            The energy of inital state
        energy_plus : float
            The energy of initial state with one more electron
        energy_minus : float
            The energy of inital state with one less electron
        '''
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus)

    @property
    def mu(self):
        '''
        Chemical potential defined as the first derivative of quadratic energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu  = {\left( {\frac{{\partial E}}{{\partial N}}} \right)_{v(r)}} =  - \frac{{I + A}}{2}$
        '''
        return -0.5 * (self._ip + self._ea)

    @property
    def chemical_potential(self):
        '''
        Chemical potential defined as the first derivative of the quadratic energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu  = {\left( {\frac{{\partial E}}{{\partial N}}} \right)_{v(r)}} =  - \frac{{I + A}}{2}$
        '''
        return self.mu

    @property
    def eta(self):
        '''
        Chemical hardness defined as the second derivative of the quadratic energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu  = {\left( {\frac{{{\partial ^2}E}}{{\partial {N^2}}}} \right)_{v(r)}} = I - A$
        '''
        return self._ip - self._ea

    @property
    def chemical_hardness(self):
        '''
        Chemical hardness defined as the second derivative of the quadratic energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu  = {\left( {\frac{{{\partial ^2}E}}{{\partial {N^2}}}} \right)_{v(r)}} = I - A$
        '''
        return self.eta

    @property
    def softness(self):
        '''Chemical softness.'''
        value = 1.0 / self.eta
        return value

    @property
    def electronegativity(self):
        '''Mulliken Electronegativity.'''
        value = -1 * self.mu
        return value

    @property
    def electrophilicity(self):
        '''Electrophilicity.'''
        value = math.pow(self.mu, 2) / (2 * self.eta)
        return value

    @property
    def n_max(self):
        '''N_max value.'''
        value = - self.mu / self.eta
        return value

    @property
    def nucleofugality(self):
        '''Nucleofugality.'''
        value = math.pow(self._ip - 3 * self._ea, 2)
        value /= (8 * (self._ip - self._ea))
        return value

    @property
    def electrofugality(self):
        '''Electrofugality.'''
        value = math.pow(3 * self._ip - self._ea, 2)
        value /= (8 * (self._ip - self._ea))
        return value


class LinearGlobalTool(BaseGlobalTool):
    '''
    Class of Global Conceptual DFT Reactivity Descriptors based on the Linear Energy Model.
    '''
    def __init__(self, energy_zero, energy_plus, energy_minus):
        '''
        Parameters
        ----------
        ip : float
            The ionization potential.
        ea : float
            The electron affinity.
        '''
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus)

    @property
    def mu_minus(self):
        '''
        Chemical potential defined as the first derivative from below of the linear energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu^{-} = -I$
        '''
        return -1 * self._ip

    @property
    def mu_plus(self):
        '''
        Chemical potential defined as the first derivative from above of the linear energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu^{+} = -A$
        '''
        return -1 * self._ea

    @property
    def mu_zero(self):
        '''
        'Chemical potential defined as the averaged first derivative of the linear energy model w.r.t.
        the number of electrons at fixed external potential.
        $\mu^{0} = \frac{\mu^{+} + \mu^{-}}{2}$
        '''
        return -0.5 * (self._ea + self._ip)


class GeneralizedGlobalTool(object):
    '''
    Generalized, Symbolic Class of Global Conceptual DFT Reactivity Descriptors.
    '''

    def __init__(self, expr, n0, n_energies, n_symbol, **kwargs):
        '''
        Initialize the GeneralizedGlobalTool instance.

        Parameters
        ----------
        expr : sp.Expr
            The energy expression.
        n0 : int
            The electron number at which to evaluate `expr`.
        n_energies : dict
            The energy values of `expr` at different electron-numbers.  The dict has int
            (electron-number) keys, float (energy) values.
        n_symbol : sp.Symbol
            The symbol in `expr` that represents electron-number.
        n0_symbol: sp.Symbol, optional
            The symbol in `expr` that represents the electron-number at which to evaluate
            `expr`.  If not specified, assume that it is already expressed numerically in
            `expr`.
        guess : dict, optional
            Guesses at the values of the parameters of `expr`.  The dict has sp.Symbol
            keys, float values.
        opts : dict, optional
            Optional keyword arguments to pass to the internal SciPy solver.

        '''

        # Handle keyword arguments
        defaults = {'n0_symbol': None, 'guess': None, 'opts': None}
        defaults.update(kwargs)
        n0_symbol = defaults['n0_symbol'] if defaults['n0_symbol'] else None
        guess = defaults['guess'] if defaults['guess'] else None
        opts = defaults['opts'] if defaults['opts'] else None

        # Assign basic attributes
        self._n0 = n0
        self._n_symbol = n_symbol
        self._n0_symbol = n0_symbol

        # Solve for the parameters of the given `expr`, update attributes
        substitutions = [(self._n_symbol, n0)]
        if self._n0_symbol:
            substitutions.append((self._n0_symbol, n0))
        solution = self._solve(expr, n_energies, n0_symbol, guess, opts)
        self.expr = solution[0]
        self.d_expr = solution[1]
        self._params = solution[2]
        self._mu = self.d_expr.subs(substitutions)
        self._eta = self.d_expr.diff(self._n_symbol).subs(substitutions)

    def __getattr__(self, attr):
        '''
        Allow any attribute corresponding to higher-order derivatives of energy
        (e.g., `gen_global.hyper_eta_2`), to be created on the fly.

        Parameters
        ----------
        attr : str
            The name of the attribute Python is attempting to access.

        Returns
        -------
        hyper_eta_x : sp.Expr
            The 'x'th derivative of the energy expression wrt electron-number.

        Raises
        ------
        AssertionError
            If the order of derivative specified by 'x' is invalid.

        '''

        order = int(attr.rsplit('_', 1)[-1])
        assert order >= 0, \
            "The order of derivative cannot be negative."
        if order == 0:
            return self._mu
        elif order == 1:
            return self._eta
        else:
            return self._eta.diff(self._n_symbol, order).subs(self._n_symbol, self._n0)

    def _solve(self, expr, n_energies, n0_symbol=None, guess=None, opts=None):
        '''
        Solve for the parameters of the tool's energy expression.

        See __init__().

        Returns
        -------
        solution : tuple of sp.Expr, dict
            The solved expression, the solved expression's first derivative wrt
            electron-number, and a dictionary of sp.Symbol keys corresponding to the
            value of the expression's solved parameters.

        Raises
        ------
        AssertionError
            If the system if underdetermined or if the expression cannot be solved.

        '''

        # Initialize variables
        expr_free = expr.copy()
        if not guess:
            guess = {}
        if not opts:
            opts = {}

        # Parse the non-electron-number parameters of `expr`
        if n0_symbol:
            expr = expr.subs(n0_symbol, self._n0)
        params = [symbol for symbol in expr.atoms() if isinstance(symbol, sp.Symbol) and
                  #(symbol is not self._n_symbol) and (symbol is not n0_symbol)]
                  symbol not in (self._n_symbol, n0_symbol)]

        # Fill in missing non-electron-number symbols from `guess`
        guess.update({symbol: 0. for symbol in params if symbol not in guess})

        # Construct system of equations to solve for `params`
        assert len(params) <= len(n_energies), \
            "The system is underdetermined.  Inlcude more (elec-number, energy) pairs."
        system_eqns = []
        for n, energy in n_energies.iteritems():
            eqn = sp.lambdify((params,), expr.subs(self._n_symbol, n) - energy, "numpy")
            system_eqns.append(eqn)

        # Solve the system of equations for `params`
        def objective(args):
            return np.array([fun(args) for fun in system_eqns])

        root_guess = np.array([guess[symbol] for symbol in params])
        roots = root(objective, root_guess, **opts)
        assert roots.success, \
            "The system of equations could not be solved."

        # Substitute in the values of `params` and differentiate wrt `n_symbol`
        expr = expr_free.subs([(params[i], roots.x[i]) for i in range(len(params))])
        d_expr = expr.diff(self._n_symbol)
        param_dict = {params[i]: roots.x[i] for i in range(len(params))}

        # Return
        return expr, d_expr, param_dict

    def energy(self, n0):
        '''
        Return the energy expression evaluated at electron-number `n0`.

        Parameters
        ----------
        n0 : int
            The electron number.

        Returns
        -------
        value : sp.Expr
            The energy expression with the value of `n0` substituted in for
            electron-number `n`.

        '''

        value = self.expr.subs(self._n_symbol, n0)
        if self._n0_symbol:
            value = value.subs(self._n0_symbol, n0)
        return value

    def d_energy(self, n0, order=1):
        '''
        Return the derivative of the energy expression wrt electron-number
        evaluated at electron-number `n0`.

        Parameters
        ----------
        n0 : int
            The electron number.
        order : int, optional
            The order of the derivative wrt electron-number to take.

        Returns
        -------
        value : sp.Expr
            The energy expression with the value of `n0` substituted in for
            electron-number `n`.

        Raises
        ------
        AssertionError
            If the order of derivative specified by 'x' is invalid.

        '''

        assert order >= 0, \
            "The order of the derivative cannot be negative."
        value = self.d_expr.diff(self._n_symbol, order - 1).subs(self._n_symbol, n0)
        if self._n0_symbol:
            value = value.subs(self._n0_symbol, n0)
        return value

    @property
    def mu(self):
        return self._mu

    @property
    def eta(self):
        return self._eta

    @property
    def params(self):
        return self._params


class ExpGlobalTool(BaseGlobalTool):
    '''
    Class of Global Conceptual DFT Reactivity Descriptors based on the Exponential Energy Model.

    The energy model has the form:
    ..math::
        E(N) = A \exp(-\gamma(N-N_0)) + B

    We can solve for :math: `\gamma`, :math:`A`, and :math:`B` with the conditions that
    ..math::
        E(N_0-1) &= A \exp(\gamma) + B\\
        E(N_0) &= A + B\\
        E(N_0+1) &= A \exp(-\gamma) + B\\

    Solving this,
    ..math::
        B &= E(N_0) A\\
        E(N_0+1) - E(N_0) &= A \exp(-\gamma) - A\\
        &= (A-1) \exp(-\gamma)\\
        A &= 1+ \frac{E(N_0+1)-E(N_0)}}{\exp(-\gamma)}\\
        &= E(N_0-1)-E(N_0) &= (A-1)\exp(\gamma)\\
        &= (E(N_0+1)) - E(N_0) \frac{\exp(\gamma)}{\exp(-\gamma)}\\
        &= (E(N_0+1)) - E(N_0) \exp(2\gamma)\\
        \exp(2\gamma) &= \frac{E(N_0-1)-E(N_0)}{E(N_0+1)-E(N_0)}\\
        \gamma &= \frac{1}{2} \ln \left( \frac{E(N_0-1)-E(N_0)}{E(N_0+1)-E(N_0)} \right)\\
        &= \frac{1}{2} \ln \frac{IP}{EA}\\
        A &= 1 + \sqrt{(E(N_0+1)-E(N_0))(E(N_0)-1)(E(N_0))}\\
        &= 1+ \sqrt{(EA)(IP)}
        B &= E(N_0) - 1 - \sqrt{(E(N_0+1)-E(N_0))(E(N_0)-1)(E(N_0))}\\
        &= 1 - A

    '''
    def __init__(self, energy_zero, energy_plus, energy_minus, num_elec_zero):
        super(self.__class__, self).__init__(energy_zero, energy_plus, energy_minus)
        self._num_elec_zero = num_elec_zero

    @property
    def num_elec_zero(self):
        """ Number of electrons of the system at ground state
        """
        return self._num_elec_zero
    
    @property
    def gamma(self):
        """ First parameter of exponential energy model

        Note
        ----
        See class documentation
        """
        return 0.5 * math.log(self.ip / self.ea)

    @property
    def A(self):
        """ Second parameter of epxonential energy model

        Note
        ----
        See class documentation
        """
        return 1 + math.sqrt(self.ip * self.ea)

    @property
    def B(self):
        """ Third parameter of exponential energy model

        Note
        ----
        See class documentation
        """
        return self.ground_state - self.A

    def energy(self, num_elec):
        """Energy for system with arbitrary number of electrons

        ..math::
            E(N) = A \exp(-\gamma(N-N_0)) + B

        Parameter
        ---------
        num_elec: float
            Number of electrons

        Returns
        -------
        energy: float
        """
        return self.A * math.exp(-self.gamma * (num_elec - self.num_elec_ground)) + self.B

    def grand_potential(self, num_elec):
        """Grand potential for system with arbitrary number of election
        
        ..math::
            \omega = E - \miu * N
        Parameters
        ----------
        num_elec : float
            Number of electrons
        
        Returns
        -------
        name : float
        """
        return self.energy(num_elec) - self.chemical_potential * num_elec

    def energy_parital_der(self, order):
        """ Returns the nth order partial derivative of the energy with respect
        to the number of electrons evaluated at the ground state number of electrons

        ..math::
            E(N) = A \exp(-\gamma(N-N_0)) + B
            \frac{\partial}{\partial N} E(N) = A (-\gamma) \exp(-\gamma(N-N_0)) = \mu
            \frac{\partial^2}{\partial N^2} E(N) = A (-\gamma)^2 \exp(-\gamma(N-N_0)) = \eta
            \frac{\partial^n}{\partial N^n} E(N) = A (-\gamma)^n \exp(-\gamma(N-N_0)) = \eta^{(n-1)}

        Parameters
        ----------
        order: int
            Order of the partial differentiation

        Returns
        -------
        hyperhardness: float
            nth order hyperhardness

        Raises
        ------
        AssertionError
            If order is negative
        """
        assert order >= 0
        if order > 0:
            return self.A * (-self.gamma)**order
        else:
            return self.A + self.B

    @property
    def chemical_potential(self):
        """ Returns the chemical potential of the system using the Exponential
        energy model at the number of electrons at ground state

        ..math::
            \frac{\partial}{\partial N} E(N_0) = A (-\gamma)
        """
        return self.energy_partial_der(1)

    @property
    def mu(self):
        """ Chemical potential
        """
        return self.chemical_potential

    @property
    def chemical_hardness(self):
        """ Returns the chemical hardness of the system using the Exponential
        energy model at the number of electrons at ground state

        ..math::
            \frac{\partial^2}{\partial N^2} E(N_0) = A (-\gamma)^2
        """
        return self.energy_partial_der(2)

    @property
    def eta(self):
        """ Chemical hardness
        """
        return self.chemical_hardness

    @property
    def softness(self):
        """ Returns chemical softness

        ..math::
            (\frac{\partial^2 \Omega}{\\partial \mu^2})_{v(\ver{r})} &= \frac{1}{\eta}\\
            &= S
        """
        return 1./self.chemical_hardness
