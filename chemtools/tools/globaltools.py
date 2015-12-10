#!/usr/bin/env python
'''Global Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


import math


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
    def __init__(self, ip, ea):
        '''
        Parameters
        ----------
        ip : float
            The ionization potential.
        ea : float
            The electron affinity.
        '''
        super(self.__class__, self).__init__(ip, ea)

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
