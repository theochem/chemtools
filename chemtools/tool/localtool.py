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
"""Local Conceptual Density Functional Theory (DFT) Reactivity Tools.

   This module contains various local tool classes corresponding to
   linear and quadratic energy models.
"""


from chemtools.utils import doc_inherit


class BaseLocalTool(object):
    """
    Base class of local conceptual DFT reactivity descriptors.
    """
    def __init__(self, density_zero, density_plus, density_minus, n0):
        r"""
        Parameters
        ----------
        density_zero : np.ndarray
            Electron density of :math:`N_0`-electron system, i.e.
            :math:`\rho_{N_0}\left(\mathbf{r}\right)`.
        density_plus : np.ndarray
            Electron density of :math:`(N_0 + 1)`-electron system, i.e.
            :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)`.
        density_minus : np.ndarray
            Electron density of :math:`(N_0 - 1)`-electron system, i.e.
            :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`.
        n0 : float
            Reference number of electrons, i.e. :math:`N_0`, which corresponds
            to the integral of density_zero over all space.
        """
        # if np.any(density_zero < 0):
        #     raise ValueError('Argument density_zero should be all positive!')
        # if np.any(density_plus < 0):
        #     raise ValueError('Argument density_plus should be all positive!')
        # if np.any(density_minus < 0):
        #     raise ValueError('Argument density_minus should be all positive!')
        self._density_zero = density_zero
        self._density_plus = density_plus
        self._density_minus = density_minus
        self._n0 = n0

    @property
    def n0(self):
        r"""
        Reference number of electrons, i.e. :math:`N_0`, corresponding to density_zero.
        """
        return self._n0

    @property
    def density_zero(self):
        r"""
        Electron density of :math:`N_0`-electron system, i.e.
        :math:`\rho_{N_0}\left(\mathbf{r}\right)`.
        """
        return self._density_zero

    @property
    def density_plus(self):
        r"""
        Electron density of :math:`(N_0 + 1)`-electron system, i.e.
        :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)`.
        """
        return self._density_plus

    @property
    def density_minus(self):
        r"""
        Electron density of :math:`(N_0 - 1)`-electron system, i.e.
        :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`.
        """
        return self._density_minus


class LinearLocalTool(BaseLocalTool):
    r"""
    Class of local conceptual DFT reactivity descriptors based on the linear energy model.

    Considering the fitted linear energy expression,

    .. math::
       E\left(N\right) =
        \begin{cases}
         \left(N - N_0 + 1\right) E\left(N_0\right) - \left(N - N_0\right) E\left(N_0 - 1\right)
            &  N \leqslant N_0 \\
         \left(N - N_0\right) E\left(N_0 + 1\right) - \left(N - N_0 - 1\right) E\left(N_0\right)
            &  N \geqslant N_0 \\
        \end{cases} \\

    and its derivative with respect to the number of electrons at constant external potential,

    .. math::
       \mu\left(N\right) =
        \begin{cases}
         \mu^- &= E\left(N_0\right) - E\left(N_0 - 1\right) = - IP &&  N < N_0 \\
         \mu^0 &= 0.5 \left(E\left(N_0 + 1\right) - E\left(N_0 - 1\right)\right) = -0.5 (IP + EA)
               && N = N_0 \\
         \mu^+ &= E\left(N_0 + 1\right) - E\left(N_0\right) = - EA &&  N > N_0 \\
        \end{cases}

    the linear local tools are obtained by taking the functional derivative of these expressions
    with respect to external potential :math:`v(\mathbf{r})` at fixed number of electrons.
    """
    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus, n0):
        super(self.__class__, self).__init__(density_zero, density_plus, density_minus, n0)
        self._ff_plus = self._density_plus - self._density_zero
        self._ff_minus = self._density_zero - self._density_minus
        self._ff_zero = 0.5 * (self._density_plus - self._density_minus)

    @property
    def ff_plus(self):
        r"""
        Fukui Function from above defined as,

        .. math::
           f^+\left(\mathbf{r}\right) = \rho_{N_0 + 1}\left(\mathbf{r}\right) -
                                        \rho_{N_0}\left(\mathbf{r}\right)
        """
        return self._ff_plus

    @property
    def ff_minus(self):
        r"""
        Fukui Function from below define as,

        .. math::
           f^-\left(\mathbf{r}\right) = \rho_{N_0}\left(\mathbf{r}\right) -
                                        \rho_{N_0 - 1}\left(\mathbf{r}\right)
        """
        return self._ff_minus

    @property
    def ff_zero(self):
        r"""
        Fukui Function from center defined as the average of :attr:`ff_plus` and :attr:`ff_minus`,

        .. math::
           f^0\left(\mathbf{r}\right) =
           \frac{f^+\left(\mathbf{r}\right) + f^-\left(\mathbf{r}\right)}{2} =
           \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}
        """
        return self._ff_zero

    def density(self, number_electrons=None):
        r"""
        Linear electron density of :math:`N`-electron system defined as the functional derivative of
        linear energy model w.r.t. external potential at fixed number of electrons, i.e.,

        .. math::
           \rho_{N}(\mathbf{r}) =
            \begin{cases}
             \rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r})
                       \right] \left(N - N_0\right) & \text{ for } N \leqslant N_0 \\
             \rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r})
                       \right] \left(N - N_0\right) & \text{ for } N \geqslant N_0 \\
            \end{cases}

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`\rho_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        rho = self._density_zero.copy()
        if (number_electrons is not None) and (number_electrons != self._n0):
            if number_electrons < self._n0:
                rho += self._ff_minus * (number_electrons - self._n0)
            elif number_electrons > self._n0:
                rho += self._ff_plus * (number_electrons - self._n0)
        return rho

    def fukui_function(self, number_electrons=None):
        r"""
        Linear Fukui function of :math:`N`-electron system defined as the functional derivative of
        linear chemical potential w.r.t. external potential at fixed number of electrons,

        .. math::
           f_{N}(\mathbf{r}) =
             \begin{cases}
               f^-(\mathbf{r}) &= \rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r}) && N < N_0 \\
               f^0\left(\mathbf{r}\right) &= 0.5 \left(\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                        \rho_{N_0 - 1}\left(\mathbf{r}\right)\right) && N = N_0 \\
               f^+(\mathbf{r}) &= \rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r}) && N > N_0 \\
             \end{cases}

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`f^0\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._ff_zero
        elif number_electrons < self._n0:
            return self._ff_minus
        elif number_electrons > self._n0:
            return self._ff_plus
        else:
            raise ValueError('Argument number_electrons={0} is not valid.'.format(number_electrons))

    def softness(self, global_softness, number_electrons=None):
        r"""
        Linear softness of :math:`N`-electron system defined as,

        .. math::
           s_N\left(\mathbf{r}\right) = S \cdot f_N\left(\mathbf{r}\right) =
             \begin{cases}
               S \cdot f^-(\mathbf{r}) & N < N_0 \\
               S \cdot f^0\left(\mathbf{r}\right) & N = N_0 \\
               S \cdot f^+(\mathbf{r}) &  N > N_0 \\
             \end{cases}

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`S \cdot f^0\left(\mathbf{r}\right)`
            is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return global_softness * self._ff_zero
        elif number_electrons < self._n0:
            return global_softness * self._ff_minus
        elif number_electrons > self._n0:
            return global_softness * self._ff_plus
        else:
            raise ValueError('Argument number_electrons={0} is not valid.'.format(number_electrons))


class QuadraticLocalTool(BaseLocalTool):
    r"""
    Class of local conceptual DFT reactivity descriptors based on the quadratic energy model.

    Considering the fitted quadratic energy expression,

    .. math::
       E\left(N\right) = E\left(N_0\right) &+ \left(\frac{E\left(N_0 + 1\right) -
                         E\left(N_0 - 1\right)}{2}\right) \left(N - N_0\right) \\
                &+ \left(\frac{E\left(N_0 - 1\right) - 2 E\left(N_0\right) +
                   E\left(N_0 + 1\right)}{2}\right) \left(N - N_0\right)^2 \\

    and its first and second derivatives with respect to the number of electrons at constant
    external potential,

    .. math::
       \mu\left(N\right) &= \left(\frac{\partial E\left(N\right)}{\partial N}
                            \right)_{v(\mathbf{r})} \\
         &= \left(\frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2}\right) +
            \left[E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)
            \right] \left(N - N_0\right) \\
       \eta\left(N\right) &= \left(\frac{\partial^2 E\left(N\right)}{\partial^2 N}
                             \right)_{v(\mathbf{r})}
         = E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)

    the quadratic local tools are obtained by taking the functional derivative of these expressions
    with respect to external potential :math:`v(\mathbf{r})` at fixed number of electrons.
    """
    @doc_inherit(BaseLocalTool)
    def __init__(self, density_zero, density_plus, density_minus, n0):
        super(self.__class__, self).__init__(density_zero, density_plus, density_minus, n0)
        # Fukui function and dual descriptor of N0-electron system
        self._ff0 = 0.5 * (self._density_plus - self._density_minus)
        self._df0 = self._density_plus - 2 * self._density_zero + self._density_minus

    def density(self, number_electrons=None):
        r"""
        Quadratic electron density of :math:`N`-electron system defined as the functional
        derivative of quadratic energy model w.r.t. external potential at fixed number of
        electrons,

        .. math::
           \rho_{N}(\mathbf{r}) = \rho_{N_0}\left(\mathbf{r}\right)
             &+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right) \\
             &+ \left(\frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
                \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)}{2}\right)
                \left(N - N_0\right)^2

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`\rho_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._density_zero
        else:
            dN = (number_electrons - self._n0)
            rho = self._density_zero + self._ff0 * dN + 0.5 * self._df0 * dN**2
            return rho

    def fukui_function(self, number_electrons=None):
        r"""
        Quadratic Fukui function of :math:`N`-electron system defined as the functional
        derivative of quadratic chemical potential w.r.t. external potential at fixed number
        of electrons,

        .. math::
           f_{N}(\mathbf{r}) = \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) -
                 \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2} \right) +
                 \left[\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
                 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)
                 \right] \left(N - N_0\right)

        Parameters
        ----------
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`f_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        if (number_electrons is None) or (number_electrons == self._n0):
            return self._ff0
        else:
            ff = self._ff0 + self._df0 * (number_electrons - self.n0)
            return ff

    def dual_descriptor(self):
        r"""
        Quadratic dual descriptor of :math:`N`-electron system defined as the functional
        derivative of quadratic chemical hardness w.r.t. external potential at fixed number
        of electrons,

        .. math::
           \Delta f_{N}(\mathbf{r}) = \rho_{N_0 - 1}\left(\mathbf{r}\right) - 2
            \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right)

        The quadratic dual descriptor is independent of the number electrons.
        """
        return self._df0

    def softness(self, global_softness, number_electrons=None):
        r"""
        Quadratic softness of :math:`N`-electron system defined as,

        .. math::
           s_N\left(\mathbf{r}\right) = S \cdot f_N\left(\mathbf{r}\right)

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        number_electrons : float, default=None
            Number of electrons. If None, the :math:`s_{N_0}\left(\mathbf{r}\right)` is returned.
        """
        s_value = global_softness * self.fukui_function(number_electrons)
        return s_value

    def hyper_softness(self, global_softness):
        r"""
        Quadratic hyper-softness of :math:`N`-electron system defined as,

        .. math::
           s_N^{(2)}\left(\mathbf{r}\right) = S^2 \cdot \Delta f_N\left(\mathbf{r}\right)

        The quadratic hyper-softness is independent of the number electrons.

        Parameters
        ----------
        global_softness : float
            The value of global softness.
        """
        s_value = global_softness**2 * self.dual_descriptor()
        return s_value
