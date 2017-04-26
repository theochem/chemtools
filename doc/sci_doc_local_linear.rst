..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2014-2015 The ChemTools Development Team
    :
    : This file is part of ChemTools.
    :
    : ChemTools is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : ChemTools is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --


.. _local_linear:

Linear Energy Model :class:`chemtools.tool.localtool.LinearLocalTool`
=====================================================================

Complementing the :ref:`linear global tools <linear_energy>`, the fitted piece-wise linear
energy expression,

 .. math::

    E\left(N\right) = \begin{cases}
             \left(N - N_0 + 1\right) E\left(N_0\right) - \left(N - N_0\right) E\left(N_0 - 1\right) & \text{ for } N \leqslant N_0 \\
	     \left(N - N_0\right) E\left(N_0 + 1\right) - \left(N - N_0 - 1\right) E\left(N_0\right) & \text{ for } N \geqslant N_0 \\
	    \end{cases}

and its derivative with respect to number of electrons :math:`N` at fixed external potential,

 .. math::

    \mu\left(N\right) = \begin{cases}
             E\left(N_0\right) - E\left(N_0 - 1\right) = - IP &= \mu^- & \text{ for } N < N_0 \\
	     E\left(N_0 + 1\right) - E\left(N_0\right) = - EA &= \mu^+ & \text{ for } N > N_0 \\
	    \end{cases}

alongside the electron density of systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons, namely
:math:`{\{\rho_{N_0 - 1}\left(\mathbf{r}\right), \rho_{N_0}\left(\mathbf{r}\right), \rho_{N_0 + 1}\left(\mathbf{r}\right)\}}`,
are used to calculate linear local descriptors. These local tools include:

 #. :ref:`Linear Electron Density <linear_density>`
 #. :ref:`Linear Fukui Function <linear_fukui_function>`
 #. :ref:`Linear Softness <linear_softness>`


.. _linear_density:

**Linear Electron Density:** According to Eq. ???, the change in linear energy expression with respect to external
potential at fixed number of electrons yields the density of :math:`N`-electron system, that is,

 .. math::

    \rho_{N}(\mathbf{r}) = \left( \frac{\delta E\left(N\right)}{\delta v(\mathbf{r})} \right)_N =
      \begin{cases}
        \left(N - N_0 + 1\right) \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N -
	\left(N - N_0\right) \left(\frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N & \text{ for } N \leqslant N_0 \\
	\left(N - N_0\right) \left(\frac{\delta E\left(N_0 + 1\right)}{\delta v(\mathbf{r})} \right)_N -
	\left(N - N_0 - 1\right) \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N & \text{ for } N \geqslant N_0 \\
      \end{cases}

So,

 .. math::

    \rho_{N}(\mathbf{r}) =
      \begin{cases}
        \left(N - N_0 + 1\right) \rho_{N_0}(\mathbf{r}) - \left(N - N_0\right) \rho_{N_0 - 1}(\mathbf{r}) & \text{ for } N \leqslant N_0 \\
	\left(N - N_0\right) \rho_{N_0 + 1}(\mathbf{r}) - \left(N - N_0 - 1\right) \rho_{N_0}(\mathbf{r}) & \text{ for } N \geqslant N_0 \\
      \end{cases}

As expected, the obtained expression for density equals :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`, :math:`\rho_{N_0}\left(\mathbf{r}\right)`
and :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)` when setting :math:`N` equal to :math:`\left(N_0-1\right)`, :math:`N_0` and
:math:`\left(N_0+1\right)`, respectively.
Also, integrating the linear electron density over all space results in :math:`N`
which confirms that the density expression is properly normalized to the number of electrons.

By rearranging the expression for linear electron density, it can be easily perceived as first-order Taylor expansion of density
around the :math:`\rho_{N_0}(\mathbf{r})` as the reference within the linear energy model, that is,

 .. math::

    \rho_{N}(\mathbf{r}) =
      \begin{cases}
        \rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r})\right] \left(N - N_0\right) & \text{ for } N \leqslant N_0 \\
	\rho_{N_0}(\mathbf{r}) + \left[\rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r})\right] \left(N - N_0\right) & \text{ for } N \geqslant N_0 \\
      \end{cases}

where,

 .. math::

    \left. \left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \right|_{N = N_0^-} &=
         \rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r}) = f_{N_0}^-(\mathbf{r}) \\
    \left. \left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \right|_{N = N_0^+} &=
         \rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r}) = f_{N_0}^+(\mathbf{r})


.. _linear_fukui_function:

**Linear Fukui Function:** According to Eq. ???, the change in linear chemical potential with respect to external potential at fixed number of electrons yields
the Fukui function of :math:`N`-electron system. Equivalently, the Fukui function of :math:`N`-electron system can be viewed as the
change in linear electron density :math:`\rho_N\left(\mathbf{r}\right)` with respect to number of electrons :math:`N` at fixed external potential.
In other words,

 .. math::

    f_{N}(\mathbf{r}) = \left( \frac{\delta \mu\left(N\right)}{\delta v(\mathbf{r})} \right)_N =
                        \left( \frac{\delta}{\delta v(\mathbf{r})} \left(\frac{\partial E\left(N\right)}{\partial N}\right)_{v(\mathbf{r})} \right)_N =
			\left( \frac{\partial}{\partial N} \left(\frac{\delta E\left(N\right)}{\delta v(\mathbf{r})}\right)_{N} \right)_{v(\mathbf{r})}  =
			\left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \\

where,

 .. math::

    \left( \frac{\delta \mu\left(N\right)}{\delta v(\mathbf{r})} \right)_N &=
      \begin{cases}
        \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N -
	\left(\frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N & \text{ for } N < N_0 \\
	\left(\frac{\delta E\left(N_0 + 1\right)}{\delta v(\mathbf{r})} \right)_N -
	\left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N & \text{ for } N > N_0 \\
      \end{cases} \\
    \left( \frac{\partial \rho_{N}(\mathbf{r})}{\partial N} \right)_{v(\mathbf{r})} &=
      \begin{cases}
        \frac{\partial}{\partial N}\left[\left(N - N_0 + 1\right) \rho_{N_0}(\mathbf{r}) - \left(N - N_0\right) \rho_{N_0 - 1}(\mathbf{r})\right] & \text{ for } N < N_0 \\
	\frac{\partial}{\partial N}\left[\left(N - N_0\right) \rho_{N_0 + 1}(\mathbf{r}) - \left(N - N_0 - 1\right) \rho_{N_0}(\mathbf{r})\right] & \text{ for } N > N_0 \\
      \end{cases}

Simplifying either of the above expressions results in the linear Fukui function of :math:`N`-electron system:

 .. math::

    f_{N}(\mathbf{r}) =
      \begin{cases}
        \rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r}) = f^-(\mathbf{r}) & \text{ for } N < N_0 \\
	\rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r}) = f^+(\mathbf{r}) & \text{ for } N > N_0 \\
      \end{cases}

Considering the fact that the linear energy model is not differentiable at :math:`N_0`,
Commonly, the average Fukui function :math:`f^0\left(\mathbf{r}\right)` is used:

 .. math::

    f^0\left(\mathbf{r}\right) = \frac{f^+\left(\mathbf{r}\right) + f^-\left(\mathbf{r}\right)}{2} =
             \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}

**(TO BE REMOVED)** Dual descriptor is defined as:

 .. math::

    d\left(\mathbf{r}\right) = f^+\left(\mathbf{r}\right) - f^-\left(\mathbf{r}\right) =
           \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0 - 1}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)

 .. todo::
    * This is not really dual descriptor for linear model. Technically the dual descriptor is zero for linear model,
      but the dual descriptor for quadratic model happens to be f+(r) - f-(r).
      Does this need to be clarified?


.. _linear_softness:

**Linear Softness:**

 .. math::

    s_{N}(\mathbf{r}) = S \cdot f_{N}(\mathbf{r}) =
      \begin{cases}
        S \cdot \left[\rho_{N_0}(\mathbf{r}) - \rho_{N_0 - 1}(\mathbf{r})\right] = S \cdot f^-(\mathbf{r}) & \text{ for } N < N_0 \\
	S \cdot \left[\rho_{N_0 + 1}(\mathbf{r}) - \rho_{N_0}(\mathbf{r})\right] = S \cdot f^+(\mathbf{r}) & \text{ for } N > N_0 \\
      \end{cases}
