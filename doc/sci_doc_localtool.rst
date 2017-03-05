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


.. _local_tools:

Local Descriptive Tools :class:`chemtools.tool.localtool`
#########################################################

Local descriptive tool :math:`p (\mathbf{r})` assigns a value to every point in space.
These tools are used to study the reactivity of a molecule at point :math:`\mathbf{r}`.

In conceptual DFT, local descriptors arise as functional derivatives of :ref:`global descriptive tools <global_tools>`
with respect to local quantities, typically the external potential :math:`v(\mathbf{r})`, at fixed number of
electrons :math:`N`. Here, :ref:`linear <local_linear>` and :ref:`quadratic <local_quadratic>` energy models are considered,
and the corresponding local descriptors will be introduced.

The functional derivatives of :ref:`energy and its derivatives <energy_derivatives>`
with respect to external potential :math:`v(\mathbf{r})` at fixed number of
electrons :math:`N` results in **electron density** :math:`\rho(\mathbf{r})`,
**Fukui function** :math:`f(\mathbf{r})`, :math:`2^{\text{nd}}` **-order Fukui function** :math:`f^{(2)}(\mathbf{r})`
commonly refered to as **dual descriptor** :math:`\Delta f(\mathbf{r})`
and :math:`n^{\text{th}}` **-order Fukui function** :math:`f^{(n)}(\mathbf{r})` as local descriptors:

 .. math::

    \rho(\mathbf{r}) = f^{(0)}(\mathbf{r}) &= {\left( \frac{\delta E}{\delta v(\mathbf{r})} \right)_N} && \\
    f(\mathbf{r}) = f^{(1)}(\mathbf{r}) &= {\left( \frac{\delta \mu}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}} \right)_N} \\
    \Delta f(\mathbf{r}) = f^{(2)}(\mathbf{r}) &= {\left( \frac{\delta \eta}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}} \right)_N} \\
    f^{(n)}(\mathbf{r}) &= {\left( \frac{\delta \eta^{(n-1)}}{\delta v(\mathbf{r})} \right)_N}
             &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                 {\left( \frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})}} \right)_N}

 .. todo::
    * What are the names for the derivatives of :math:`\Omega` listed below?

On the other hand, the functional derivative of :ref:`grand potential and its derivatives <grand_potential_derivatives>`
with respect to external potential
:math:`v(\mathbf{r})` at fixed chemical potential :math:`\mu` results in **electron density** :math:`\rho(\mathbf{r})`, **local softness** :math:`s(\mathbf{r})`,
and :math:`2^{\text{nd}}` **-order local softness** :math:`s^{(2)}(\mathbf{r})`, and
:math:`n^{\text{th}}` **-order local softness** :math:`s^{(n)}(\mathbf{r})`
as local descriptors for a grand canonical ensemble:

 .. math::

    \rho(\mathbf{r}) = s^{(0)}(\mathbf{r}) &= {\left( \frac{\delta \Omega}{\delta v(\mathbf{r})} \right)_{\mu}}  \\
    s(\mathbf{r}) = s^{(1)}(\mathbf{r}) &= {\left( \frac{\delta S}{\delta v(\mathbf{r})} \right)_{\mu}}
              = {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial \Omega}{\partial \mu} \right)_{v(\mathbf{r})}} \right)_{\mu}}
                = S \cdot f(\mathbf{r})  \\
    s^{(2)}(\mathbf{r}) &= {\left( \frac{\delta S}{\delta v(\mathbf{r})} \right)_{\mu}}
              = {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial^2 \Omega}{\partial {\mu}^2} \right)_{v(\mathbf{r})}} \right)_{\mu}}
		= S^{2} \cdot f^{(2)}(\mathbf{r}) + S^{(2)} \cdot f(\mathbf{r}) \\
    s^{(n)}(\mathbf{r}) &= {\left( \frac{\delta S^{(n-1)}}{\delta v(\mathbf{r})} \right)_{\mu}}
             = {\left( \frac{\delta}{\delta v(\mathbf{r})}
                 {\left( \frac{\partial^n \Omega}{\partial {\mu}^n} \right)_{v(\mathbf{r})}} \right)_{\mu}} \\
               &= -\sum_{k=1}^n f^{(k)}(\mathbf{r}) \cdot B_{n,k}\left(S^{(1)}, S^{(2)}, ..., S^{(n-k+1)} \right)  \\

For details of expressing local softness expressions in terms of Fukui functions using Faà di Bruno's formula,
please refer to :ref:`derivation_local_softness`.

To derive a working expression for these local quantities, the functional
derivative of :ref:`energy and its derivatives <energy_derivatives>` with respect to external potential
:math:`v(\mathbf{r})` at fixed number of electrons :math:`N` should be calculated.
Considering that,

 .. math::

    \frac{\delta E_{N_0 \pm k}}{\delta v(\mathbf{r})} = \rho_{N_0 \pm k}(\mathbf{r})

for a general energy model the Fukui function expressions can be derived using chain-rule:

 .. math::

    f^{(n)}(\mathbf{r}) &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n E_{\text{model}}}
                           {\partial N^n}\right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n}
       {\partial N^n} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\frac{\delta}
       {\delta v(\mathbf{r})} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_N} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}} \frac{\delta \alpha_i}{\delta v(\mathbf{r})}
       \right)_N} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \left( \frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}} \cdot \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}} \frac{\delta E_{N_0 \pm k}}{\delta v(\mathbf{r})}
       \right)_N\right)} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \left( \left.\frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}}\right|_{N=N_0} \cdot \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}} \rho_{N_0 \pm k}(\mathbf{r})
       \right)\right)} \right)_{v(\mathbf{r})}} \\
    &= \sum_{i=1}^n \left(\left. \frac{\partial^{n+1} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial N^n \partial\alpha_{i}} \right|_{N=N_0} \cdot
       \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}}
       \rho_{N_0 \pm k}(\mathbf{r})\right)_{v(\mathbf{r})}

..    &= \sum_{i=1}^n \left(\left. \underbrace {\frac{\partial^{n+1} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial N^n \partial\alpha_{i}}}_{\text {evaluated analytically}}\right|_{N=N_0} \cdot
       \sum_k \underbrace {\frac{\partial \alpha_i}{\partial E_{N_0 \pm k}}}_{\text{evaluated}\atop\text{numerically}}
       \rho_{N_0 \pm k}(\mathbf{r})\right)_{v(\mathbf{r})}


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
and :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)` when setting :math:`N` equal to :math:`N_0-1`, :math:`N_0` and  :math:`N_0+1`, respectively.
Also, integrating the linear electron density over all space results in :math:`N`
which confirms that the density expression is properly normalized to the number of electrons.

By rearranging the expression for linear electron density, it can be easily perceived as first-order Taylor expanssion of density
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


.. _local_quadratic:

Quadratic Energy Model :class:`chemtools.tool.localtool.QuadraticLocalTool`
===========================================================================

Complementing the :ref:`quadratic global tools <quadratic_energy>`, the fitted quadratic energy expression

 .. math::

    E\left(N\right) = E\left(N_0\right) &+ \left(\frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2}\right) \left(N - N_0\right) \\
                  &+ \left(\frac{E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)}{2}\right) \left(N - N_0\right)^2

and its first and second derivatives with respect to :math:`N` at fixed external potential,

 .. math::

    \mu\left(N\right) &= \left(\frac{\partial E\left(N\right)}{\partial N}\right)_{v(\mathbf{r})} \\
      &= \left(\frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2}\right) +
         \left(E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)\right) \left(N - N_0\right) \\
    \eta\left(N\right) &= \left(\frac{\partial^2 E\left(N\right)}{\partial^2 N}\right)_{v(\mathbf{r})} \\
      &= E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)

alongside the electron density of systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons, namely
:math:`{\{\rho_{N_0 - 1}\left(\mathbf{r}\right), \rho_{N_0}\left(\mathbf{r}\right), \rho_{N_0 + 1}\left(\mathbf{r}\right)\}}`,
are used to calculate quadratic local descriptors. These local tools include:

 #. :ref:`Quadratic Electron Density <quadratic_density>`
 #. :ref:`Quadratic Fukui Function <quadratic_fukui_function>`
 #. :ref:`Quadratic Dual Descriptor <quadratic_dual_descriptor>`
 #. :ref:`Quadratic Softness <quadratic_softness>`
 #. :ref:`Quadratic Hyper-Softness <quadratic_hyper_softness>`


.. _quadratic_density:

**Quadratic Electron Density:** According to Eq. ???, the change in quadratic energy expression with respect to external
potential at fixed number of electrons yields the density of :math:`N`-electron system, that is,

 .. math::

    \rho_{N}(\mathbf{r}) =& \left( \frac{\delta E\left(N\right)}{\delta v(\mathbf{r})} \right)_N \\
     =& \left( \frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N &&+
	\left[\left( \frac{\delta E\left(N_0 +1\right)}{\delta v(\mathbf{r})} \right)_N -
	      \left( \frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N \right] \frac{N - N_0}{2} \\
      & &&+ \left[\left(\frac{\delta E\left(N_0 + 1\right)}{\delta v(\mathbf{r})} \right)_N - 2
                          \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N +
	                  \left(\frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N \right] \frac{\left(N - N_0\right)^2}{2}

So,

 .. math::
    \rho_{N}(\mathbf{r}) = \rho_{N_0}\left(\mathbf{r}\right)
     &+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right) \\
     &+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right)^2

As expected, the obtained expression for density equals :math:`\rho_{N_0 - 1}\left(\mathbf{r}\right)`, :math:`\rho_{N_0}\left(\mathbf{r}\right)`
and :math:`\rho_{N_0 + 1}\left(\mathbf{r}\right)` when setting :math:`N` equal to :math:`N_0-1`, :math:`N_0` and  :math:`N_0+1`, respectively.
Also, integrating the quadratic electron density over all space gives :math:`N` which confirms that the density expression properly integrates to the number of electrons.

It is important to note that the obtained expression for quadratic electron density in Eq. ???  can
be perceived as the second-order Taylor expansion of density around the :math:`\rho_{N_0}(\mathbf{r})`
as the reference within the quadratic energy model, that is,

 .. math::

    \rho_{N}(\mathbf{r}) = \rho_{N_0}\left(\mathbf{r}\right) +
         \left. \left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \right|_{N = N_0} \left(N - N_0\right) + \frac{1}{2}
         \left. \left(\frac{\partial^2 \rho_{N}(\mathbf{r})}{\partial N^2}\right)_{v(\mathbf{r})} \right|_{N = N_0} \left(N - N_0\right)^2

where,

 .. math::

    \left. \left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \right|_{N = N_0} &=
         \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) = f_{N_0}(\mathbf{r})  \\
    \left. \left(\frac{\partial^2 \rho_{N}(\mathbf{r})}{\partial N^2}\right)_{v(\mathbf{r})} \right|_{N = N_0} &=
         \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right) = \Delta f_{N_0}(\mathbf{r})

The first and second derivatives of electron density with respect to number of electrons evaluated at :math:`N_0` are the quadratic
Fukui function and dual descriptor of :math:`N_0`-electron system.
The quadratic density of :math:`N`-electron system is implemented in :class:`chemtools.tool.localtool.QuadraticLocallTool.density`.


.. _quadratic_fukui_function:

**Quadratic Fukui Function:** According to Eq. ???, the change in quadratic chemical potential with respect to external potential at fixed number of electrons yields
the Fukui function of :math:`N`-electron system. Equivalently, the Fukui function of :math:`N`-electron system can be viewed as the
change in quadratic electron density :math:`\rho_N\left(\mathbf{r}\right)` with respect to number of electrons :math:`N` at fixed external potential.
In other words,

 .. math::

    f_{N}(\mathbf{r}) = \left( \frac{\delta \mu\left(N\right)}{\delta v(\mathbf{r})} \right)_N =
                        \left( \frac{\delta}{\delta v(\mathbf{r})} \left(\frac{\partial E\left(N\right)}{\partial N}\right)_{v(\mathbf{r})} \right)_N =
			\left( \frac{\partial}{\partial N} \left(\frac{\delta E\left(N\right)}{\delta v(\mathbf{r})}\right)_{N} \right)_{v(\mathbf{r})}  =
			\left(\frac{\partial \rho_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} \\

where,

 .. math::

    \left( \frac{\delta \mu\left(N\right)}{\delta v(\mathbf{r})} \right)_N = \frac{1}{2}
         && \left[\left( \frac{\delta E\left(N_0 +1\right)}{\delta v(\mathbf{r})} \right)_N -
                  \left( \frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N \right] + \\
         && \left[\left(\frac{\delta E\left(N_0 + 1\right)}{\delta v(\mathbf{r})} \right)_N - 2
              \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N +
	      \left(\frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N \right] \left(N - N_0\right) \\

    \left( \frac{\partial \rho_{N}(\mathbf{r})}{\partial N} \right)_{v(\mathbf{r})} = \frac{\partial}{\partial N}
         \rho_{N_0}\left(\mathbf{r}\right)
      &&+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right) \\
      &&+ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}\right) \left(N - N_0\right)^2

**(Note: Fix the missing bracket in the last expression)**
Simplifying either of the above expressions results in the quadratic Fukui function of :math:`N`-electron system:

 .. math::

    f_{N}(\mathbf{r}) = \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2} \right) +
	\left[\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right) \right] \left(N - N_0\right)

Integrating the Fukui function expression confirms that it is normalized to one for any number of electrons :math:`N`.
For :math:`N=N_0`, the familiar expression of Fukui function in obtained:

 .. math::

    f_{N_0}\left(\mathbf{r}\right) = \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2} \\

It is important to note that the obtained expression for quadratic Fukui function in Eq. ???  can be perceived as the first-order Taylor expansion
of Fukui function around the :math:`f_{N_0}(\mathbf{r})` as the reference within the quadratic energy model, that is,

 .. math::

    f_{N}(\mathbf{r}) = f_{N_0}\left(\mathbf{r}\right) + \left. \left(\frac{\partial f_N(\mathbf{r})}{\partial N}\right) \right|_{N = N_0} \left(N - N_0\right)

where,

 .. math::

    \left. \left(\frac{\partial f_N(\mathbf{r})}{\partial N}\right) \right|_{N = N_0} =
    \rho_{N_0 - 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 + 1}\left(\mathbf{r}\right) = \Delta f_{N_0}(\mathbf{r})

the derivative of Fukui function with respect to the number of electrons is the dual descriptor.
The quadratic Fukui function is implemented in :class:`chemtools.tool.localtool.QuadraticLocallTool.fukui_function`.


.. _quadratic_dual_descriptor:

**Quadratic Dual Descriptor:** According to Eq. ???, the change in quadratic chemical hardness with respect to external potential at fixed number of electrons yields
the dual descriptor of :math:`N`-electron system. Equivalently, the dual descriptor of :math:`N`-electron system can be viewed as the
change in quadratic Fukui function with respect to the number of electrons :math:`N` at fixed external potential. That is,

 .. math::

    \Delta f_{N}(\mathbf{r}) = \left( \frac{\delta \eta\left(N\right)}{\delta v(\mathbf{r})} \right)_N =
                        \left( \frac{\delta}{\delta v(\mathbf{r})} \left(\frac{\partial \mu\left(N\right)}{\partial N}\right)_{v(\mathbf{r})} \right)_N =
			\left( \frac{\partial}{\partial N} \left(\frac{\delta \mu\left(N\right)}{\delta v(\mathbf{r})}\right)_{N} \right)_{v(\mathbf{r})}  =
			\left(\frac{\partial f_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})}

where,

 .. math::

    \left( \frac{\delta \eta\left(N\right)}{\delta v(\mathbf{r})} \right)_N &=
        \left(\frac{\delta E\left(N_0 + 1\right)}{\delta v(\mathbf{r})} \right)_N - 2
        \left(\frac{\delta E\left(N_0\right)}{\delta v(\mathbf{r})} \right)_N +
	\left(\frac{\delta E\left(N_0 - 1\right)}{\delta v(\mathbf{r})} \right)_N \\
    \left(\frac{\partial f_{N}(\mathbf{r})}{\partial N}\right)_{v(\mathbf{r})} &=
        \frac{\partial}{\partial N} \left[ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2} \right) +
	\left(\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right) \right) \left(N - N_0\right) \right]

Simplifying either of the above expressions results in the quadratic dual descriptor of :math:`N` -electron system:

 .. math::

    \Delta f_{N}(\mathbf{r}) = \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)

The dual descriptor does not depend on :math:`N` as one expects for the quadratic energy model.
Also, the obtained expression properly integrates to zero.
The dual descriptor is implemented in :class:`chemtools.tool.localtool.QuadraticLocallTool.dual_descriptor`.

**Mention that higher order local descriptors do not exist.**


.. _quadratic_softness:

**Quadratic Softness:** The quadratic local softness is easily found by substituting the quadratic Fukui functions in Eq. (????):

 .. math::

    s_N\left(\mathbf{r}\right) &= S \cdot f_N\left(\mathbf{r}\right) \\
      &= \frac{1}{\eta} \cdot \left[ \left(\frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2} \right) +
         \left(\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right) \right) \left(N - N_0\right) \right]

For :math:`N=N_0`,

 .. math::

     s_{N_0}\left(\mathbf{r}\right) = \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2 \eta} =
     \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)}


.. _quadratic_hyper_softness:

**Quadratic HyperSoftness:**

 .. math::

    s^{(2)}\left(\mathbf{r}\right) &= S^{2} \cdot f^{(2)}(\mathbf{r}) + S^{(2)} \cdot f(\mathbf{r}) \\
     &= \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
        \rho_{N_0 - 1}\left(\mathbf{r}\right)}{\eta^2} =
        \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
        \rho_{N_0 - 1}\left(\mathbf{r}\right)}{\left[E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)\right]^2} \\


Analytical
==========

Here the analytical evaluation of fukui function, dual descriptor, etc. will be described!


