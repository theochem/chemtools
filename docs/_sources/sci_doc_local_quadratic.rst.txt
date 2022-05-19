..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2016-2019 The ChemTools Development Team
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


.. _local_quadratic:

Quadratic Energy Model
======================

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

