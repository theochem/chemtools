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


.. _quadratic_energy:

Quadratic Energy Model :class:`chemtools.tool.globaltool.QuadraticGlobalTool`
=============================================================================

In this model, the energy is approximated as a quadratic function of the number of electrons:

 .. math::
    :nowrap:
    :label: quadratic

    \begin{eqnarray}
     E\left(N\right) = a + b N + c {N^2}
    \end{eqnarray}

As it contains three parameters, :math:`a`, :math:`b` and :math:`c`, this model requires
three values of :math:`E\left(N\right)` to interpolate the energy. Commonly, the energies of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations,

 .. math::

    \begin{cases}
          E\left(N_0 - 1\right) &= a + b \left(N_0 - 1\right) + c {\left(N_0 - 1\right) ^2} \\
             E \left(N_0\right) &= a + b \left(N_0\right) + c {\left(N_0\right) ^2} \\
          E\left(N_0 + 1\right) &= a + b \left(N_0 + 1\right) + c {\left(N_0 + 1\right) ^2}
    \end{cases}

which can be solved for the three unknowns,

 .. math::

    a &= E\left(N_0\right) - b N_0 - c {N_0 ^2} \\
    b &= \frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2} - 2 N_0 c \\
    c &= \frac{E\left(N_0 - 1\right) -2 E\left(N_0\right) + E\left(N_0 + 1\right)}{2} \\

Substituting the obtained parameters :math:`a`, :math:`b` and :math:`c` into the energy expression,
Eq. :eq:`quadratic`, gives the fitted energy model as:

 .. math::

    E\left(N\right) = E\left(N_0\right) &+ \left(\frac{E\left(N_0 + 1\right) - E\left(N_0 - 1\right)}{2}\right) \left(N - N_0\right) \\
                  &+ \left(\frac{E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)}{2}\right) \left(N - N_0\right)^2

or equivalently,

 .. math::

    E\left(N\right) = E\left(N_0\right) - \left(\frac{IP + EA}{2}\right) \left(N - N_0\right) + \left(\frac{IP - EA}{2}\right) \left(N - N_0\right)^2

At this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.energy`. By rearranging
the obtained quadratic energy expression, the energy change :math:`\Delta E = E(N) - E(N_0)` due to
the electron transfer :math:`\Delta N = N - N_0`, when the external potential :math:`v(\mathbf{r})`
is fixed, is given by:

 .. math::

    \Delta E = -\left(\frac{IP + EA}{2}\right) \Delta N + \left(\frac{IP - EA}{2}\right) (\Delta N)^2

As detailed below, the prefactor of :math:`\Delta N` is the first derivative of energy with respect to :math:`N`
and the prefactor of :math:`(\Delta N)^2` is one-half the second order derivatives of the energy with
respect to :math:`N` at fixed external potential
:math:`v(\mathbf{r})` evaluated at :math:`N = N_0`. As a result, this energy model is equivalent
to the second-order Taylor expansion of the energy as a function of :math:`N` around the reference
state :math:`N_0`.

To obtain the :ref:`fundamental global reactivity indicators <global_fundamental_indicators>` for the
quadratic energy model, the derivatives of the energy with respect to the number of electrons at
fixed external potential :math:`v(\mathbf{r})` should be calculated. These are given by:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}
         &= b + 2cN \\
	 &= \frac{E(N_0 + 1) - E(N_0 - 1)}{2} + \left(E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1)\right) \left(N - N_0\right) \\
	 &= -\frac{IP + EA}{2} + (IP - EA) \left(N - N_0\right) \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}
         &= 2c \\
	 &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) \\
	 &= IP - EA \\
    \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})}
         &= 0 \text{   for   } n \geq 2

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.energy_derivative`.
In the quadratic model, evaluating the first-, second-, and higher-order derivatives of
energy evaluated at :math:`N_0` gives the following expressions for the chemical potential,
chemical hardness, and hyper-hardnesses,

 .. math::

    \mu = \left. \left(\frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0}
       &= \frac{E(N_0 + 1) - E(N_0 - 1)}{2}  = - \frac{{IP + EA}}{2} \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
        &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) = IP - EA \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
              &= 0 \text{   for   } n \geq 2

These are implemented in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_hardness`.

Accordingly, within the quadratic energy model, the chemical softness and hyper-softnesses are given
by the expressions,

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} = \frac{1}{IP - EA} \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= 0 \text {     for } n \geq 2

To obtain the :ref:`derived global reactivity indicators <global_derived_indicators>` for
the quadratic energy model, the maximum number of electrons to saturate the system should be calculated.
This is obtained by setting the first derivative of the energy with respect to the number of electrons equal
to zero,

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} = 0 &= b + 2cN = -\frac{IP + EA}{2} + (IP - EA)(N - N_0) \\
    & \to N_{\text{max}} = \frac{-b}{2c} = N_{0} + \frac{IP + EA}{2 \left(IP - EA \right)} = N_{0} - \frac{\mu}{\eta} \\
    & \to \Delta N_{\text{max}} = N_0 - N_{\text{max}} = \frac{IP + EA}{2 \left(IP - EA \right)} = - \frac{\mu}{\eta}

The related :ref:`derived global reactivity indicators <global_derived_indicators>` for the quadratic energy model are:

 .. todo:: include the generalized signed definitions.

 .. math::

    \omega_{\text{electrophilicity}} &= E\left(N_0\right) - E\left(N_{\text{max}}\right)
                        &&= \frac{\left(IP + EA\right)^2}{8\left(IP - EA\right)}
		       &&&= \frac{\mu^2}{2 \eta} \\
    \nu_{\text{nucleofugality}} &= E\left(N_0 + 1\right) - E\left(N_{\text{max}}\right)
                                &&= \frac{\left(IP - 3 \cdot EA \right)^2}{8 \left(IP - EA \right)}
			       &&&=  \frac{\left(\mu + \eta\right)^2}{2\eta} = -EA + \omega_{\text{electrophilicity}} \\
    \nu_{\text{electrofugality}} &= E\left(N_0 - 1\right) - E\left(N_{\text{max}}\right)
                                 &&= \frac{\left(3 \cdot IP - EA \right)^2}{8 \left(IP - EA \right)}
				&&&= \frac{\left(\mu - \eta\right)^2}{2\eta} = IP + \omega_{\text{electrophilicity}}



**References:**

.. bibliography:: ../data/references.bib
   :style: unsrt
   :start: continue
   :list: bullet
   :filter: key % "Parr1983JACS"
