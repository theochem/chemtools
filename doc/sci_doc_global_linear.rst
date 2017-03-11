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


.. _linear_energy:

Linear Energy Model :class:`chemtools.tool.globaltool.LinearGlobalTool`
=======================================================================

In this model, the energy is approximated as a piece-wise linear function of the number of electrons:

 .. math:: E(N) = a + b N

 .. TODO::
    Technically, linear model has two parameters, so providing two E values is enough to fit the model.

 .. math::

    \text{For } N \leq N_0: E\left(N\right) &= a + b N \\
         a &= E\left(N_0\right) - N_0 \left(E\left(N_0\right) - E\left(N_0 - 1\right)\right) \\
         b &= E\left(N_0\right) - E\left(N_0 - 1\right)

 .. math::

    \text{For } N \geq N_0: E\left(N\right) &= a + b N \\
         a &= E\left(N_0\right) - N_0 \left(E\left(N_0 + 1\right) - E\left(N_0\right)\right) \\
         b &= E\left(N_0 + 1\right) - E\left(N_0\right)

The model requires three values of :math:`E(N)` to interpolate the energy. Commonly, the energies of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given data points results in two equations:

 .. math::

    E\left(N\right) &= \begin{cases}
             \left(N_0 - N\right) E\left(N_0 - 1\right) + \left(N - \left(N_0 - 1\right)\right) E\left(N_0\right) & \text{ for } N < N_0 \\
	     \left(N_0 + 1 + N\right) E\left(N_0 - 1\right) + \left(N - N_0\right) E\left(N_0 + 1\right) & \text{ for } N \geqslant N_0 \\
	    \end{cases} \\

or equivalently,

 .. math::

    E\left(N\right) &= \begin{cases}
	     E\left(N_0\right) + \left(N_0 - N\right) \cdot IP & \text{ for } N < N_0 \\
	     E\left(N_0\right) + \left(N_0 - N\right) \cdot EA & \text{ for } N \geqslant N_0 \\
	    \end{cases} \\

At this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.LinearGlobalTool.energy`.

The energy model is not differentiable at integer number of electrons, so the chemical potential
is not defined. Instead one calculates the chemical potential from above, below and averaged:

 .. math::

    \mu^{+} &= \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}^+ = -EA \\
    \mu^{-} &= \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}^- = -IP \\
    \mu^{0} &= \frac{\mu^{+} + \mu^{-}}{2} = \frac{-\left(IP + EA\right)}{2} \\

In this model, second and higher order derivatives of the energy with respect to the number of electrons is zero
(they are not defined at integer number of electrons, because energy model is not differentiable).
So, chemical hardness and hyper-hardness are zero, and softness and hyper-softness are not defined.

 .. TODO::
    Is it better to skip derived global tools for this model?
    How the code should handle these?

To obtain the :ref:`derived global reactivity indicators <derived_indicators>` for
the linear energy model, the maximum number of electrons accepted by the system should be calculated.
This is obtained by setting the first order derivative of energy equal to zero, however, in this model
the first derivative of energy is not defined.

The related :ref:`derived global reactivity indicators <derived_indicators>` for the linear energy model are:

 .. math::

    \omega_{\text {electrophilicity}} &= E\left(N_0\right) - E\left(N_{\text max}\right) &&= 0 \\
    \omega_{\text {nucleophilicity}} &= ? \\
    \nu_{\text {nucleofugality}} &= E\left(N_0 + 1\right) - E\left(N_{\text max}\right)
                                &&=  \\
    \nu_{\text {electrofugality}} &= E\left(N_0 - 1\right) - E\left(N_{\text max}\right)
