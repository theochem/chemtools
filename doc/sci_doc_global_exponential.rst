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


.. _exponential_energy:

Exponential Energy Model :class:`chemtools.tool.globaltool.ExponentialGlobalTool`
=================================================================================

In this model, energy is approximated by an exponential function of the number of electrons:

 .. math::

    E(N) = A \exp(-\gamma(N-N_0)) + B

Containing three parameters, :math:`A`, :math:`B` and :math:`\gamma`, this model requires
three values of :math:`E(N)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations:

 .. math::

    \begin{cases}
          E(N_0-1) &= A \exp(\gamma) + B \\
          E(N_0)   &= A + B \\
          E(N_0+1) &= A \exp(-\gamma) + B \\
    \end{cases}

This allows us to solve for the three unknonws:

 .. math::

    A      &= \frac{(E\left(N_0 - 1\right) - E\left(N_0\right))(E\left(N_0\right) - E\left(N_0 + 1\right))}
                   {E\left(N_0 - 1\right) - 2 E\left(N_0\right) + E\left(N_0 + 1\right)}
            = \frac{IP \cdot EA}{IP - EA} \\
    B      &= E\left(N_0\right) - A  \\
    \gamma &= \ln \left( 1 - \frac{E\left(N_0 - 1\right) - 2E\left(N_0\right) + E\left(N_0 + 1\right)}
                                  {E\left(N_0 + 1\right) - E\left(N_0\right)} \right) \\

Due to the complexity of the obtained parameters, we skip substituting them into the energy expression.
However, knowing the parameters of the model, at this stage the energy expression can be evaluated for
any given number of electrons as implemented in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.energy`.

The derivatives of the energy model with respect to the number of electrons at
fixed external potential :math:`v(\mathbf{r})` are:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}
         &= A \left(-\gamma\right) \exp\left(-\gamma \left(N - N_0\right)\right) \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}
         &= A {\left(-\gamma\right)^2} \exp\left(-\gamma \left(N - N_0\right)\right) \\
    \left( \frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})}
         &= A {\left(-\gamma\right)^n} \exp\left(-\gamma \left(N - N_0\right)\right) \text{   for   } n \geq 1

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0}
       &= -A \gamma \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
        &= A {\gamma ^2} \\
    \eta^{(2)} = \left. \left( \frac{\partial^{3} E}{\partial N^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
              &= -A \gamma^3 \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
              &= A {(-\gamma)^{(n+1)}} \text{  for  } n \geq 2

These are implemented in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.ExponentialGlobalTool.chemical_hardness`.

Accordingly, given the exponential energy model, chemical softness and :math:`2^{\text{nd}}` and
:math:`3^{\text{rd}}` -order hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &=  \frac{1}{\eta} = \frac{1}{A \gamma^2} \\
    S^{(2)} = - \left. \left( \frac{\partial^{3}\Omega}{\partial\mu^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 = - \left(-A\gamma^3\right) \left(\frac{1}{A \gamma^2} \right)^3 = \frac{1}{A^2\gamma^3} \\
    S^{(3)} = - \left. \left( \frac{\partial^{4}\Omega}{\partial\mu^{4}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
	   &= - \left(A\gamma^4\right) \left(\frac{1}{A\gamma^2}\right)^4 +
	      3 \left(\frac{1}{-A\gamma^3}\right)^2 \left(\frac{1}{A\gamma^2}\right)^5 = \frac{-4}{A^3\gamma^4}\\

The higher order hyper-softness exists and can be evaluated through Eq. ???, as implemented in
:meth:`chemtools.tool.globaltool.ExponentialGlobalTool.hyper_softness`.

To obtain the :ref:`derived global reactivity indicators <derived_indicators>` for
the exponential energy model, the maximum number of electrons accepted by the system should be calculated.

 .. TODO::
    #. Write down the value of N_max and derived global reactivity tools

**References:**

 .. TODO::
    #. Add references
