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

.. _global_tools:

Global Descriptive Tools :mod:`chemtools.tool.globaltool`
#########################################################

Global descriptive tools assign a single value to the entire molecule. In other words,
the property value :math:`P_{\text{global}}` describes the intrinsic reactivity of the molecule
as a whole.

Very simple examples of such global reactivity descriptors are the first **ionization potential**
(:math:`IP`) and **electron affinity** (:math:`EA`) of the system which allow one to measure its
propensity to donate or accept one electron. For a system with :math:`N_0` electrons,
these are defined as,

 .. math::
    IP = E\left(N_0 - 1\right) - E\left(N_0\right) \\
    EA = E\left(N_0\right) - E\left(N_0 + 1\right)

 .. TODO::
    Are there any other global reactivity descriptors independent from conceptual DFT?

The global reactivity indicators from conceptual DFT are either :ref:`fundamental <fundamental_indicators>`
or :ref:`derived <derived_indicators>`. Obtaining these global tools requires selecting
an energy model, :math:`E_{\text{model}} = E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})`,
representing the dependence of energy on the
number electrons :math:`N`. The set :math:`{\{\alpha_1, \alpha_2, ..., \alpha_n\}}` denotes the parameters
of the energy model which are determined by fitting the energy expression to the known values of
the energy for :math:`n` different numbers of electrons,
:math:`{\{E(N_i)\}}_{i=1}^n`. Commonly, the system with :math:`N_0` electrons denotes the reference state,
and the values of energy for systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1`,
i.e. :math:`{\{E(N_0 - 1), E(N_0), E(N_0 + 1)\}}`, are used to fit the parameters in the energy model.
However, other values of energy can be used to parametrized the model as well. Needless to say, the number of
required energy values to solve for the parameters depends on the complexity of the energy model,
i.e. the number of parameters in the model.

The implemented energy models and the associated global reactivity descriptors in ChemTools include:

 #. :ref:`Linear Energy Model <linear_energy>`
 #. :ref:`Quadratic Energy Model <quadratic_energy>`
 #. :ref:`Exponential Energy Model <exponential_energy>`
 #. :ref:`Rational Energy Model <rational_energy>`
 #. :ref:`General Energy Model <general_energy>`

.. _fundamental_indicators:

**Fundamental Global Reactivity Descriptors**

.. _energy_derivatives:

In the canonical ensemble, the fundamental global reactivity descriptors include the derivatives
of the energy model with respect to the number of electrons :math:`N` at fixed external potential
:math:`v(\mathbf{r})`:

 .. math:: \left( \frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}{\partial N} \right)_{v(\mathbf{r})}

More specifically, the fundamental global reactivity indicators of the :math:`N_0` electron reference state
are the first, second and higher order derivatives of energy
evaluated at :math:`N=N_{0}`, which are called **chemical potential** denoted by :math:`\mu`,
**chemical hardness** denoted by :math:`\eta`, and :math:`n^{\text{th}}` **-order hyper-hardness**
denoted by :math:`\eta ^{(n)} \text{for } n \geq 2`, respectively:

 .. math::

    \mu \equiv \left. \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \\
    \eta \equiv \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \\
    \eta^{(2)} \equiv \left. \left( \frac{\partial^{3} E}
                {\partial N^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0} \\
    \eta^{(n)} \equiv \left. \left( \frac{\partial^{n+1} E}
                {\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \text{for } n \geq 2

 .. TODO::
    #. Which one is correct: :math:`\equiv` or :math:`=`?
    #. Mention electronegativity; is it negative chemical potential for all energy models?

.. _grand_potential_derivatives:

In the grand canonical ensemble, the fundamental global reactivity descriptors include the derivatives of the grand
potential model :math:`\Omega = E \left(\left\langle N \right\rangle\right) - \mu \left\langle N \right\rangle`
with respect to the chemical potential :math:`\mu` at fixed external potential :math:`v(\mathbf{r})`.

 .. TODO::
    #. Elaborate on grand potential

 .. math::

    - \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})}
          = - \left( \frac{\partial^n}{\partial\mu^n} \frac{\partial\Omega}{\partial\mu} \right)_{v(\mathbf{r})}
          = \left( \frac{\partial^n N}{\partial \mu^n} \right)_{v(\mathbf{r})}

More specifically, the fundamental global indicators of the :math:`N_0` electron reference state
are the first, second and higher order derivatives evaluated at :math:`N_0` which result in the number of electrons
denoted by :math:`N`, **chemical softness** denoted by :math:`S`, and :math:`n^{\text{th}}`
**-order hyper-softness** denoted by :math:`S^{(n)} \text{for } n \geq 2`, respectively:

 .. math::

    - \left. \left( \frac{\partial\Omega}{\partial\mu} \right)_{v(\mathbf{r})} \right|_{N = N_0} &= N \\
    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} \\
    S^{(2)} = - \left. \left( \frac{\partial^3\Omega}{\partial\mu^3} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 \\
    S^{(3)} = - \left. \left( \frac{\partial^4\Omega}{\partial\mu^4} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= \frac{-\sum_{k=1}^{n-1} S^k \cdot B_{n,k}
              \left(\eta^{(1)}, \eta^{(2)}, ..., \eta^{(n-k+1)} \right)}{B_{n,n}\left( \eta^{(1)}\right)}

The explicit formulas for softness and hyper-softness are obtained using chain rule;
please reafer to :ref:`derivation_global_softness` for details.

 .. TODO::
    #. Work on the derivation so explicit formula for hyper_softness :ref:`derivation_global_softness`

It is clear that the key quantity required for calculating various fundamental global reactivity
descriptors is the derivatives of the energy model with respect to the number of electrons :math:`N`
at fixed external potential :math:`v(\mathbf{r})`. In what follows, these derivatives are calculated
for various energy models.


.. _derived_indicators:

**Derived Global Reactivity Descriptors**

These reactivity indicators are derived based on some handwaving analysis,
or merely based on correlation. The most important one in the maximum number of electrons that can be
accepted by the system denoted by :math:`N_{\text{max}}`.

 .. math:: N_{\text{max}} &= \underbrace {\min }_N E(N)

**Electrophilicity index** :math:`\omega_{\text{electrophilicity}}` measures the capability
of an agent to accept electrons from the environment. However, in contrast to electron affinity :math:`EA`
which measures the energy loweing due to adding one electron to the system, electrophilicity index
:math:`\omega_{\text{electrophilicity}}` measures the energy lowering due to maximal electron flow
(which may be either less or more than one) from the environment,

 .. math:: \omega_{\text{electrophilicity}} = E(N_0) - E(N_0 + \Delta N_{\text{max}}) = E(N_0) - E(N_{\text{max}})

 .. math:: \omega_{\text{nucleophilicity}} = ?

 .. TODO::
    #. Talk about nucleophilicity; is it related to IP in the same way that electrophilicity is realted to EA.

**Nucleofugality index** :math:`\nu_{\text{nucleofugality}}` measures the susceptibility/quality/ability of the system for being
a nucleofuge (a leaving group which takes an electron with it) which is quantified by the energy penalty associated with forcing
a molecular fragment to accept an electron; the lower values of :math:`\nu_{\text{nucleofugality}}` are associated with high nucleofugality.

 .. math:: \nu_{\text{nucleofugality}} = E(N_0 + 1) - E(N_0 + \Delta N_{\text{max}}) = E(N_0 + 1) - E(N_{\text{max}})

**Electrofugality index** :math:`\nu_{\text{electrofugality}}` measures the susceptibility/quality/ability of the system for being
a electroguge (a leaving group which leaves an electron behind),

 .. math:: \nu_{\text{electrofugality}} = E(N_0 - 1) - E(N_0 + \Delta N_{\text{max}}) = E(N_0 - 1) - E(N_{\text{max}})

 .. TODO::
    #. Elaborate on derived tools
    #. Paul had sgn function used in his notes for defining these tools, is it required?

 .. math::

    {\Delta E}_{\text {electrophile}} &= \text {sgn}(N_0 - N_{max}) (E(N_0) - E(N_{max})) \\
    {\Delta E}_{\text {nucleophile}} &= ? \\
    {\Delta E}_{\text {nucleofuge}} &= \text {sgn}(N_0 + 1 - N_{max}) (E(N_0 + 1) - E(N_{max})) \\
    {\Delta E}_{\text {electrofuge}} &= \text {sgn}(N_0 - 1 - N_{max}) (E(N_0 - 1) - E(N_{max}))


.. _linear_energy:

Linear Energy Model :class:`chemtools.tool.globaltool.LinearGlobalTool`
=======================================================================

In this model, energy is approximated as a piece-wise linear function of the number of electrons:

 .. math:: E(N) = a + b N

 .. TODO::
    Technically, linear model has two parameters, so providing two E values is enough to fit the model.

The model requires three values of :math:`E(N)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given data points results in three equations:

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

In this model, second and higher order derivatives of energy with respect to the numbr of electrons is zero
(or not defined?).
So, chemical hardness and hyper-hardness are zero, and softness and hyper-softness are not defined.

 .. TODO::
    Is it better to skip derived global tools for this model?
    How the code should handel these?

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


.. _quadratic_energy:

Quadratic Energy Model :class:`chemtools.tool.globaltool.QuadraticGlobalTool`
=============================================================================

In this model, energy is approximated as a quadratic function of the number of electrons:

 .. TODO::
    #. Fix Equation number here, and assign number to other equations

 .. math::
    :nowrap:
    :label: quadratic

    \begin{eqnarray}
     E\left(N\right) = a + b N + c {N^2}
    \end{eqnarray}

Containing three parameters, :math:`a`, :math:`b` and :math:`c`, this model requires
three values of :math:`E\left(N\right)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations:

 .. math::

    \begin{cases}
          E\left(N_0 - 1\right) &= a + b \left(N_0 - 1\right) + c {\left(N_0 - 1\right) ^2} \\
             E \left(N_0\right) &= a + b \left(N_0\right) + c {\left(N_0\right) ^2} \\
          E\left(N_0 + 1\right) &= a + b \left(N_0 + 1\right) + c {\left(N_0 + 1\right) ^2}
    \end{cases}

This allows us to solve for the three unknowns:

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

    \Delta E = \left(\frac{IP + EA}{2}\right) \Delta N + \left(\frac{IP - EA}{2}\right) (\Delta N)^2

As detailed below, the prefactor of :math:`\Delta N` is the first derivative of energy with respect to :math:`N`
and the prefactor of :math:`(\Delta N)^2` is one-half the second order derivatives of energy with
respect to :math:`N` at fixed external potential
:math:`v(\mathbf{r})` evaluated at :math:`N = N_0`. As a result, this energy model is equivalent
to the second-order Taylor expansion of the energy as a function of :math:`N` around the reference
state :math:`N_0`.

 .. TODO::
    Check! The sings would not match the Tylor series...

To obtain the :ref:`fundamental global reactivity indicators <fundamental_indicators>` for the
quadratic energy model, the derivatives of the energy with respect to the number of electrons at
fixed external potential :math:`v(\mathbf{r})` should be calculated. These are given by:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}
         &= b + 2cN \\
	 &= \frac{E(N_0 + 1) - E(N_0 - 1)}{2} + \left(\frac{E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1)}{2}\right) \left(N - N_0\right) \\
	 &= -\frac{IP + EA}{2} + (IP - EA) \left(N - N_0\right) \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}
         &= 2c \\
	 &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) \\
	 &= IP - EA \\
    \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})}
         &= 0 \text{   for   } n \geq 2

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left(\frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0}
       &= \frac{E(N_0 + 1) - E(N_0 - 1)}{2}  = - \frac{{IP + EA}}{2} \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
        &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) = IP - EA \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
              &= 0 \text{   for   } n \geq 2

These are implemented in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_hardness`.

Accordingly, given the quadratic energy model, chemical softness and hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} = \frac{1}{IP - EA} \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= 0 \text {     for } n \geq 2

To obtain the :ref:`derived global reactivity indicators <derived_indicators>` for
the quadratic energy model, the maximum number of electrons to saturate the system should be calculated.
This is obtained by setting the first order derivative of energy, derived in Eq. ???, equal to zero:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} = 0 &= b + 2cN = -\frac{IP + EA}{2} + (IP - EA)(N - N_0) \\
    & \to N_{max} = \frac{-b}{2c} = N_{0} + \frac{IP + EA}{2 \left(IP - EA \right)} = N_{0} - \frac{\mu}{\eta} \\
    & \to \Delta N_{\text{max}} = N_0 - N_{\text{max}} = \frac{IP + EA}{2 \left(IP - EA \right)} = - \frac{\mu}{\eta}

The related :ref:`derived global reactivity indicators <derived_indicators>` for the quadratic energy model are:

 .. TODO::
    #. Show in more detail where these equations are coming from!!!
    #. Add **Electrodonating power** and **Electroaccepting power** for only quadratic model (in each interval)
    #. Add chemical potential defined by GV (quadratic model in each interval)

 .. math::

    \omega_{\text{electrophilicity}} &= E\left(N_0\right) - E\left(N_{\text max}\right)
                        &&= \frac{\left(IP + EA\right)^2}{8\left(IP - EA\right)}
		       &&&= \frac{\mu^2}{2 \eta} \\
    \omega_{\text{nucleophilicity}} &= ? \\
    \nu_{\text{nucleofugality}} &= E\left(N_0 + 1\right) - E\left(N_{\text max}\right)
                                &&= \frac{\left(IP - 3 \cdot EA \right)^2}{8 \left(IP - EA \right)}
			       &&&=  \frac{\left(\mu + \eta\right)^2}{2\eta} = -EA + \omega_{\text{electrophilicity}} \\
    \nu_{\text{electrofugality}} &= E\left(N_0 - 1\right) - E\left(N_{\text max}\right)
                                 &&= \frac{\left(3 \cdot IP - EA \right)^2}{8 \left(IP - EA \right)}
				&&&= \frac{\left(\mu - \eta\right)^2}{2\eta} = IP + \omega_{\text{electrophilicity}}

 .. TODO::
    #. Add references

**References:**
  * `Parr R. G., Pearson R. G., J. Am. Chem. Soc. (1983), 105, 7512 <http://pubs.acs.org/doi/abs/10.1021/ja00364a005>`_.

Sample Code:

 .. TODO::
    #. It would be nice to have the actual values showing up; something like IPython, or at least comment the results that should
       be obtained.

 .. code-block:: python
    :linenos:
    :emphasize-lines: 6

    import chemtools
    # H2O molecule with N0=10 electrons, & E(9)= , E(10)= , E(11)=
    energy_zero = 0.0  # E(N0) = E(10) =
    energy_plus = 0.0
    energy_minus = 0.0
    model = QuadraticGlobalTool(energy_zero, energy_plus, energy_minus, 10)
    # Retrieve global descriptors
    print model.chemical_potential
    print model.mu
    print model.chemical_hardness
    print model.eta
    print model.softness
    print model.hyper_hardness(2)


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

Sample Code:

 .. TODO::
    #. Add sample code!


.. _rational_energy:

Rational Energy Model :class:`chemtools.tool.globaltool.RationalGlobalTool`
===========================================================================

In this model, energy is approximated by a rational function of the number of electrons.
In the most general form, this model can be written as:

 .. math::

    E^{(m,n)}\left(N\right) = \left( \frac{a_0 + a_1N + a_2{N^2} + ... + a_m{N^m}}{1 + b_1N + b_2{N^2} + ... + b_n{N^n}} \right)
                 = \frac{\sum_{j=0}^{m} a_j N^j}{1 + \sum_{i=1}^{n} b_i N^i}

The number of unknown parameters in this model depends on the :math:`m` and :math:`n` values.
Having a set of :math:`m+n` values of :math:`N` for which the energy is known, the model can be parametrized
by solving a system of linear equations. By rearranging the rational energy expression above,
the equations can be written as:

 .. math::

    \sum_{j=0}^{m} \left(N^j\right) a_j - \sum_{i=1}^{n} \left(N^i \cdot E^{(m,n)}\left(N\right) \right) b_i = E^{(m,n)}\left(N\right)

Having the parameters :math:`\{a_j\}_{j=0}^m` and :math:`\{b_i\}_{i=1}^n`, the energy model is known,
and the derivatives of the rational energy model with respect to the number of electrons at fixed external
potential can be calculated.

However, in order to solve for the parameters in this model analytically, a simpler form of the rational energy model
containing three parameters, :math:`E^{(2,1)}\left(N\right) = E\left(N\right)`, is considered. For implementing more
complex rational energy models, please refer to the :ref:`general energy model <general_energy>`.

 .. math:: E\left(N\right) = E^{(2,1)}\left(N\right) = \frac{a_0 + a_1 N}{1 + b_1 N}

Containing three parameters, :math:`a_0`, :math:`a_1` and :math:`b_1`, this model requires
three values of :math:`E\left(N\right)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations:

 .. math::

    \begin{cases}
     \left(1 + b_1 \left(N_0 - 1\right)\right) & E\left(N_0-1\right) &&= a_0 + a_1 \left(N_0 - 1\right)  \\
     \left(1 + b_1 N_0\right) & E\left(N_0\right) &&= a_0 + a_1 N_0 \\
     \left(1 + b_1 \left(N_0 + 1\right)\right) & E\left(N_0+1\right) &&= a_0 + a_1 \left(N_0 + 1\right) \\
    \end{cases}

This allows us to solve for the three unknonws:

 .. math::

    b_1 &= -\frac{E\left(N_0 + 1\right) - 2 E\left(N_0\right) + E\left(N_0 - 1\right)}
                 {\left(N_0 + 1\right) E\left(N_0 + 1\right) - 2 N_0 E\left(N_0\right) + \left(N_0 - 1\right) E\left(N_0 - 1\right)} \\
    a_1 &= \left(1 + b_1 N_0\right) \left(E\left(N_0 + 1\right) - E\left(N_0\right)\right) + b_1 E\left(N_0 + 1\right) \\
    a_0 &= - a_1 N_0 + \left(1 + b_1 N_0\right) E\left(N_0\right)

Due to the complexity of the obtained parameters, we skip substituting them into the energy expression.
However, at this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.RationalGlobalTool.energy`.

The derivatives of the energy model with respect to the number of electrons at
fixed external potential are:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}
	 &= \frac{a_1 - a_0 b_1}{\left(1 + b_1 N\right)^2} \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}
         &= \frac{-2 b_1 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N\right)^3} \\
    \left( \frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})}
         &= \frac{(-b_1)^{n - 1} \left(a_1 - a_0 b_1\right) n!}{\left(1 + b_1 N\right)^{n+1}}

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.RationalGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0}
       &= \frac{a_1 - a_0 b_1}{\left(1 + b_1 N_0\right)^2} \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
        &= \frac{-2 b_1 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^3} \\
    \eta^{(2)} = \left. \left( \frac{\partial^3 E}{\partial N^3} \right)_{v(\mathbf{r})} \right|_{N = N_0}
         &= \frac{6 b_1^2 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^4} \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
         &= \frac{(-b_1)^n \left(a_1 - a_0 b_1\right) \left(n+1\right)!}{\left(1 + b_1 N_0\right)^{n+2}} \text{   for } n\geq2

These are implemented in :class:`chemtools.tool.globaltool.RationalGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.RationalGlobalTool.chemical_hardness`.

Accordingly, given the rational energy model, chemical softness and hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} = \frac{\left(1 + b_1 N_0\right)^3}{2 b_1 \left(a_1 - a_0 b_1\right)} \\
    S^{(2)} = - \left. \left( \frac{\partial^{3}\Omega}{\partial\mu^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 \\
           &= -\frac{6 b_1^2 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^4}
	       \frac{\left(1 + b_1 N_0\right)^9}{2^3 b_1^3 \left(a_1 - a_0 b_1\right)^3}
           = \frac{-3 \left(1 + b_1 N_0\right)^5}{2 b_1 \left(a_1 - a_0 b_1\right)^2} \\
    S^{(3)} = - \left. \left( \frac{\partial^{4}\Omega}{\partial\mu^{4}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
	   &= -\frac{24 b_1^3 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^5}
	       \frac{\left(1 + b_1 N_0\right)^12}{2^4 b_1^4 \left(a_1 - a_0 b_1\right)^4} \\
	   &  + 3\frac{6^2 b_1^4 \left(a_1 - a_0 b_1\right)^2}{\left(1 + b_1 N_0\right)^8}
	       \frac{\left(1 + b_1 N_0\right)^15}{2^5 b_1^5 \left(a_1 - a_0 b_1\right)^5} \\
	   &= \frac{15 \left(1 + b_1 N_0\right)^7}{8 b_1 \left(a_1 - a_0 b_1\right)^3}

The higher order hyper-softness exists and can be evaluated through Eq. ???, as implemented in
:meth:`chemtools.tool.globaltool.RationalGlobalTool.hyper_softness`.

To obtain the :ref:`derived global reactivity indicators <derived_indicators>` for
the exponential energy model, the maximum number of electrons accepted by the system should be calculated.

 .. TODO::
    #. Write down the value of N_max and derived global reactivity tools

**References:**

 .. TODO::
    #. Add references

Sample Code:

 .. TODO::
    #. Add sample code!


.. _general_energy:

General Energy Model :class:`chemtools.tool.globaltool.GeneralGlobalTool`
=========================================================================

In this model, energy is approximated by an user-specified energy model. Given the
known energy values, this model is parametrized and the energy expression can be evaluated
for any number of electrons.
Being a generic models, this model can reproduce the results of
:ref:`linear <linear_energy>`, :ref:`quadratic <quadratic_energy>`, :ref:`exponential <exponential_energy>`,
and :ref:`rational <rational_energy>` energy models as special cases.

The energy expression should be specified symbolically through `Sympy <http://www.sympy.org/en/index.html>`_.


 .. TODO::
    12. Elaborate more on this model.
    #. Add sample code!

Example: Build a quadratic energy model:

  .. ipython:: python

     import chemtools
     import sympy

  .. ipython:: python

     # define symbols used in the energy expression
     n, a, b, c = sympy.symbols('N, a, b, c')
     # define the energy expression
     expression = a + b * n + c * (n**2)
     # dictionary {N : E(N)}
     energies = {}
     # parametrize energy model
     model = GeneralizedGlobalTool(expression, energies, n)
     # ready to retrieve any global tool
     print model.mu


Analytical
==========

Here the analytical evaluation of chemical potential and hardness, etc. will be discussed!
