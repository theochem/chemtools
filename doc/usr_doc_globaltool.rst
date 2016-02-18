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


Global Descriptive Tools :mod:`chemtools.tool.globaltool`
#########################################################

Global descriptive tools assign a single value to the entire molecule.
The global reactivity indicators from conceptual DFT are either :ref:`fundamental <fundamental_indicators>`
or :ref:`derived <derived_indicators>`. As explained below, obtaining these global tools requires selecting
an energy model, :math:`E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})`,  which represents the dependence of energy on the
number electrons, :math:`N`. The set :math:`{\{\alpha_1, \alpha_2, ..., \alpha_n\}}` denotes the parameters
of the energy model which are determined by fitting the energy expression to the known values of the energy for :math:`n` different numbers of electrons,
:math:`{\{E(N_i)\}}_{i=1}^n`. Commonly, the values of energy for systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` are used, :math:`{\{E(N_0 - 1), E(N_0), E(N_0 + 1)\}}`; however, other values of energy can be used to parametrized the model as well. Needless to say, the number of required energy values to solve for the parameters depends on the complexity of the energy model, i.e. the number of parameters in the model.
The implemented energy models and the corresponding global reactivity descriptors in ChemTools include:

 #. :ref:`Linear Energy Model <linear_energy>`
 #. :ref:`Quadratic Energy Model <quadratic_energy>`
 #. :ref:`Exponential Energy Model <exponential_energy>`
 #. :ref:`Rational Energy Model <rational_energy>`
 #. :ref:`General Energy Model <general_energy>`

.. _fundamental_indicators:

**Fundamental Global Reactivity Descriptors:**
In the canonical ensemble, these include the derivatives of the energy model with respect to the number of electrons :math:`N` at fixed external potential :math:`v(r)` :

 .. math:: \left( \frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}{\partial N} \right)_{v(r)}

More specifically, these fundamental global indicators are the first, second and higher order derivatives
evaluated at :math:`N=N_{0}` which are called **chemical potential** denoted by :math:`\mu`,
**chemical hardness** denoted by :math:`\eta`, and :math:`n^{(\text th)}` **-order hyper-hardness**
denoted by :math:`\eta^{(n)} \text{for } n \geq 2`, respectively:

 .. math::

    \mu \equiv \left. \left( \frac{\partial E}{\partial N} \right)_{v(r)} \right|_{N = N_0} & \\
    \eta \equiv \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)} \right|_{N = N_0} & \\
    \eta^{(n)} \equiv \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(r)} \right|_{N = N_0} & \text{for } n \geq 2

In the grand canonical ensemble, these include the derivatives of the grand
potential model :math:`\Omega = E (\left\langle N \right\rangle) - \mu \left\langle N \right\rangle`
with respect to the chemical potential :math:`\mu` at fixed external potential :math:`v(r)`.

 .. TODO::
    #. Elaborate on grand potential

 .. math::

    - \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)}
          = - \left( \frac{\partial^n}{\partial\mu^n} \frac{\partial\Omega}{\partial\mu} \right)_{v(r)}
          = \left( \frac{\partial^n N}{\partial \mu^n} \right)_{v(r)}

More specifically, these fundamental global indicators are the first, second and
higher order derivatives evaluated at :math:`N_0` which result in the number of electrons
denoted by :math:`N`, **chemical softness** denoted by :math:`S`, and :math:`n^{(\text th)}`
**-order hyper-softness**, respectively:

 .. math::

    - \left. \left( \frac{\partial\Omega}{\partial\mu} \right)_{v(r)} \right|_{N = N_0} &= N \\
    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(r)} \right|_{N = N_0}
     &= \frac{1}{\eta} \\
    S^{(2)} = - \left. \left( \frac{\partial^3\Omega}{\partial\mu^3} \right)_{v(r)} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 \\
    S^{(3)} = - \left. \left( \frac{\partial^4\Omega}{\partial\mu^4} \right)_{v(r)} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)} \right|_{N = N_0}
           &= \text {complicated formula for Faa di Bruno identity}

Using the Faa di Bruno formula, :ref:`hyper-softness can be solved explicitly <derivation_softness>`
using the recursive formula,

 .. TODO::
    #. Work on the derivation so explicit formula for hyper_softness :ref:`derivation_softness`

 .. math::

    S^{(n)} = \frac{-\sum_{k=1}^{n-1} S^k \cdot B_{n,k}
          \left(\eta^{(1)}, \eta^{(2)}, ..., \eta^{(n-k+1)} \right)}{B_{n,n}\left( \eta^{(1)}\right)}

It is clear that fundamental global descriptors based in the grand canonical ensemble can be calculated
based on the fundamental global descriptors in the canonical ensemble.

.. _derived_indicators:

**Derived Global Reactivity Descriptors:**
These reactivity indicators are derived based on some handwaving analysis,
or merely based on correlation. The most important one is the maximum number
of electrons that can be accepted and the related energetic quantities like
**electrophilicity**, **nucleofugality**, and **electrofugality**.

 .. TODO::
    #. Elaborate on derived tools

 .. math::

    N_{max} &= \underbrace {\min }_N E(N) \\
    \omega_{\text {electrophile}} &= \text {sgn}(N_0 - N_{max}) (E(N_0) - E(N_{max})) \\
    \omega_{\text {nucleophile}} &= ? \\
    \nu_{\text {nucleofuge}} &= \text {sgn}(N_0 + 1 - N_{max}) (E(N_0 + 1) - E(N_{max})) \\
    \nu_{\text {electrofuge}} &= \text {sgn}(N_0 - 1 - N_{max}) (E(N_0 - 1) - E(N_{max}))


.. _linear_energy:

Linear Energy Model: :class:`chemtools.tool.globaltool.LinearGlobalTool`
========================================================================

In this model, energy is approximated as a piece-wise linear function of the number of electrons:

 .. math:: E(N) = a + b N

The model requires three values of :math:`E(N)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given data points results in three equations:

 .. math::

    E(N) &= \begin{cases}
             (N_0 - N) E(N_0 - 1) + (N - (N_0 - 1)) E(N_0) =
	     E(N_0) + (N_0 - N) \cdot IP & \text{  for  } N < N_0 \\
	     (N_0 + 1 + N) E(N_0 - 1) + (N - N_0) E(N_0 + 1) =
	     E(N_0) + (N_0 - N) \cdot EA & \text{  for  } N \geqslant N_0 \\
	    \end{cases} \\

At this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.LinearGlobalTool.energy`.

Because the energy model is not differentiable at integer number of electrons, the chemical potential
is not defined and is instead calculated from above, below and averaged:

 .. math::

    \mu^{-} &= -I \\
    \mu^{0} &= \frac{\mu^{+} + \mu^{-}}{2} \\
    \mu^{+} &= -A \\

 .. todo::

    This still needs some work!


.. _quadratic_energy:

Quadratic Energy Model: :class:`chemtools.tool.globaltool.QuadraticGlobalTool`
==============================================================================

In this model, energy is approximated as a quadratic function of the number of electrons:

 .. TODO::
    #. Fix Equation number here, and assign number to other equations

 .. math::
    :nowrap:
    :label: quadratic

    \begin{eqnarray}
     E(N) = a + b N + c {N^2}
    \end{eqnarray}

Containing three parameters, :math:`a`, :math:`b` and :math:`c`, this model requires
three values of :math:`E(N)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations:

 .. math::

    \begin{cases}
          E(N_0 - 1) &= a + b (N_0 - 1) + c {(N_0 - 1) ^2} \\
             E (N_0) &= a + b (N_0) + c {(N_0) ^2} \\
          E(N_0 + 1) &= a + b (N_0 + 1) + c {(N_0 + 1) ^2}
    \end{cases}

This allows us to solve for the three unknowns:

 .. math::

    a &= E(N_0) - b N_0 - c {N_0 ^2} \\
    b &= \frac{E(N_0 + 1) - E(N_0 - 1)}{2} - 2 N_0 c \\
    c &= \frac{E(N_0 - 1) -2 E(N_0) + E(N_0 + 1)}{2} \\

Substituting the obtained parameters :math:`a`, :math:`b` and :math:`c` into the energy expression,
Eq. :eq:`quadratic`, gives the fitted energy model as:

 .. math::

          E(N) &=& E(N_0) + (N - N_0) \left( \frac{E(N_0 + 1) - E(N_0 - 1)}{2} \right)  \\
                       & & +    {(N - N_0) ^2} \left( \frac{E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1)}{2} \right) \\
	       &=& E(N_0) - (N - N_0) \left( \frac{IP + EA}{2} \right) + {(N - N_0) ^2} \left( \frac{IP - EA}{2} \right) \\

where :math:`IP = E(N_0 - 1) - E(N_0)` and :math:`EA = E(N_0) - E(N_0 + 1)` are the
ionization potential and electron affinity of the :math:`N_0` -electron system, respectively.
At this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.energy`.

To obtain the :ref:`fundamental global reactivity indicators <fundamental_indicators>` for the
quadratic energy model, the derivatives of the energy with respect to the number of electrons at
fixed external potential should be calculated. These are given by:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(r)}
         &= b + 2cN \\
	 &= \left(\frac{E(N_0 + 1) - E(N_0 - 1)}{2}\right) + (N - N_0) \left(\frac{E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1)}{2}\right) \\
	 &= -\left( \frac{IP + EA}{2} \right) + (N - N_0) (IP - EA) \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)}
         &= 2c \\
	 &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) \\
	 &= IP - EA \\
    \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(r)}
         &= 0 \text{   for   } n \geq 2

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left( \frac{\partial E}{\partial N} \right)_{v(r)} \right|_{N = N_0}
       &= \left(\frac{E(N_0 + 1) - E(N_0 - 1)}{2}\right)  = - \left(\frac{{IP + EA}}{2}\right) \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)} \right|_{N = N_0}
        &= E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1) = IP - EA \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(r)} \right|_{N = N_0}
              &= 0 \text{   for   } n \geq 2

These are implemented in :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.QuadraticGlobalTool.chemical_hardness`.

Accordingly, given the quadratic energy model, chemical softness and hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(r)} \right|_{N = N_0}
     &= \frac{1}{\eta} = \frac{1}{IP - EA} \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)} \right|_{N = N_0}
           &= 0 \text {     for } n \geq 2

To obtain the :ref:`derived global reactivity indicators <derived_indicators>` for
the quadratic energy model, the maximum number of electrons accepted by the system should be calculated.
This is obtained by setting the first order derivative of energy, derived in Eq. ???, equal to zero:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(r)} = 0 &= b + 2cN = -\left( \frac{IP + EA}{2} \right) + (N - N_0) (IP - EA) \\
    & \to N_{max} = \frac{-b}{2c} = N_{0} + \frac{IP + EA}{2 \left(IP - EA \right)} = N_{0} - \frac{\mu}{\eta}

The related derived global reactivity indicators for the quadratic energy model are:

 .. TODO::
    #. Show in more detail where these equations are coming from!!!

 .. math::

    \omega_{\text {electrophile}} &= \frac{\mu^2}{2 \cdot \eta} \\
    \omega_{\text {nucleophile}} &= ? \\
    \nu_{\text {nucleofuge}} &= \frac{(IP - 3 \cdot A)^2}{8 (IP - EA)} \\
    \nu_{\text {electrofuge}} &= \frac{(3 \cdot IP - A)^2}{8 (IP - EA)}

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

Exponential Energy Model: :class:`chemtools.tool.globaltool.ExponentialGlobalTool`
==================================================================================

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

    A      &= \frac{(E(N_0 - 1) - E(N_0))(E(N_0) - E(N_0 + 1))}{E(N_0 - 1) - 2 E(N_0) + E(N_0 + 1)}
            = \frac{IP \cdot EA}{IP - EA} \\
    B      &= E(N_0) - A  \\
    \gamma &= \ln \left( 1 - \frac{E(N_0 - 1) - 2E(N_0) + E(N_0 + 1)}{E(N_0 + 1) - E(N_0)} \right) \\

Due to the complexity of the obtained parameters, we skip substituting them into the energy expression.
However, at this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.energy`.

The derivatives of the energy model with respect to the number of electrons at
fixed external potential are:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(r)}
         &= A (-\gamma) \exp(-\gamma (N - N_0)) \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)}
         &= A {(-\gamma) ^2} \exp(-\gamma (N - N_0)) \\
    \left( \frac{\partial^n E}{\partial N^n} \right)_{v(r)}
         &= A {(-\gamma) ^n} \exp(-\gamma (N - N_0)) \text{   for   } n \geq 1

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left( \frac{\partial E}{\partial N} \right)_{v(r)} \right|_{N = N_0}
       &= -A \gamma \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)} \right|_{N = N_0}
        &= A {\gamma ^2} \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n} E}{\partial N^{n}} \right)_{v(r)} \right|_{N = N_0}
              &= A {(- \gamma) ^n} \text{   for   } n \geq 1

These are implemented in :class:`chemtools.tool.globaltool.ExponentialGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.ExponentialGlobalTool.chemical_hardness`.

Accordingly, given the exponential energy model, chemical softness and hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(r)} \right|_{N = N_0}
     &=  \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)} \right|_{N = N_0}
           &=  \text {     for } n \geq 2

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

    E^{(m,n)}(N) = \left( \frac{a_0 + a_1N + a_2{N^2} + ... + a_m{N^m}}{1 + b_1N + b_2{N^2} + ... + b_n{N^n}} \right)
                 = \frac{\sum_{j=0}^{m} a_j N^j}{1 + \sum_{i=1}^{n} b_i N^i}

The number of unknown parameters in this model depends on the :math:`m` and :math:`n` values.
Having a set of :math:`m+n` values of :math:`N` for which the energy is known, the model can be parametrized
by solving a system of linear equations. By rearranging the rational energy expression above,
the equations can be written as:

 .. math::

    \sum_{j=0}^{m} (N^j) a_j - \sum_{i=1}^{n} \left(N^i \cdot E^{(m,n)}(N) \right) b_i = E^{(m,n)}(N)

Having the parameters :math:`\{a_j\}_{j=0}^m` and :math:`\{b_i\}_{i=1}^n`, the energy model is known,
and the derivatives of the rational energy model with respect to the number of electrons at fixed external
potential can be calculated.

However, in order to solve for the parameters in this model analytically, a simpler form of the rational energy model
containing three parameters, :math:`E^{(2,1)}(N) = E(N)`, is considered. For implementing more
complex rational energy models, please refer to the :ref:`general energy model <general_energy>`.

 .. math:: E(N) = E^{(2,1)}(N) = \frac{a_0 + a_1 N}{1 + b_1 N}

Containing three parameters, :math:`a_0`, :math:`a_1` and :math:`b_1`, this model requires
three values of :math:`E(N)` to interpolate energy. Commonly, the energy of the system
with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons are provided.
Fitting the energy expression to the given energy values results in three equations:

 .. math::

    \begin{cases}
     (1 + b_1 (N_0 - 1)) E(N_0-1) &= a_0 + a_1 (N_0 - 1)  \\
     (1 + b_1 N_0) E(N_0-1) &= (a_0 + a_1 N_0) \\
     (1 + b_1 (N_0 + 1)) E(N_0-1) &= (a_0 + a_1 (N_0 + 1)) \\
    \end{cases}

This allows us to solve for the three unknonws:

 .. math::

    a_0 &= \frac{N_0 E(N_0) - 2 N_0^2 E(N_0 + 1)}{N_0 + 1} \\
    a_1 &= \frac{2 N_0 E(N_0 + 1) + (2 N_0 - 1) E(N_0)}{N_0 + 1} \\
    b_1 &= - \frac{1}{N_0 + 1}

Due to the complexity of the obtained parameters, we skip substituting them into the energy expression.
However, at this stage, the energy expression can be evaluated for any given number of electrons as
implemented in :class:`chemtools.tool.globaltool.RationalGlobalTool.energy`.

The derivatives of the energy model with respect to the number of electrons at
fixed external potential are:

 .. math::

    \left( \frac{\partial E}{\partial N} \right)_{v(r)}
	 &= \frac{a_1 - a_0 b_1}{(1 + b_1 N)^2} \\
    \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)}
         &= 2 \frac{- b_1 (a_1 - a_0 b_1)}{(1 + b_1 N)^3} \\
    \left( \frac{\partial^n E}{\partial N^n} \right)_{v(r)}
         &= \frac{b_1^{n - 1} (a_1 - a_0 b_1) n!}{(1 + b_1 N)^{2n}}

These derivatives can be evaluated for any number of electrons as implemented
in :class:`chemtools.tool.globaltool.RationalGlobalTool.energy_derivative`.
In this model, the first, second and higher order derivatives of energy evaluated at :math:`N_0`,
the so-called chemical potential and chemical hardness and hyper-hardness, equal:

 .. math::

    \mu = \left. \left( \frac{\partial E}{\partial N} \right)_{v(r)} \right|_{N = N_0}
       &= \frac{a_1 - a_0 b_1}{(1 + b_1 N_0)^2} \\
    \eta = \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(r)} \right|_{N = N_0}
        &= 2 \frac{- b_1 (a_1 - a_0 b_1)}{(1 + b_1 N_0)^3} \\
    \eta^{(n)} = \left. \left( \frac{\partial^{n+1} E}{\partial N^{n+1}} \right)_{v(r)} \right|_{N = N_0}
        &= \frac{b_1^{n - 1} (a_1 - a_0 b_1)}{(1 + b_1 N_0)^2n} \text{   for   } n \geq 0

These are implemented in :class:`chemtools.tool.globaltool.RationalGlobalTool.chemical_potential`
and :class:`chemtools.tool.globaltool.RationalGlobalTool.chemical_hardness`.

Accordingly, given the rational energy model, chemical softness and hyper-softness equal:

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(r)} \right|_{N = N_0}
     &=  \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(r)} \right|_{N = N_0}
           &=  \text {     for } n \geq 2

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
    #. Elaborate more on this model.
    #. Add sample code!

Example: Build a quadratic energy model:

  .. code-block:: python
     :linenos:

     import chemtools
     import sympy

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
