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


.. _rational_energy:

Rational Energy Model :class:`chemtools.tool.globaltool.RationalGlobalTool`
===========================================================================

In this model, energy is approximated by a rational function of the number of electrons.
In the most general form, this model can be written as:

 .. math::

    E^{(m,n)}\left(N\right) = \left( \frac{a_0 + a_1N + a_2{N^2} + ... + a_m{N^m}}{1 + b_1N + b_2{N^2} + ... + b_n{N^n}} \right)
                 = \frac{\sum_{j=0}^{m} a_j N^j}{1 + \sum_{i=1}^{n} b_i N^i}

The number of unknown parameters in this model depends on the values of :math:`m` and :math:`n`.
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
three values of :math:`E\left(N\right)` to interpolate the energy. Commonly, the energy of the system
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

 .. math::

    a_0 &=  \frac{E\left(N_0\right) E\left(N_0-1\right) N_{0} + E\left(N_0\right) E\left(N_0-1\right) + E\left(N_0\right) E\left(N_0+1\right) N_{0} -
                E\left(N_0\right) E\left(N_0+1\right) - 2 E\left(N_0-1\right) E\left(N_0+1\right) N_{0}}{2 E\left(N_0\right) N_{0} - E\left(N_0-1\right) N_{0} + E\left(N_0-1\right) - E\left(N_0+1\right) N_{0} - E\left(N_0+1\right)} \\
    a_1 &=  \frac{- E\left(N_0\right) E\left(N_0-1\right) - E\left(N_0\right) E\left(N_0+1\right) + 2 E\left(N_0-1\right) E\left(N_0+1\right)}{2 E\left(N_0\right) N_{0} - E\left(N_0-1\right) N_{0} + E\left(N_0-1\right) - E\left(N_0+1\right) N_{0} - E\left(N_0+1\right)} \\
    b_1 &=  \frac{- 2 E\left(N_0\right) + E\left(N_0-1\right) + E\left(N_0+1\right)}{2 E\left(N_0\right) N_{0} - E\left(N_0-1\right) N_{0} + E\left(N_0-1\right) - E\left(N_0+1\right) N_{0} - E\left(N_0+1\right)}

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

In the 3-point rational model, evaluating the first-, second-, and higher-order derivatives of energy evaluated
at :math:`N_0` gives the following expressions for the chemical potential, chemical hardness, and hyper-hardnesses,

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

Using these expressions, one can derive the following expressions for the chemical softness and the low-order
hyper-softnesses,

 .. math::

    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} = \frac{-\left(1 + b_1 N_0\right)^3}{2 b_1 \left(a_1 - a_0 b_1\right)} \\
    S^{(2)} = - \left. \left( \frac{\partial^{3}\Omega}{\partial\mu^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 \\
           &= -\frac{6 b_1^2 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^4}
	       \frac{\left(1 + b_1 N_0\right)^9}{2^3 b_1^3 \left(a_1 - a_0 b_1\right)^3}
           = \frac{3 \left(1 + b_1 N_0\right)^5}{4 b_1 \left(a_1 - a_0 b_1\right)^2} \\
    S^{(3)} = - \left. \left( \frac{\partial^{4}\Omega}{\partial\mu^{4}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
	   &= -\frac{24 b_1^3 \left(a_1 - a_0 b_1\right)}{\left(1 + b_1 N_0\right)^5}
	       \frac{\left(1 + b_1 N_0\right)^12}{2^4 b_1^4 \left(a_1 - a_0 b_1\right)^4} \\
	   &  + 3\frac{6^2 b_1^4 \left(a_1 - a_0 b_1\right)^2}{\left(1 + b_1 N_0\right)^8}
	       \frac{\left(1 + b_1 N_0\right)^15}{2^5 b_1^5 \left(a_1 - a_0 b_1\right)^5} \\
	   &= \frac{-15 \left(1 + b_1 N_0\right)^7}{8 b_1 \left(a_1 - a_0 b_1\right)^3}


ChemTools can also compute higher-order hyper-softnesses, using the (extended) inverse function theorem for
derivatives. Please refer to :ref:`derivation_global_softness` for details.

To obtain the :ref:`derived global reactivity indicators <global_derived_indicators>` for
the rational energy model, the maximum number of electrons accepted by the system should be calculated.

 .. TODO::
    #. Incude :math:`N_{\text{max}}=\infty` and derived global reactivity tools

**References:**

 .. TODO::
    #. Add references
