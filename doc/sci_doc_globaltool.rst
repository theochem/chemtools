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

 .. math:: N_{\text{max}} = \underbrace {\min }_N E(N)

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

In the following sections, various energy models are considered and the corresponding global reactivity descriptors are obtained.

.. toctree::
   :maxdepth: 2

   sci_doc_global_linear
   sci_doc_global_quadratic
   sci_doc_global_rational
   sci_doc_global_exponential
   sci_doc_global_general

