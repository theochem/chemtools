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

The simplest global reactivity descriptors are the first **ionization potential**
, :math:`IP`, and **electron affinity**, :math:`EA`, of the system, measuring the system’s
propensity to donate or accept one electron. For a system with :math:`N_0` electrons,
these are defined as,

 .. math::
    IP = E\left(N_0 - 1\right) - E\left(N_0\right) \\
    EA = E\left(N_0\right) - E\left(N_0 + 1\right)

The global reactivity indicators from conceptual DFT can be classified as either
:ref:`fundamental <global_fundamental_indicators>` or :ref:`derived <global_derived_indicators>`.
Before one can evaluate any of these indicators, one must select
an energy model, :math:`E_{\text{model}} = E\left(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}\right)`,
representing the dependence of the energy on the
number of electrons :math:`N`. The :math:`{\{\alpha_1, \alpha_2, ..., \alpha_n\}}` set denotes the parameters
of the energy model which are determined by interpolating the energy expression to the known values of
the energy for :math:`n` different numbers of electrons,
:math:`{\{E(N_i)\}}_{i=1}^n`. Commonly, the system with :math:`N_0` electrons denotes the reference state,
and the energy values for systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons,
i.e. :math:`{\{E(N_0 - 1), E(N_0), E(N_0 + 1)\}}`, are used to find the parameters in the energy model.
However, other values of energy can be used to parametrized the model as well. Needless to say, the number of
energy values required to solve for the parameters depends on the complexity of the energy model,
i.e. the number of parameters in the model.

For global reactivity tools, ChemTools support the built-in and user-defined energy models listed below.
Using the code of these models as templates, users can develop new energy models in ChemTools.

 #. :ref:`Linear Energy Model <linear_energy>`
 #. :ref:`Quadratic Energy Model <quadratic_energy>`
 #. :ref:`Exponential Energy Model <exponential_energy>`
 #. :ref:`Rational Energy Model <rational_energy>`
 #. :ref:`General Energy Model <general_energy>`


.. _global_fundamental_indicators:

**Fundamental Global Reactivity Descriptors**

.. _energy_derivatives:

In the canonical ensemble, the fundamental global reactivity descriptors are the derivatives
of the energy model with respect to the number of electrons :math:`N` at fixed external potential
:math:`v(\mathbf{r})`:

 .. math:: \left( \frac{\partial^n E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}{\partial N^n} \right)_{v(\mathbf{r})} \qquad n=1,2,\dots

More specifically, the fundamental global reactivity indicators of the :math:`N_0-` electron reference state
are the first-, second- and higher-order derivatives of energy
evaluated at :math:`N=N_{0}`, which are called **chemical potential** denoted by :math:`\mu`,
**chemical hardness** denoted by :math:`\eta`, and :math:`n^{\text{th}}` **-order hyper-hardness**
denoted by :math:`\eta ^{(n)} \text{for } n \geq 2`, respectively:

 .. math::

    \mu &= \left. \left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \\
    \eta &= \left. \left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \\
    \eta^{(2)} &= \left. \left( \frac{\partial^{3} E}
                {\partial N^{3}} \right)_{v(\mathbf{r})} \right|_{N = N_0} \\
    \eta^{(n)} &= \left. \left( \frac{\partial^{n+1} E}
                {\partial N^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0} & \qquad \text{for } n \geq 2

By analogy to thermodynamics, the electronic chemical potential, :math:`\mu`, is identified as the first derivative
of the energy with respect to the number of electrons. The electronic chemical potential is always negative and it
represents the intrinsic Lewis acidity/basicity of a system.
Systems with very negative chemical potentials are strong Lewis acids; systems with nearly zero (but still negative)
chemical potential are strong Lewis bases. Electrons flow from places with higher chemical potential to places with
lower chemical potential until the chemical potential is everywhere equal.
This motivates the definition of electronegativity as the negative of the chemical potential, :math:`\chi_e=-\mu`.

The second derivative of the energy with respect to electron number is the chemical hardness, :math:`\eta`.
The chemical hardness is always positive. Hard molecules have properties (e.g., the chemical potential and electron density)
that are relatively insensitive to perturbations. In addition, according to the principle of maximum hardness,
it is frequently (but far from universally) the case that hard systems are more stable than soft systems.
According to the hard/soft acid/base (HSAB) principle: if all other effects are (nearly) equal, hard acids will prefer
to react with hard bases and soft acids will prefer to react with soft bases.

Higher-order derivatives of the energy with respect to electron number are called hyper-hardnesses, and are not often used.
The third derivative, :math:`\eta^{(2)}`, which is often simply called the hyper-hardness, appears in some contexts,
but is typically small.

ChemTools will evaluate energy derivatives to arbitrary order if those derivatives exist in the energy model you have selected.

Reactivity indicators in the canonical ensemble are appropriate when the number of electrons in the system is clearly defined.
For reactions in, or on the surface of, liquid or solid solutions, the number of electrons in the reacting molecule is determined
by the electronic chemical potential of its environment. In these cases, and also when comparing the reactivity of molecules of
disparate size, it is better to use reactivity descriptors from the grand canonical ensemble.


.. _grand_potential_derivatives:

In the grand canonical ensemble, the fundamental global reactivity descriptors are the derivatives of the grand
potential model :math:`\Omega = E \left(\left\langle N \right\rangle\right) - \mu \left\langle N \right\rangle`
with respect to the chemical potential :math:`\mu` at fixed external potential :math:`v(\mathbf{r})`.

 .. math::

    - \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})}
          = - \left( \frac{\partial^n}{\partial\mu^n} \frac{\partial\Omega}{\partial\mu} \right)_{v(\mathbf{r})}
          = \left( \frac{\partial^n N}{\partial \mu^n} \right)_{v(\mathbf{r})}  \qquad n=0,1,\dots

More specifically, the fundamental global indicators of the :math:`N_0-` electron reference state
are the first-, second- and higher-order derivatives of :math:`\Omega` with respect to :math:`\mu`
evaluated at the chemical potential of the :math:`N_0-` electron system, :math:`\mu_0=\mu\left(N_0\right)`.
These include the number of electrons denoted by :math:`N_0`, **chemical softness** denoted by :math:`S`, and :math:`n^{\text{th}}`
**-order hyper-softness** denoted by :math:`S^{(n)} \text{for } n \geq 2`, respectively:

 .. math::

    - \left. \left( \frac{\partial\Omega}{\partial\mu} \right)_{v(\mathbf{r})} \right|_{N = N_0} &= N_0 \\
    S = - \left. \left( \frac{\partial^2\Omega}{\partial\mu^2} \right)_{v(\mathbf{r})} \right|_{N = N_0}
     &= \frac{1}{\eta} \\
    S^{(2)} = - \left. \left( \frac{\partial^3\Omega}{\partial\mu^3} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(2)} \cdot S^3 \\
    S^{(3)} = - \left. \left( \frac{\partial^4\Omega}{\partial\mu^4} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= -\eta^{(3)} \cdot S^4 + 3 \left(\eta^{(2)}\right)^2 \cdot S^5 \\
    S^{(n)} = - \left. \left( \frac{\partial^{n+1}\Omega}{\partial\mu^{n+1}} \right)_{v(\mathbf{r})} \right|_{N = N_0}
           &= \frac{-\sum_{k=1}^{n-1} S^k \cdot B_{n,k}
              \left(\eta^{(1)}, \eta^{(2)}, ..., \eta^{(n-k+1)} \right)}{B_{n,n}\left( \eta^{(1)}\right)}

The first derivative just gives back the number of electrons in the reference state. The second derivative of
the grand potential with respect to the chemical potential gives the global chemical softness, which is just
the reciprocal of the chemical hardness we obtained in the canonical ensemble.
The third- and higher-order derivatives of the grand potential with respect to the chemical potential define the
global hyper-softnesses.
The first global hyper-softness, :math:`S^{(2)}`, appears in a few applications, but none of these higher-order
reactivity indicators have not been explored much in literature.

ChemTools will evaluate global hyper-softnesses to all orders if those derivatives exist in the energy model you
have selected.
Also, in ChemTools, the expressions for reactivity indicators in the grand canonical ensemble are evaluated from the
reactivity indicators in the canonical ensemble. (This is done so that the user does not have to select a model
for grand potential, :math:`\Omega`, analogous to the energy model, :math:`E\left(N\right)`, that users must
specify in the canonical ensemble.)
This leads to rather complicated expressions for the global hyper-softnesses, which are derived from the inverse
function theorem for second- and higher-order derivatives.
Please refer to :ref:`derivation_global_softness` for details.


.. _global_derived_indicators:

**Derived Global Reactivity Descriptors**

In contrast to the fundamental reactivity descriptors presented above, derived reactivity indicators are not
response functions. Instead, they are motivated by handwaving arguments or inspired by empirical correlations.
Many of the most important derived reactivity indicators are based on the  maximum number of electrons that can
be accepted by the system, :math:`N_{\text{max}}`. That is,

 .. math:: N_{\text{max}} = \underbrace {\min }_N E(N)

It is also convenient to define :math:`\Delta N_{\text{max}}` as the difference between the number of electrons
in the reference system, :math:`N_0`, and the maximum number of bound electrons,

 .. math:: \Delta N_{\text{max}} = N_{\text{max}} - N_0


**Electrophilicity Index** :math:`\omega_{\text{electrophilicity}}` measures the ability
of an agent to accept electrons from the environment. However, in contrast to electron affinity, :math:`EA`,
which measures the energy lowering from adding one electron to the system, the electrophilicity index
measures the energy lowering due to maximal electron flow
(which may be either less or more than one) from the environment,

 .. math:: \omega_{\text{electrophilicity}} = E(N_0) - E(N_0 + \Delta N_{\text{max}}) = E(N_0) - E(N_{\text{max}})

This definition is sensible for systems that can bind additional electrons, i.e. :math:`\Delta N_{\text{max}} > 0`,
and have a positive electrophilicity as a result.
By convention, we decided that for systems with :math:`\Delta N_{\text{max}} < 0`, the electrophilicity should be defined
as negative. This gives an extended definition for the electrophilicity

 .. math::
    {\Delta E}_{\text {electrophile}} =
    \text{sgn}\left(N_{\text{max}} - N_0\right) \times \left(E(N_0) - E(N_{\text{max}})\right)

Here :math:`\text{sgn}(x)` is the sign of the argument :math:`x`. I.e.,

 .. math::
    \text{sgn}(x) =
     \begin{cases}
      1 \qquad & x > 0 \\
      0 \qquad & x = 0 \\
      -1 \qquad & x < 0
     \end{cases}

With this extended definition, larger values of :math:`\omega_{\text{electrophilicity}}` are associated with stronger
electrophiles, that is,

 .. math::
    \omega_{\text{electrophilicity}} =
     \begin{cases}
      \gg 0 \qquad & \text{strong electrophile} \\
      \approx 0 \qquad & \text{(very) weak electrophile} \\
      < 0 \qquad & \text{nucleophile/not electrophile}
     \end{cases}


**Nucleophilicity Index** :math:`\omega_{\text{nucleophilicity}}` measures the capacity of a reagent to donate
electrons to its environment. Unlike the electrophilicity, there is not a strong consensus about how one should
define the nucleophilicity. None of the existing definitions is very convincing, perhaps because electrostatic
interactions are often important for nucleophiles, which tend to be locally harder (because nucleophilic sites
tend to be positively charged) than electrophiles. A reasonable definition for nucleophilicity, perhaps first
suggested by Contreras et al., is the effective ionization potential. The effective ionization potential is
closely related to the adiabatic ionization potential: the idea is to include the effects of the environment
(e.g. solvation, surrounding molcules) and geometric relaxation. Given the simplifications inherent in this
description, it is usually reasonable to  simply use the vertical ionization potential,

 .. math:: \omega_{\text{nucleophilicity}} \sim I_{\text{eff}} \sim I

Sometimes the electrofugality (which is closely related) is also identified as the nucleophilicity,

 .. math:: \omega_{\text{nucleophilicity}} \sim \frac{\left(3I - A\right)^2}{8\left(I - A\right)}


**Nucleofugality Index** :math:`\nu_{\text{nucleofugality}}` measures the ability of the system to be a nucleofuge,
i.e., to be a leaving group which takes an electron with it.
Nucleofugality is defined as the energetic penalty that a nucleofuge with :math:`\Delta N_{\text{max}} < 1` must pay
when it receives an entire electron, :math:`\Delta N = 1`, from the molecule it leaves behind,

.. by the energy penalty associated with forcing a molecular fragment to accept an electron; the lower
.. values of :math:`\nu_{\text{nucleofugality}}` are associated with high nucleofugality.

 .. math:: \nu_{\text{nucleofugality}} = E(N_0 + 1) - E(N_0 + \Delta N_{\text{max}}) = E(N_0 + 1) - E(N_{\text{max}})

Nucleofuges with :math:`\Delta N_{\text{max}} > 1` are exceptional leaving groups that readily accept a full electron
from their environment under almost any circumstance.
We therefore define them so that they have negative values of nucleofugality index, motivating the extended definition

 .. math::
    {\Delta E}_{\text {nucleofuge}} = \text{sgn}\left(N_0 + 1 - N_{\text{max}}\right) \times \left(E(N_0 + 1) - E(N_{\text{max}})\right)

In general, small (ideally negative) values of :math:`\nu_{\text{nucleofugality}}` are associated with better leaving groups.
Ergo,

 .. math::
    \nu_{\text{nucleofugality}} =
     \begin{cases}
      \gg 0 \qquad & \text{poor nucleofuge} \\
      \approx 0 \qquad & \text{good nucleofuge} \\
      < 0 \qquad & \text{exceptional nucleofuge}
     \end{cases}


**Electrofugality Index** :math:`\nu_{\text{electrofugality}}` measures the ability
of the system for being a electrofuge, i.e., to be a leaving group which leaves an electron behind with the molecular remnant.
Electrofugality is defined as the energetic penalty that an electrofuge with, :math:`\Delta N_{\text{max}} > -1` must
pay when it donates an entire electron, :math:`\Delta N = -1`, to the molecule it leaves behind,

 .. math:: \nu_{\text{electrofugality}} = E(N_0 - 1) - E(N_0 + \Delta N_{\text{max}}) = E(N_0 - 1) - E(N_{\text{max}})

Electrofuges with :math:`\Delta N_{\text{max}} < -1` are exceptional leaving groups that eagerly donate a full electron
to their environment. We therefore define them so that they have negative values of electrofugality index,
motivating the extended definition

 .. math::
    {\Delta E}_{\text {electrofuge}} =
    \text{sgn}\left(N_{\text{max}} - N_0 + 1\right) \times \left(E(N_0 - 1) - E(N_{\text{max}})\right)

In general, small (ideally negative) values of electrofugality index are associated with better leaving groups.
Ergo,

 .. math::
    \nu_{\text{electrofugality}} =
     \begin{cases}
      \gg 0 \qquad & \text{poor electrofuge} \\
      \approx 0 \qquad & \text{good electrofuge} \\
      < 0 \qquad & \text{exceptional electrofuge}
     \end{cases}


.. _global_energy_models:

**Energy Models**

Calculating fundamental and derived global reactivity descriptors requires that one choose a model for the dependence of the
electronic energy on the number of electrons at fixed external potential, :math:`v\left(\mathbf{r}\right)`.
Details about how ChemTools evaluates global reactivity descriptors using its built-in and user-defined energy models can be found below.

.. toctree::
   :maxdepth: 2

   sci_doc_global_linear
   sci_doc_global_quadratic
   sci_doc_global_rational
   sci_doc_global_exponential
   sci_doc_global_general

