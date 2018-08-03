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

Local Descriptive Tools
#######################

Local descriptors, :math:`p_{\text{local}} (\mathbf{r})`, assign a value to every point in space.
These values then provide insight into the reactivity of the molecule at point :math:`\mathbf{r}` in space.
As with the :ref:`global_tools`, local reactivity descriptors may be classified as either
:ref:`fundamental <local_fundamental_indicators>` or :ref:`derived <local_derived_indicators>`.
Also similar to global reactivity descriptors, evaluating local descriptors requires specification of an
appropriate model for the dependence of the electronic energy on the number of electrons, i.e., :math:`E\left(N\right)`.
ChemTools has built-in support for the linear and quadratic energy models, as well as the capability to compute
the local descriptors associated with a user-defined energy model.


.. _local_fundamental_indicators:

**Fundamental Local Reactivity Descriptors**

In conceptual DFT, the fundamental local descriptors are functional derivatives of
:ref:`fundamental global descriptive tools <global_fundamental_indicators>` with respect to local quantities, typically the
external potential :math:`v(\mathbf{r})`, but occasionally also the electron density :math:`\rho(\mathbf{r})`.

In the canonical ensemble, fundamental local descriptors arise by differentiating
:ref:`energy and its derivatives <energy_derivatives>` (with respect to the number of electrons),
with respect to the external potential :math:`v(\mathbf{r})` at constant number of electrons :math:`N`.

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

The derivative of the energy with respect to the external potential is, by the Hellmann-Feynman theorem,
equal to the ground-state electron density, :math:`\rho(\mathbf{r})`. :cite:`Hellmann1937,Feynman1939PR`
The derivative of the electronic chemical potential with respect to the external potential,
at fixed electron number, is the **Fukui function**, :math:`f(\mathbf{r})`. :cite:`Fukui1952JCP,Parr1984JACS`
The Fukui function is the fundamental regioselectivity descriptor in conceptual DFT, and indicates the
areas of the molecule that are best able to accept/donate electrons.
The derivative of the chemical hardness with respect to external potential at constant :math:`N` is the
**dual descriptor**, :math:`\Delta f(\mathbf{r})`. :cite:`Morell2005JPC`
The dual descriptor is positive in electrophilic regions negative in nucleophilic regions. It is especially
useful for reactions where electrons are simultaneously accepted/donated electrons (like concerted pericyclic
reactions) or for describing ambiphilic reagents. Higher-order local reactivity descriptors—corresponding to
the derivatives of hyper-hardnesses with respect to the external potential—are called **hyper-Fukui functions** or
:math:`\mathbf{n^{\text{th}}}` **-order Fukui function**, :math:`f^{(n)}(\mathbf{r})`.
While the computational utility of hyper-Fukui functions has not been established, they can be computed with ChemTools.

In the grand canonical ensemble, fundamental local descriptors arise by differentiating
:ref:`grand potential and its derivatives <grand_potential_derivatives>` (with respect to the electronic chemical potential)
with respect to the external potential :math:`v(\mathbf{r})` at constant chemical potential :math:`\mu`.

 .. math::

    \rho(\mathbf{r}) = s^{(0)}(\mathbf{r}) &= {\left( \frac{\delta \Omega}{\delta v(\mathbf{r})} \right)_{\mu}}  \\
    s(\mathbf{r}) = s^{(1)}(\mathbf{r}) &= -{\left( \frac{\delta N}{\delta v(\mathbf{r})} \right)_{\mu}}
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

The derivative of the grand potential with respect to the external potential at constant chemical potential
is the ground-state electron density, :math:`\rho(\mathbf{r})`.
The derivative of minus the number of electrons with respect to the external potential, at fixed electron number,
gives the **local softness**, :math:`s(\mathbf{r})`. :cite:`Yang1985PNAS`
The local softness is the fundamental regioselectivity descriptor for open systems, and is especially useful when
comparing the relative electrophilicity/nucleophilicity of reactive sites of different molecules, especially if
those molecules are somewhat different in size. The derivative of the global softness with respect to external
potential at constant :math:`\mu` is the **dual local softness**, :math:`s^{(2)}(\mathbf{r})`.
The dual local softness is positive in electrophilic regions negative in nucleophilic regions.
It has similar properties and amplicability to the dual descriptor. Higher-order local reactivity
descriptors—corresponding to the derivatives of global hyper-softnesses with respect to the external potential—are
called local **hyper-softnesses** or :math:`\mathbf{n^{\text{th}}}` **-order local softness**, :math:`s^{(n)}(\mathbf{r})`.
While the computational utility of local hyper-softnesses have not been established,
they can be computed with ChemTools.

In ChemTools, the expressions for reactivity indicators in the grand canonical ensemble are evaluated from the
reactivity indicators in the canonical ensemble. This leads to rather complicated expressions for the local
hyper-softnesses. For details about how expressions for the local hyper-softnesses can be expressed in terms of the
hyper-Fukui functions using Faà di Bruno’s formula, please refer to :ref:`derivation_local_softness`.


.. _local_derived_indicators:

**Derived Local Reactivity Descriptors**

Many derived local reactivity descriptors are defined as the local response of a global derived descriptors.
I.e., the derivative of any global reactivity descriptor with respect to the external potential, :math:`v(\mathbf{r})`
at constant :math:`N` or :math:`\mu` is a local reactivity descriptor.
Most other derived descriptors are defined as products of several global and local descriptors.


**Energy Models**

Calculating fundamental and derived global reactivity descriptors requires that one choose a model for the dependence
of the electronic energy upon the number of electrons at fixed external potential, :math:`v(\mathbf{r})`.
For the user-selected reference states with electron numbers :math:`N_0,N_1,N_2,...,N_n` and electronic energies
:math:`E\left(N_0\right),E\left(N_1\right),E\left(N_2\right),...,E\left(N_n\right)`, a model for the electron density
and Fukui functions as a function of the number of electrons can be derived using chain-rule:

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

where the known electron densities for the reference states are

 .. math::

    \rho_{N_k}(\mathbf{r}) = \left(\frac{\delta E\left(N_k\right)}{\delta v(\mathbf{r})}\right)_{N=N_k} \qquad k=1,2,\dots,n


ChemTools provide explicit implementations for local reactivity descriptors computed using the

.. toctree::
   :maxdepth: 2

   sci_doc_local_linear
   sci_doc_local_quadratic

The rational model and the exponential model do not have built-in support because the local reactivity
indicators (even the electron density!) do not satisfy the appropriate normalization constraints for these models.
Please see :cite:`Heidar-Zadeh2016JCTC` for more information.

 .. However, the local reactivity descriptors for these, and other, models can be computed using the general linear model.
 .. In the general linear model, the (hyper) Fukui functions are derived using the chain rule,
