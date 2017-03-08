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

Local Descriptive Tools :class:`chemtools.tool.localtool`
#########################################################

Local descriptive tool :math:`p (\mathbf{r})` assigns a value to every point in space.
These tools are used to study the reactivity of a molecule at point :math:`\mathbf{r}`.

In conceptual DFT, local descriptors arise as functional derivatives of :ref:`global descriptive tools <global_tools>`
with respect to local quantities, typically the external potential :math:`v(\mathbf{r})`, at fixed number of
electrons :math:`N`. Here, :ref:`linear <local_linear>` and :ref:`quadratic <local_quadratic>` energy models are considered,
and the corresponding local descriptors will be introduced.

The functional derivatives of :ref:`energy and its derivatives <energy_derivatives>`
with respect to external potential :math:`v(\mathbf{r})` at fixed number of
electrons :math:`N` results in **electron density** :math:`\rho(\mathbf{r})`,
**Fukui function** :math:`f(\mathbf{r})`, :math:`2^{\text{nd}}` **-order Fukui function** :math:`f^{(2)}(\mathbf{r})`
commonly refered to as **dual descriptor** :math:`\Delta f(\mathbf{r})`
and :math:`n^{\text{th}}` **-order Fukui function** :math:`f^{(n)}(\mathbf{r})` as local descriptors:

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

 .. todo::
    * What are the names for the derivatives of :math:`\Omega` listed below?

On the other hand, the functional derivative of :ref:`grand potential and its derivatives <grand_potential_derivatives>`
with respect to external potential
:math:`v(\mathbf{r})` at fixed chemical potential :math:`\mu` results in **electron density** :math:`\rho(\mathbf{r})`, **local softness** :math:`s(\mathbf{r})`,
and :math:`2^{\text{nd}}` **-order local softness** :math:`s^{(2)}(\mathbf{r})`, and
:math:`n^{\text{th}}` **-order local softness** :math:`s^{(n)}(\mathbf{r})`
as local descriptors for a grand canonical ensemble:

 .. math::

    \rho(\mathbf{r}) = s^{(0)}(\mathbf{r}) &= {\left( \frac{\delta \Omega}{\delta v(\mathbf{r})} \right)_{\mu}}  \\
    s(\mathbf{r}) = s^{(1)}(\mathbf{r}) &= {\left( \frac{\delta S}{\delta v(\mathbf{r})} \right)_{\mu}}
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

For details of expressing local softness expressions in terms of Fukui functions using Faà di Bruno's formula,
please refer to :ref:`derivation_local_softness`.

To derive a working expression for these local quantities, the functional
derivative of :ref:`energy and its derivatives <energy_derivatives>` with respect to external potential
:math:`v(\mathbf{r})` at fixed number of electrons :math:`N` should be calculated.
Considering that,

 .. math::

    \frac{\delta E_{N_0 \pm k}}{\delta v(\mathbf{r})} = \rho_{N_0 \pm k}(\mathbf{r})

for a general energy model the Fukui function expressions can be derived using chain-rule:

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

In the following sections, various energy models are considered and their corresponding local reactivity descriptors are obtained.
It is import to note that, in contrast to global energy models, the local counterparts of rational and exponential energy models
are not included, because these would **not** be properly normalized.
Please see `J. Chem. Theory Comput., 2016, 12 (12), pp 5777–5787  <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00494>`_ for more information.

.. toctree::
   :maxdepth: 2

   sci_doc_local_linear
   sci_doc_local_quadratic
