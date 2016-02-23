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

Local descriptive tools :math:`p (\mathbf{r}) = p_{\text local} \left(\mathbf{r}\right)`
assign a value to every point in space.
In other words, these tools describe how a molecule react at point :math:`\mathbf{r}`.


Conceptual DFT Local Descriptors: :class:`chemtools.tool.localtool`
===================================================================

In conceptual DFT, local descriptors arise as functional derivatives of :ref:`global descriptive tools <global_tools>`
with respect to local quantities, typically the external potential :math:`v(\mathbf{r})` at fixed number of
electrons :math:`N`.

The functional derivatives of energy :math:`E`, chemical potential :math:`\mu`, chemical hardness :math:`\eta`
and :math:`n^{\text{th}}` -order hyper-hardness :math:`\eta^{(n)}` with respect to external potential
:math:`v(\mathbf{r})` at fixed number of
electrons :math:`N` results in **electron density** :math:`\rho(\mathbf{r})`, **Fukui function**
:math:`f(\mathbf{r})`, :math:`2^{\text{nd}}` **-order Fukui function** :math:`f^{(2)}(\mathbf{r})`
and :math:`n^{\text{th}}` **-order Fukui function** as local descriptors:

 .. math::

    \rho(\mathbf{r}) &= {\left( \frac{\delta E}{\delta v(\mathbf{r})} \right)_N} && \\
    f(\mathbf{r}) &= {\left( \frac{\delta \mu}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}} \right)_N} \\
    f^{(2)}(\mathbf{r}) &= {\left( \frac{\delta \eta}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}} \right)_N} \\
    f^{(n)}(\mathbf{r}) &= {\left( \frac{\delta \eta^{(n-1)}}{\delta v(\mathbf{r})} \right)_N}
             &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                 {\left( \frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})}} \right)_N} \text{for } n\geq2

On the other hand, the functional derivative of grand potential :math:`\Omega`, number of electrons :math:`N`,
softness :math:`S`, and :math:`n^{\text{th}}` -order hyper-softness :math:`\S^{(n)}` with respect to external potential
:math:`v(\mathbf{r})` at fixed chemical potential :math:`\mu` results in **???**, **local softness** :math:`s(\mathbf{r})`,
and :math:`2^{\text{nd}}` **-order local softness** :math:`S^{(n)}(\mathbf{r})` as local descriptors
for a grand canonical ensemble:

 .. math::

    ?(\mathbf{r}) &= {\left( \frac{\delta \Omega}{\delta v(\mathbf{r})} \right)_{\mu}} && \\
    s(\mathbf{r}) &= {\left( \frac{\delta S}{\delta v(\mathbf{r})} \right)_{\mu}}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial \Omega}{\partial \mu} \right)_{v(\mathbf{r})}} \right)_{\mu}}
                = S \cdot f(\mathbf{r})  \\
    s^{(2)}(\mathbf{r}) &= {\left( \frac{\delta S}{\delta v(\mathbf{r})} \right)_{\mu}}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial^2 \Omega}{\partial {\mu}^2} \right)_{v(\mathbf{r})}} \right)_{\mu}}
		= S^{2} \cdot f^{(2)}(\mathbf{r}) + S^{(2)} \cdot f(\mathbf{r}) \\
    s^{(n)}(\mathbf{r}) &= {\left( \frac{\delta S^{(n-1)}}{\delta v(\mathbf{r})} \right)_{\mu}}
             &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                 {\left( \frac{\partial^n \Omega}{\partial {\mu}^n} \right)_{v(\mathbf{r})}} \right)_{\mu}}
	       = -\sum_{k=1}^n f^{(k)}(\mathbf{r}) \cdot B_{n,k}\left(S^{(1)}, S^{(2)}, ..., S^{(n-k+1)} \right)

For details of expressing local softness in terms of Fukui functions using chain rule,
please refer to :ref:`derivation_local_softness`. Chain rule has been used to
:ref:`express local softness in terms of Fukui functions <derivation_local_softness>`.

To derive working expressions for these local quantities, the functional derivative of energy model and its derivatives

 .. math::

    f^{(n)}(\mathbf{r}) &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n E_{\text{model}}}
                           {\partial N^n}\right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n}
       {\partial N^n} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\frac{\delta}
       {\delta v(\mathbf{r})} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_N} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \left(\frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}} \cdot \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}} \frac{\delta E_{N_0 \pm k}}
       {\delta v(\mathbf{r})}\right)\right)_N} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \left(\frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}} \cdot \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}} \rho_{N_0 \pm k}(\mathbf{r})
       \right)\right)} \right)_{v(\mathbf{r})}} \\
    &= \sum_{i=1}^n \left( \underbrace {\frac{\partial^{(n+1)} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial N^n \partial\alpha_{i}}}_{\text {evaluated analytically}} \cdot
       \sum_k \underbrace {\frac{\partial \alpha_i}{\partial E_{N_0 \pm k}}}_{\text{evaluated}\atop\text{numerically}}
       \rho_{N_0 \pm k}(\mathbf{r})\right)_{v(\mathbf{r})}


Linear Local Model:
-------------------


Quadratic Local Model:
----------------------


Density-Based Local Descriptors: :class:`chemtools.tool.densitytool`
====================================================================

All the tools for calculating wich the electron density :math:`\rho\left(\mathbf{r}\right)`, gradianet and hessian
of the :math:`N` electron reference state is enough.

**Electron density** :math:`\rho\left(\mathbf{r}\right)` represents ...

**Gradient of electron density** :math:`\nabla \rho\left(\mathbf{r}\right)` represents the first-order partial
derivatives of electron density with respect to coordinates:

 .. math:: \nabla \rho\left(\mathbf{r}\right) =
           \left( \frac{\partial}{\partial x}\mathbf{i}, \frac{\partial}{\partial y}\mathbf{j}, \frac{\partial}{\partial z}\mathbf{k}\right) \rho\left(\mathbf{r}\right)

**Hessian of electron density** :math:`\nabla^2 \rho\left(\mathbf{r}\right)` represents the second-order
partial derivative of electron density with respect to coordinates:



Orbital-Based Local Descriptors: :class:`chemtools.tool.orbitaltool`
====================================================================

All the tools for calculating which the orbital information of the :math:`N` electron reference state is enough.
