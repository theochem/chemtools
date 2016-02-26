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

The functional derivatives of :ref:`energy and its derivatives <energy_derivatives>`
with respect to external potential :math:`v(\mathbf{r})` at fixed number of
electrons :math:`N` results in **electron density** :math:`\rho(\mathbf{r})`,
**Fukui function** :math:`f(\mathbf{r})`, **dual descriptor** :math:`d(\mathbf{r})`
and :math:`n^{\text{th}}` **-order Fukui function** :math:`f^{(n)}(\mathbf{r})` as local descriptors:

 .. math::

    \rho(\mathbf{r}) = f^{(0)}(\mathbf{r}) &= {\left( \frac{\delta E}{\delta v(\mathbf{r})} \right)_N} && \\
    f(\mathbf{r}) = f^{(1)}(\mathbf{r}) &= {\left( \frac{\delta \mu}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial E}{\partial N} \right)_{v(\mathbf{r})}} \right)_N} \\
    d(\mathbf{r}) = f^{(2)}(\mathbf{r}) &= {\left( \frac{\delta \eta}{\delta v(\mathbf{r})} \right)_N}
              &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                  {\left( \frac{\partial^2 E}{\partial N^2} \right)_{v(\mathbf{r})}} \right)_N} \\
    f^{(n)}(\mathbf{r}) &= {\left( \frac{\delta \eta^{(n-1)}}{\delta v(\mathbf{r})} \right)_N}
             &&= {\left( \frac{\delta}{\delta v(\mathbf{r})}
                 {\left( \frac{\partial^n E}{\partial N^n} \right)_{v(\mathbf{r})}} \right)_N}

 .. todo::
    * Is 2nd derivative always called dual descriptor? Is it denoted with d(r)?
    * Is there any name for higher otders?
    * What are the names for the derivatives of :math:`\Omega` listed below?

On the other hand, the functional derivative of :ref:`grand potential and its derivatives <grand_potential_derivatives>`
with respect to external potential
:math:`v(\mathbf{r})` at fixed chemical potential :math:`\mu` results in **???**, **local softness** :math:`s(\mathbf{r})`,
and :math:`2^{\text{nd}}` **-order local softness** :math:`S^{(2)}(\mathbf{r})`, and
:math:`n^{\text{th}}` **-order local softness** :math:`S^{(n)}(\mathbf{r})`
as local descriptors for a grand canonical ensemble:

 .. math::

    ?(\mathbf{r}) = s^{(0)}(\mathbf{r}) &= {\left( \frac{\delta \Omega}{\delta v(\mathbf{r})} \right)_{\mu}}  \\
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

For details of expressing local softness expressions in terms of Fukui functions using chain rule,
please refer to :ref:`derivation_local_softness`.

To derive a working expression for these local quantities, the functional
derivative of :ref:`energy and its derivatives <energy_derivatives>` with respect to external potential
:math:`v(\mathbf{r})` at fixed number of electrons :math:`N` should be calculated:


 .. math::

    f^{(n)}(\mathbf{r}) &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n E_{\text{model}}}
                           {\partial N^n}\right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\delta}{\delta v(\mathbf{r})}{\left(\frac{\partial^n}
       {\partial N^n} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_{v(\mathbf{r})}}\right)_N} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\frac{\delta}
       {\delta v(\mathbf{r})} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\}) \right)_N} \right)_{v(\mathbf{r})}} \\
    &= {\left(\frac{\partial^n}{\partial N^n}{\left(\sum_{i=1}^n \left( \left.\frac{\partial E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial \alpha_{i}}\right|_{N=N_0} \cdot \sum_k \frac{\partial \alpha_i}{\partial E_{N_0 \pm k}} \rho_{N_0 \pm k}(\mathbf{r})
       \right)\right)} \right)_{v(\mathbf{r})}} \\
    &= \sum_{i=1}^n \left(\left. \underbrace {\frac{\partial^{n+1} E(N; \{\alpha_1, \alpha_2, ..., \alpha_n\})}
       {\partial N^n \partial\alpha_{i}}}_{\text {evaluated analytically}}\right|_{N=N_0} \cdot
       \sum_k \underbrace {\frac{\partial \alpha_i}{\partial E_{N_0 \pm k}}}_{\text{evaluated}\atop\text{numerically}}
       \rho_{N_0 \pm k}(\mathbf{r})\right)_{v(\mathbf{r})}

In the following sections, this expression will be derived for :ref:`linear <local_linear>`
and :ref:`quadratic <local_quadratic>` energy models,
and the corresponding local descriptors will be introduced.

.. _local_linear:

Linear Local Model: :class:`chemtools.tool.localtool.LinearLocalTool`
---------------------------------------------------------------------

Complementing the :ref:`linear global model <linear_energy>`, the piece-wise linear
energy model :math:`E\left(N\right) = a + b N` alongside the electron density of
systems with :math:`N_0 - 1`, :math:`N_0` and :math:`N_0 + 1` electrons,
:math:`{\{\rho_{N_0 - 1}\left(\mathbf{r}\right), \rho_{N_0}\left(\mathbf{r}\right), \rho_{N_0 + 1}\left(\mathbf{r}\right)\}}`
, are used to calculate linear local descriptors. Due to non-differentiability
of energy model at integer number of electrons, these tools are calculated for :math:`N \geq N_0`
and :math:`N \leq N_0` separately.

Substituting the derivatives of :math:`E\left(N\right)` with respect to parameters :math:`a`
and :math:`b`, as well as the derivatives of :math:`a` and :math:`b` with respect to
:math:`E\left(N_0 - 1\right)`, :math:`E\left(N_0\right)` and :math:`E\left(N_0 + 1\right)` in
Eq. (????), results in Fukui functions. For detailed derivation, please refer to
:ref:`derivation_linear_fukui_function`.

For :math:`N \leq N_0`, the Fukui function from below :math:`f^-\left(\mathbf{r}\right)` is:

 .. math::

    f^{(0)}(\mathbf{r}) &= \rho_{N_0}\left(\mathbf{r}\right) \\
    f^-\left(\mathbf{r}\right) = f^{(1)}(\mathbf{r}) &=
       \rho_{N_0}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right) \\
    f^{(n)}(\mathbf{r}) &= 0 \text{  for   } n \geq 2

For :math:`N \geq N_0`, the Fukui function from above :math:`f^+\left(\mathbf{r}\right)` is:

 .. math::

    f^{(0)}(\mathbf{r}) &= \rho_{N_0}\left(\mathbf{r}\right) \\
    f^+\left(\mathbf{r}\right) = f^{(1)}(\mathbf{r}) &=
       \rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0}\left(\mathbf{r}\right) \\
    f^{(n)}(\mathbf{r}) &= 0 \text{  for   } n \geq 2


Commonly, the average Fukui function :math:`f^0\left(\mathbf{r}\right)` is used:

 .. math::

    f^0\left(\mathbf{r}\right) = \frac{f^+\left(\mathbf{r}\right) + f^-\left(\mathbf{r}\right)}{2} =
             \frac{\rho_{N_0 + 1}\left(\mathbf{r}\right) - \rho_{N_0 - 1}\left(\mathbf{r}\right)}{2}


Dual descriptor is defined as:

 .. math::

    d\left(\mathbf{r}\right) = f^+\left(\mathbf{r}\right) - f^-\left(\mathbf{r}\right) =
           \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0 - 1}\left(\mathbf{r}\right) + \rho_{N_0 - 1}\left(\mathbf{r}\right)

 .. todo::
    * This is not really dual descriptor for linear model. Technically the dual descriptor is zero for linear model,
      but the dual descriptor for quadratic model happens to be f+(r) - f-(r).
      Does this need to be clarified?


.. _local_quadratic:

Quadratic Local Model: :class:`chemtools.tool.localtool.QuadraticLocalTool`
---------------------------------------------------------------------------

Complementing the :ref:`quadratic global model <quadratic_energy>`, the :math:`E\left(N\right) = a + b N + c N^2`
energy model alongside the electron density of systems with :math:`N_0 - 1`, :math:`N_0` and
:math:`N_0 + 1` electrons,
:math:`{\{\rho_{N_0 - 1}\left(\mathbf{r}\right), \rho_{N_0}\left(\mathbf{r}\right), \rho_{N_0 + 1}\left(\mathbf{r}\right)\}}`
, are used to calculate quadratic local descriptors.

Substituting the derivatives of :math:`E\left(N\right)` with respect to parameters :math:`a`,
:math:`b` and :math:`c`, as well as the derivatives of :math:`a`, :math:`b` and :math:`c` with respect to
:math:`E\left(N_0 - 1\right)`, :math:`E\left(N_0\right)` and :math:`E\left(N_0 + 1\right)` in
Eq. (????), results in Fukui functions. For detailed derivation, please refer to
:ref:`derivation_quadratic_fukui_function`.


 .. math::

    f^{(0)}(\mathbf{r}) &= \rho_{N_0}\left(\mathbf{r}\right) \\
    f\left(\mathbf{r}\right) = f^{(1)}\left(\mathbf{r}\right) &=
     \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2} \\
    \eta\left(\mathbf{r}\right) = f^{(2)}\left(\mathbf{r}\right) &=
      \rho_{N_0 + 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
      \rho_{N_0 - 1}\left(\mathbf{r}\right) \\
    f^{(n)}(\mathbf{r}) &= 0 \text{  for   } n \geq 3

The local softness is easily found by substituting the Fukui functions in Eq. (????):

 .. math::

    s\left(\mathbf{r}\right) = s^{(1)}\left(\mathbf{r}\right) &= S \cdot f\left(\mathbf{r}\right) =
     \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2 \eta} =
     \frac{\rho_{N_0+1}\left(\mathbf{r}\right) - \rho_{N_0-1}\left(\mathbf{r}\right)}{2 \left(IP - EA\right)} \\
    s^{(2)}\left(\mathbf{r}\right) &= S^{2} \cdot f^{(2)}(\mathbf{r}) + S^{(2)} \cdot f(\mathbf{r}) \\
     &= \frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
        \rho_{N_0 + 1}\left(\mathbf{r}\right)}{\eta^2} =
        \frac{\rho_{N_0 - 1}\left(\mathbf{r}\right) - 2 \rho_{N_0}\left(\mathbf{r}\right) +
        \rho_{N_0 + 1}\left(\mathbf{r}\right)}{\left(IP - EA\right)^2} \\
    s^{(n)}(\mathbf{r}) &= 0 \text{  for   } n \geq 3

where :math:`\eta` represents the gloabl chemical hardness in quadratic energy model.


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
