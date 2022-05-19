..
    : ChemTools is a collection of interpretive chemical tools for
    : analyzing outputs of the quantum chemistry calculations.
    :
    : Copyright (C) 2016-2019 The ChemTools Development Team
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


.. _dftbased_ked:

Kinetic Energy Density (KED)
============================


The kinetic-energy density, or local kinetic energy, is inherently ambiguous, because specifying
the expectation value for the kinetic energy of the electrons, :math:`\tfrac{1}{2m} p^2`, at a
specified point in space requires knowing the momentum and location of the electrons simultaneously,
which is impossible according to the Heisenberg Uncertainty Principle. This has not stopped
scientists from defining a myriad exact, and approximate, definitions for the kinetic energy
density. As with other DFT-based tools, these can be defined in either spin-resolved or
spin-unresolved manner, depending on how the class is initialized.

While kinetic energy densities are not commonly employed on their own, they are common ingredients
in descriptors that are used to elucidate molecular electronic structure and bonding, as well as in
topological analysis and partitioning techniques.


**Positive-definite Kinetic Energy Density:** The manifestly nonnegative kinetic energy density,
which is the key ingredient in many descriptors and also in the

 .. math::
    \tau_\text{PD}(\mathbf{r}) = \tfrac{1}{2}\sum_i^N n_i \rvert\nabla\phi_i(\mathbf{r})\lvert^2

where :math:`\phi_i(\mathbf{r})` and :math:`n_i` are the molecule (spin) orbitals and their
occupation numbers.


**General(ish) Kinetic Energy Density**

 .. math::
    \tau_\text{G}(\mathbf{r}, \alpha) =
    \tau_\text{PD}(\mathbf{r}) + \tfrac{1}{4} (a - 1) \nabla^2\rho(\mathbf{r})

The integral of this kinetic energy density gives the exact non-interacting kinetic energy.
Several of the more common kinetic energy densities in the literature arise as special sub-cases,
including the form of Ghosh, Berkowitz, and Parr (GBP, :math:`a=\tfrac{1}{2}`) which can be
justified by maximizing the entropy of the underlying quasiprobability distribution function and
the form of Yang, Liu, and Wang (YLW, :math:`a=0`), which can be justified by appealing to the
local energy. The form   corresponds to the “Schrödinger” kinetic energy that Bader denoted
:math:`K(\mathbf{r})`, while :math:`a=1` is the positive-definite kinetic energy that Bader denoted
:math:`G(\mathbf{r})`. While only the form :math:`a=0` is is non-negative, but the integral gives
the correct global kinetic energy for any value of :math:`a`.


**“Unambiguous” Kohn-Sham Kinetic Energy Density:** Using the virial theorem, one can derivate a
kinetic-energy density from the Kohn-Sham potential.
The usual expression is not invariant to translation/rotation of the electron density, but an
alternative expression, which requires solving a Poisson equation, avoids this issue:

 .. math::
    \tfrac{-1}{4\pi} \nabla^2 \tau_\text{unambiguous}(\mathbf{r}) =
    \tfrac{3}{8\pi} \nabla \cdot \rho(\mathbf{r}) \nabla v_s(\mathbf{r})

Here, :math:`v_s(\mathbf{r})` is the Kohn-Sham potential. The integral of this kinetic energy
density gives the exact non-interacting kinetic energy.


**Thomas-Fermi Kinetic Energy Density** :math:`\tau_\text{TF}\left(\mathbf{r}\right)`
represents non-interacting kinetic energy density that is exact for a uniform electron gas:

 .. math::
    \tau_\text{TF}(\mathbf{r}) = \begin{cases}
    \tfrac{3}{10} \left(3\pi\right)^{2/3} \rho(\mathbf{r})^{5/3} \text{  for spin}=\alpha + \beta\\
    \tfrac{3}{10} \left(6\pi\right)^{2/3} \rho(\mathbf{r})^{5/3} \text{  for spin}=\alpha, \beta,
    \alpha - \beta \end{cases}


**Von Weizsacker Kinetic Energy Density** :math:`\tau_\text{W}\left(\mathbf{r}\right)`
The kinetic energy density for a system of bosons with the same density as the electron density.
It is exact for one- and two-electron systems (with nondegenerate ground states).

 .. math::
    \tau_\text{W}(\mathbf{r}) =
    \frac{\lvert \nabla \rho(\mathbf{r}) \rvert ^2}{8 \rho(\mathbf{r})} =
    \frac{\nabla\rho(\mathbf{r}) \cdot \nabla\rho(\mathbf{r})}{8 \rho(\mathbf{r})} =
    \tfrac{1}{2} \nabla\sqrt{\rho(\mathbf{r})} \cdot \nabla\sqrt{\rho(\mathbf{r})}


**Gradient Expansion Approximation:** This expression is exact for the kinetic energy density of a
slowly-varying electron gas. It is derived by performing a perturbation expansion about the uniform
electron gas limit, and truncating the expansion at second order. The fourth- and sixth-order terms
in the expansion are known, but rarely used because they are strongly singular near the nucleus.

 .. math::
    \tau_\text{GEA}(\mathbf{r}) = \tau_\text{TF}(\mathbf{r}) +
    \tfrac{1}{9} \tau_\text{W}(\mathbf{r}) + \tfrac{1}{6} \nabla^2\rho(\mathbf{r})


**Empirical Gradient Expansion Approximation:** The gradient expansion approximation with a higher
Weizsacker factor, which often gives more accurate results in practice because atoms and molecules
are far from the slowly-varying limit.

 .. math::
    \tau_\text{empGEA}(\mathbf{r}) = \tau_\text{TF} (\mathbf{r}) +
    \tfrac{1}{5} \tau_\text{W}(\mathbf{r}) + \tfrac{1}{6} \nabla^2 \rho(\mathbf{r})


**Nuclear-Corrected Kinetic Energy Density:** For any kinetic energy expression, this expression
makes a correction near the nucleus, based on the fact the Weizsacker kinetic energy density is
very accurate there. This is particularly useful for functionals that include Laplacian
contributions, which tend to be bad there.

 .. math::
    w(\mathbf{r}) = \sum_{A=1}^{N_\text{atoms}} \exp\left(-
    \frac{\left(Z_A \lvert \mathbf{r} - \mathbf{R}_{A} \rvert\right)^4}{\left(\ln2\right)^3}\right)

 .. math::
    \tau_\text{mod}(\mathbf{r}) = w(\mathbf{r}) \tau_\text{W}(\mathbf{r}) +
    (1 + w(\mathbf{r})) \tau(\mathbf{r})

where :math:`\tau(\mathbf{r})` can be any kinetic energy density.


While the kinetic energy density from this equation tends to resemble the positive-definite kinetic
energy density more strongly than the preceding approximations, its integrated value is not
significantly more accurate.

TODO: (Eventually we should also provide some options for correlation-kinetic energies.)


**Global Kinetic Energy:** We do not recommend usually evaluating the global kinetic energy from
the preceding kinetic energy densities, but we support this by allowing one to evaluate the kinetic
energy corresponding to any kinetic energy density (by simple integration):

 .. math::
    T = \int \tau(\mathbf{r}) d\mathbf{r}
