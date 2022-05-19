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


.. _dftbased_it:

Information-Theoretic Descriptors
=================================


Information-theoretic descriptors are often used in quantitative structure property relationships
and to analyze bonding in molecules. One of the most useful tools related to information-theoretic
analysis are the steric energy, steric potential, and steric force.
These quantities are closely linked to both the non-covalent interactions analysis (NCI) and the
density-overlaps regions indicator (DORI), and can be used to perform similar analysis.

A few information-theoretic descriptors are so fundamental that they are given special names are
described here.


**Shannon Entropy** is the "uncertainty" associated with the electron density. It is larger when
the electron density is more uniform.

 .. math::
    S = \int \rho\left(\mathbf{r}\right) \ln \rho\left(\mathbf{r}\right) d\mathbf{r}


**Shannon Entropy Density** is the local Shannon entropy, equal to the integrand of Shannon entropy:

 .. math::
    s\left(\mathbf{r}\right) = \rho\left(\mathbf{r}\right) \ln \rho\left(\mathbf{r}\right)


**Shannon Shape Entropy** is the Shannon entropy of electron density per particle or shape function
:math:`\sigma\left(\mathbf{r}\right)`:

 .. math::
    S = \int \sigma\left(\mathbf{r}\right) \ln \sigma\left(\mathbf{r}\right) d\mathbf{r}


**Shannon Shape Entropy Density** is the local Shannon entropy, equal to the integrand of Shannon
entropy:

 .. math::
    s(\mathbf{r}) = \sigma(\mathbf{r}) \ln\sigma(\mathbf{r})


**Fisher Information** is the Fisher entropy of locality is (by abuse/extension of mathematical
equations) is also associated with the electron density,

 .. math::
    I = \int \frac{\nabla\rho(\mathbf{r}) \cdot \rho(\mathbf{r})}{\rho(\mathbf{r})} d\mathbf{r}

The Fisher information is smaller when the electron density is more uniform, and is closely related
to the Weizsacker kinetic energy and the steric energy.


**Fisher Information Density** is the local Fisher information, which is the integrand of the
Fisher information:

 .. math::
    i(\mathbf{r}) = \frac{\nabla\rho(\mathbf{r}) \cdot \rho(\mathbf{r})}{\rho(\mathbf{r})}


**Fisher Shape Information** is the Fisher entropy of locality for the shape function, or
density-per-particle:

 .. math::
    I = \int \frac{\nabla\sigma(\mathbf{r})\cdot\sigma(\mathbf{r})}{\sigma(\mathbf{r})} d\mathbf{r}


**Fisher Shape Information Density** is the local Fisher shape information, defined as the
integrand in the Fisher shape information:

 .. math::
    i(\mathbf{r}) = \frac{\nabla\sigma(\mathbf{r}) \cdot \sigma(\mathbf{r})}{\sigma(\mathbf{r})}


**Ghosh-Berkowitz-Parr Entropy** is based on a local thermodynamic treatment of the electron
density:

 .. math::
    S_\text{GBP} = \int \tfrac{3}{2} \rho(\mathbf{r}) \ln\left(
    \frac{\tau_\text{PD}(\mathbf{r})}{\tau_\text{TF}(\mathbf{r})}\right) d\mathbf{r}

To accommodate the range of descriptors that have been associated with the Ghosh-Berkowitz-Parr
entropy, the user can pass any local kinetic energy to the expression for the Ghosh-Berkowitz-Parr
expression, but the default expression is the one they derived in their paper, which corresponds to
a factor :math:`-\tfrac{1}{8}` of the Laplacian term.


**Ghosh-Berkowitz-Parr Local Entropy** is the integrand of the preceding Ghosh-Berkowitz-Parr
expression for the global entropy.


**Ghosh-Berkowitz-Parr Shape Entropy** is the shape-function analogue of the Ghosh-Berkowitz-Parr
entropy (up to a constant additive factor)

 .. math::
    S_\text{GBP} = \int \tfrac{3}{2} \sigma(\mathbf{r}) \ln\left(N_\text{spin}^{\tfrac{2}{3}}
    \frac{\tau_\text{PD}(\mathbf{r})}{\tau_\text{TF}(\mathbf{r})}\right) d\mathbf{r}

To accommodate the range of descriptors that have been associated with the Ghosh-Berkowitz-Parr
entropy, the user can pass any local kinetic energy to the expression for the Ghosh-Berkowitz-Parr
expression, though the dependence on the number of electrons in the spin-channel one considers will
be incorrect for some (uncommon) choices.


**Ghosh-Berkowitz-Parr Local Shape Entropy** is the integrand of the preceding Ghosh-Berkowitz-Parr
expression for the global shape entropy.


**Steric Energy** is the kinetic energy of a system if the electrons were bosons; ergo it quantifies
the inherent “kinetic energy pressure” with quantum (Pauli exclusion) effects are omitted.
It is equal to the Weizsacker kinetic energy:

 .. math::
    E_\text{steric} = \int \tau_\text{W}(\mathbf{r}) d\mathbf{r} = \int
    \frac{\lvert \nabla \rho(\mathbf{r}) \rvert ^2}{8 \rho(\mathbf{r})} d\mathbf{r}


**Steric Potential** is the functional derivative of the steric energy.
The steric potential plus the Pauli potential equals the Kohn-Sham potential.

 .. math::
    v_\text{steric}(\mathbf{r}) = \frac{\delta E_\text{steric}}{\delta \rho(\mathbf{r})} &=
    -\frac{1}{8} \frac{\nabla\rho(\mathbf{r}) \cdot \nabla\rho(\mathbf{r})}{\rho^2(\mathbf{r})}
    -\frac{1}{4} \nabla \cdot \frac{\nabla\rho(\mathbf{r})}{\rho(\mathbf{r})} \\ &=
    -\frac{1}{8} \frac{\nabla\rho(\mathbf{r}) \cdot \nabla\rho(\mathbf{r})}{\rho^2(\mathbf{r})}
    +\frac{1}{4} \frac{\nabla\rho(\mathbf{r}) \cdot \nabla\rho(\mathbf{r})}{\rho^2(\mathbf{r})}
    -\frac{1}{4} \frac{\nabla^2\rho(\mathbf{r})}{\rho(\mathbf{r})} \\ &=
    \frac{1}{8} \frac{\nabla\rho(\mathbf{r}) \cdot \nabla\rho(\mathbf{r})}{\rho^2(\mathbf{r})}
    -\frac{1}{4} \frac{\nabla^2\rho(\mathbf{r})}{\rho(\mathbf{r})}


**Steric Charge Density** is the classical charge distribution that generates the steric potential,

 .. math::
    q_\text{steric}(\mathbf{r}) = \tfrac{-1}{4 \pi} \nabla^2 v_\text{steric}(\mathbf{r})


This is an interesting example of the power of ChemTools, because it requires very high-order
differentiation, which most programs do not support (and where numerical differential is often
inaccurate).


**Steric Charge** is the total steric charge:

 .. math::
    Q_\text{steric} = \int q_\text{steric}(\mathbf{r}) d\mathbf{r}
