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


.. _condensed_tools:

Condensed Descriptive Tools :mod:`chemtools.tool.condensedtool`
###############################################################

Local reactivity indicators indicate the susceptibility of a particular point in space to reactions.
Chemists usually think, however, in terms of the reactivity of atoms and functional groups.
The coarse-graining of pointwise local reactivity indicators into atomic and/or functional group
contributions gives condensed reactivity indicators.

In conceptual DFT, the fundamental local reactivity indicators are derivatives of the electron density with
respect to either the number of electrons or the chemical potential,



Sometimes (e.g., the local electrophilicity), one will multiply one of the above reactivity indicators
by a global reactivity indicator, or consider the sum of two or more local reactivity indicators.
To coarse-grain these descriptors, we introduce a method for partitioning the molecule into atoms/functional groups.
This partitioning is expressed in terms of atomic weighting function (or, occasionally, atomic weighting operators)
in real space,



In the **fragment-of-molecular-response (FMR) approach**, local properties are divided directly,



In the **response-of-molecular-fragment (RMF) approach**, the electron density is condensed into atomic populations,



and the atomic populations are then differentiated,



The FMR and RMF approaches are just two among many different methods for atom-condensing local reactivity indicators;
they give the same results only if the atomic weight functions do not depend on the number of electrons in the molecule,
as is true for the (ordinary) Hirshfeld partitioning, the Voronoi/Becke partitioning, and the Mulliken partitioning.
For more sophisticated partitioning methods like iterative Hirshfeld charges, the FMR and RMF methods give different results,
though we know of no compelling formal or practical reasons to favor one approach over the other.
In some contexts, the RMF approach is easier to compute, since it only requires performing population analysis on several
different charge states of the molecule being studied. For higher-order reactivity indicators, corresponding to :math:`k > 2`
in Eq. (1), the FMR approach seems somewhat simpler.

The fundamental nonlocal reactivity indicators,



can be condensed in a similar way, producing a matrix that expresses the change in the reactivity of one portion of
the molecule in response to a perturbation of a different portion of the molecule.
As before, in the fragment-of-molecular-response (FMR) approach one condenses the nonlocal reactivity indicator directly,



while the response-of-molecular-fragment (RMF) approach one condenses the linear response function,



or the softness kernel,



and then differentiates these quantities with respect to either the number of electrons or the chemical potential,



Nonlocal reactivity indicators depending on three or more points in space, :math:`v\left(\mathbf{r},\mathbf{r'},\mathbf{r"},...\right)`,
can be condensed into tensors, :math:`v_{ABC...}`, using the same strategy.



Evaluating the FMR and RMF condensed reactivity indicators requires that one select an appropriate model for the
dependence of the energy upon the number of electrons. ChemTools can evaluate condensed reactivity indicators for a
general energy model using the formulation in ref. 7, but the most common choices are the linear model and the quadratic model.

In the linear model, the condensed Fukui functions are



in the FMR approach and



in the RMF approach. Here we have used the notation :math:`\rho\left(N;\mathbf{r}\right)` to indicate the :math:`N-` electron
ground-state density and :math:`w_A\left(N;\mathbf{r}\right)` to indicate the atomic weighting function for atom :math:`A`
in the :math:`N-` electron molecule. Similarly we use :math:`N_A\left(N\right)` and :math:`q_A\left(N\right)` to indicate
the population and charge, respectively, of atom :math:`A` in the :math:`N-` electron molecule. In the linear model, the
condensed dual descriptor is technically undefined. In the quadratic model, the condensed Fukui function and condensed
dual descriptor are defined as,



in the FMR approach and



in the RMF approach.

Condensed reactivity indicators corresponding to derivatives with respect to the chemical potential are computed through the
condensed reactivity indicators corresponding to derivatives with respect to electron number.  For example, the condensed local
softness is defined as



and the condensed dual local softness is defined as



These reactivity indicators can be computed using Fukui functions and/or dual descriptors from either the FMR or RMF approaches.
