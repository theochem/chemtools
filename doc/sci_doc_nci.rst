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


.. _nci:

Non-Covalent Interactions (NCI) :class:`chemtools.tool.analysis.nci`
####################################################################

NCI (Non-Covalent Interactions) is a tool to visualize an identify non-covalent interactions,
such as van der Waals interactions, hydrogen bonds, and steric clashes. The index is based on
the density and its derivatives.

The reduced density gradient, which is calculated using the density and its first derivative,

 .. math::
    s\left(\mathbf{r}\right) = \frac{1}{3\left(2\pi ^2 \right)^{1/3}}
    \frac{\lvert \nabla \rho\left(\mathbf{r}\right) \rvert}{\rho\left(\mathbf{r}\right)^{4/3}}

is a dimensionless function, used to describe the deviation from a homogeneous electron distribution.
In regions far from the molecule, in which the density decays to zero exponentially, the reduced gradient
will have very large positive values. In regions of covalent bonding and non-covalent interactions, the
reduced gradient have very small values, close to zero.

To study the non-covalent interactions the plot of the reduced gradients, :math:`s\left(\mathbf{r}\right)`,
versus :math:`\rho\left(\mathbf{r}\right)` plot is examined at low reduced gradient regions (Figure 1).
For the separate molecules, the top left points, corresponding to a small density and large reduced density gradient,
correspond to the tail regions of the density, far from the nuclei. The points at the lower right correspond to covalent
bonds, with higher density, but a low reduced density gradient.

to be continued...

**References:**
  * `Erin R. Johnson, E. R., Keinan S., Mori-Sanchez, P. Contreras-Garcia J., Cohen, A. J.and Yang, J. Am. Chem. Soc.(2010), 132, 6498 <http:/pubs.acs.org/doi/abs/10.1021/ja100936w>`_.
