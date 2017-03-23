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


Conceptual Density Functional Theory
####################################

Conceptual Density Functional Theory (DFT) provides chemists with a hierarchy of well-defined chemical
concepts that contribute to the qualitative understanding and quantitative prediction of chemical reactivity.
When a molecule undergoes a reaction, its number of electrons increases (nucleophilic attack) or decreases
(electrophilic attack). Also, the external potential felt by the electrons in the molecule changes, because
now the electrons are not only attracted to the nuclei of the molecule, but also attracted to the nuclei and
repelled by the electrons of the attacking reagent. Therefore, the susceptibility of a molecule to chemical
reactions is determined by its response to changes in the (a) number of electrons and (b) the external potential.
This hierarchy of molecular responses, which are the fundamental reactivity indicators of conceptual DFT,
is summarized in the following diagram.

 .. image:: ./_static/CDFT.png
     :align: center

:ref:`Global descriptors <global_tools>` measure the overall susceptibility of a system to different types
of reactions, e.g., electrophilic or nucleophilic attacks. The leftmost entries in each row of the above figure
are global descriptors.
The second-from-the-left entries in the figure are local descriptors. :ref:`Local descriptors <local_tools>`
show where a molecule is most susceptible to different types of reagents; they are regioselectivity indicators.
Coarse-graining of local descriptors (by integrating their values over atomic or functional-group regions) gives
:ref:`condensed descriptors <condensed_tools>`, which identify the atoms, functional groups,
or bonds that are most reactive.
The remaining entries in the figure are non-local descriptors. **Non-local descriptors** provide information about how the
properties/reactivity of a molecule at one point change in response to changes (e.g., due to an attacking reagent)
elsewhere in the molecule. They can also be condensed, giving response matrices (2nd-order non-local
descriptors) and response tensors (high-order responses).

When computing conceptual DFT descriptors associated with electron transfer, one must choose an energy model
for the dependence of the energy upon the number of electrons. At the simplest level, one can choose to describe
these changes using the **frontier molecular orbital (FMO)** energies. Alternatively, one can compute the change in
energy due to electron donation/acceptance directly, as **finite differences (FD)** between the :math:`N-` electron
system’s energy and the energies of the :math:`\left(N-1\right)` and :math:`\left(N+1\right)` electron systems,

 .. math::
    IP = E\left(N - 1\right) - E\left(N\right) \approx -\varepsilon_{\text{HOMO}}  \\
    EA = E\left(N\right) - E\left(N + 1\right) \approx -\varepsilon_{\text{LUMO}}

where HOMO denotes the highest occupied molecular orbital and LUMO denotes the lowest unoccupied molecular orbital.

Once energies for the systems with integer electron number have been modeled using either the FMO or FD approaches,
a continuous model for the energy as a function of the number of electrons should be chosen.
Popular choices include the piecewise linear model, the quadratic model, the exponential model, and the rational
model. Of these, the piecewise linear model is the most mathematically rigorous and the quadratic model is the most
popular (and perhaps the most useful). The exponential and rational models have undesirable mathematical properties
and one should be especially wary about using them for local descriptors.

.. toctree::
   :maxdepth: 2

   sci_doc_globaltool
   sci_doc_localtool
   sci_doc_condensedtool






