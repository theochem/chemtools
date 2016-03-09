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


.. _orbital_tools:

Orbital-Based Local Descriptors :class:`chemtools.tool.orbitaltool`
###################################################################

All the tools for calculating which the orbital information of the :math:`N` electron reference state is enough.

Positive Definite Kinetic Energy Density: :class:`chemtools.tool.orbitaltool.kinetic_energy_density`
----------------------------------------------------------------------------------------------------

.. math:: 

    \tau_{\sigma} (\mathbf{r}) = 
        \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2 


Electron Localisation Function: :class:`chemtools.tool.orbitaltool.elf`
-----------------------------------------------------------------------

The concepts of chemical bonds and electron pairs are fundamental in our 
understanding of chemistry. Quantum chemical calculations, however, 
yield canonical orbitals which are delocalized in space and do not contain this 
*chemical* information at first sight. There are different schemes to generate
*localized* orbitals by means of an unitary transformation of the canonical orbitals,
but none of these schemes are unique.

The Electron Localization Function (ELF), introduced by A. D. Becke and K. E. Edgecombe, 
was introduced as an alternative to these unitary transformations.
It measures the likelihood of finding an electron with the same spin 
in the neighbourhood of a reference electron.
The derivation of ELF starts from the conditional probability 
:math:`P^{\sigma \sigma}_{cond} (\mathbf{r}, \mathbf{r}')` 
to find an electron with :math:`\sigma` spin at position :math:`\mathbf{r}'` 
when a reference electron with the same spin is found at :math:`\mathbf{r}`.
A. D. Becke showed that, when :math:`\mathbf{r}'` approaches :math:`\mathbf{r}`, 
the leading term in the *Taylor expansion* of the *spherically averaged* pair probability 
is:

 .. math:: 

    P^{\sigma \sigma}_{cond} (\mathbf{r},s) = \frac{1}{3} 
        \lbrack \tau_{\sigma} (\mathbf{r}) - 
	\frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}} \rbrack s^2

where :math:`(\mathbf{r},s)` denotes the average on a sphere of radius :math:`s` 
around the reference electron at :math:`(\mathbf{r}` and :math:`\tau_{\sigma}` is the positive definite kinetic energy density:

.. math:: 

    \tau_{\sigma} (\mathbf{r}) = 
        \sum_i^{\sigma} \lvert \nabla \phi_i (\mathbf{r}) \rvert^2 

The smaller the probability of finding a second electron with the same spin 
near the reference electron, the more localized the reference electron, 
so the information on the electron localisation can thus be found by:

 .. math:: 

    D_{\sigma} (\mathbf{r}) =  \tau_{\sigma} (\mathbf{r}) - 
	\frac{1}{4} \frac{(\nabla \rho_{\sigma})^2}{\rho_{\sigma}} .

There is an inverse relation between the electron localisation and :math:`D_{\sigma}`: 
highly localized electrons can be found in regions of small :math:`D_{\sigma}`. 
Moreover, :math:`D_{\sigma}` is not bounded from above. 
For these reasons A. D. Becke and K. E. Edgecombe proposed the following formula
for the Electron Localization Function:

 .. math:: 

    ELF (\mathbf{r}) = 
        \frac{1}{\left( 1 + \left(\frac{D_{\sigma}(\mathbf{r})}
        {D_{\sigma}^0 (\mathbf{r})} \right)^2\right)} , 

where :math:`D_{\sigma}^0` corresponds to the *uniform electron gas*:

 .. math:: 

    D_{\sigma}^0 (\mathbf{r}) =  
        \frac{3}{5} (6 \pi^2)^{2/3} \rho_{\sigma}^{5/3} (\mathbf{r}) ,

making the ratio of :math:`D_{\sigma} / D_{\sigma}^0` dimensionless and scaled with 
respect to the uniform electron gas. As such, the ELF is always in the range

 .. math:: 

    0 \leq ELF \leq 1 , 

with ELF = 1 corresponding to the perfectly localised reference electron and 
ELF = :math:`^1/_2` to the uniform electron gas-like correlation.


**References:**
  * `Becke A. D. and Edgecombe K. E., J. Chem. Phys. (1990), 92, 5397 <http://scitation.aip.org/content/aip/journal/jcp/92/9/10.1063/1.458517>`_.
  * `Savin A., Becke A. D., Flad J., Nesper R., Preuss H. and von Schnering H. G., Angew. Chem. Int. Ed. Engl. (1991), 30, 409 <http://onlinelibrary.wiley.com/doi/10.1002/anie.199104091/full>`_.
