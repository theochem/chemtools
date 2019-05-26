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


.. _dftbased_esp:

Electrostatic Potential (ESP)
=============================


The electrostatic potential (ESP) is used to analyze the positively and negatively charged regions
of a molecule. In atomic units, it is defined as the interaction energy of the molecule with an
infinitesimal positive point charge, per unit charge:

 .. math::
    \Phi \left(\mathbf{r}\right) =
    \sum_{A=1}^{N_\text{atoms}} \frac{Z_A}{\rvert \mathbf{R}_A - \mathbf{r} \lvert} -
    \int \frac{\rho \left(\mathbf{r^\prime}\right)}{\rvert \mathbf{r} - \mathbf{r^\prime} \lvert}
    d\mathbf{r^\prime}

Here, :math:`Z_A` is the charge on an atomic nucleus, :math:`\mathbf{R}_A` is the coordinates of an
atomic nucleus, and :math:`\rho(\mathbf{r})` is the electron density. The ESP can also be written
in terms of the external potential, :math:`\mathit{v}(\mathbf{r})`, used in DFT:

 .. math::
    \Phi \left(\mathbf{r}\right) = - \left(\mathit{v}(\mathbf{r}) +
    \int \frac{\rho \left(\mathbf{r^\prime}\right)}{\rvert \mathbf{r} - \mathbf{r^\prime} \lvert}
    d\mathbf{r^\prime}\right)

Typically, hard acids/electrophiles attack a molecule where the electrostatic potential is most
negative, and the hard bases/nucleophiles attack a molecule where it is most positive.
To predict how a molecule would interact with an approaching reagent, one commonly plots the
electrostatic potential on a van der Waals surface of the molecule or an iso-electron-density
surface of about :math:`\rho(\mathbf{r})=0.002`.

 .. TODO::
    Add links to example gallery.

The value of the electrostatic potential at a nucleus is often interesting, because it shows how
the energy changes when the atomic number of the nucleus changes (to first order).
Therefore it is relevant for alchemical changes as an atom changes to an adjacent atom in the
periodic table,

 .. math::
    E(Z_A \pm 1) - E(Z_A) \stackrel{\text{1st order}}{\approx} \pm \frac{\partial E}{\partial Z_A}
    = \pm \left(
    \sum_{B=1 \\ B \neq A}^{N_\text{atoms}} \frac{Z_B}{\rvert \mathbf{R}_B - \mathbf{r} \lvert} -
    \int \frac{\rho \left(\mathbf{r^\prime}\right)}{\rvert \mathbf{r} - \mathbf{r^\prime} \lvert}
    d\mathbf{r^\prime}\right)

For example, at a hydrogen nucleus, the electrostatic potential correlates tightly with the proton
affinity (ergo pKa).
