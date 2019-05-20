.. _api:
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

*****************
API Documentation
*****************

.. module:: chemtools


Toolbox Module
==============

* Conceptual DFT Tools

  * :class:`Global Conceptual DFT Tools <toolbox.conceptual.GlobalConceptualDFT>`
  * :class:`Local Conceptual DFT Tools <toolbox.conceptual.LocalConceptualDFT>`
  * :class:`Condensed Conceptual DFT Tools <toolbox.conceptual.CondensedConceptualDFT>`

* Density Based Tools

  * :class:`Density Local Tool <toolbox.densbased.DensityLocalTool>`
  * :class:`Noncovalent Interaction (NCI) <toolbox.interactions.NCI>`
  * :class:`Electron Localization Function (ELF) <toolbox.interactions.ELF>`
  * :class:`Localized orbital Locator (LOL) <toolbox.interactions.LOL>`
  * :class:`Kinetic Energy Density (KED) <toolbox.kinetic.KED>`

* Orbital Based Tools

  * :class:`Orbital Local Tool <toolbox.orbsbased.OrbitalLocalTool>`


Molecular Orbital (MO) Theory Module
====================================

* :class:`Molecular Orbtital Theory Based Tool <toolbox.motbased.MOTBasedTool>`



Conceptual Module
=================

* Global Conceptual DFT Tools

  * :class:`Base Global Tool <conceptual.base.BaseGlobalTool>`
  * :class:`Linear Global Tool <conceptual.linear.LinearGlobalTool>`
  * :class:`Quadratic Global Tool <conceptual.quadratic.QuadraticGlobalTool>`
  * :class:`Exponential Global Tool <conceptual.exponential.ExponentialGlobalTool>`
  * :class:`Rational Global Tool <conceptual.rational.RationalGlobalTool>`
  * :class:`Cubic Global Tool <conceptual.cubic.CubicGlobalTool>`
  * :class:`General Global Tool <conceptual.general.GeneralGlobalTool>`
  * :class:`Mixed Global Tool <conceptual.mixed.MixedGlobalTool>`

* Local Conceptual DFT Tools

  * :class:`Base Local Tool <conceptual.base.BaseLocalTool>`
  * :class:`Linear Local Tool <conceptual.linear.LinearLocalTool>`
  * :class:`Quadratic Local Tool <conceptual.quadratic.QuadraticLocalTool>`
  * :class:`Mixed Local Tool <conceptual.mixed.MixedLocalTool>`

* Condensed Conceptual DFT Tools

  * :class:`Base Condensed Tool <conceptual.base.BaseCondensedTool>`
  * :class:`Linear Condensed Tool <conceptual.linear.LinearCondensedTool>`
  * :class:`Quadratic Condensed Tool <conceptual.quadratic.QuadraticCondensedTool>`
  * :class:`Mixed Condensed Tool <conceptual.mixed.MixedCondensedTool>`


Density-Based Module
====================

* Density-Based Tools

  * :class:`Density Based Tool <denstools.densbased.DensTool>`
  * :class:`Density & Gradient Based Tool <denstools.densbased.DensGradTool>`
  * :class:`Density, Gradient & Laplacian Based Tool <denstools.densbased.DensGradLapTool>`
  * :class:`Density, Gradient, Laplacian & KED Based Tool <denstools.densbased.DensGradLapKedTool>`


Topological Analysis
====================

  * :class:`Eigenvalue Descriptors <topology.eigenvalues.EigenValueTool>`


Wrappers Module
===============

* :class:`Molecule <wrappers.molecule.Molecule>`
* :class:`MolecularGrid <wrappers.grid.MolecularGrid>`
* :class:`UniformGrid <utils.cube.UniformGrid>`


Output Module
=============

* VMD Scripts

  * :func:`print_vmd_script_nci <outputs.vmd.print_vmd_script_nci>`
  * :func:`print_vmd_script_isosurface <outputs.vmd.print_vmd_script_isosurface>`
  * :func:`print_vmd_script_multiple_cube <outputs.vmd.print_vmd_script_multiple_cube>`
  * :func:`print_vmd_script_vector_field <outputs.vmd.print_vmd_script_vector_field>`

* 2-D Plots

  * :func:`plot_scatter <outputs.plot.plot_scatter>`


Utilities
=========

* :func:`plane_mesh <utils.mesh.plane_mesh>`



.. Silent api generation
    .. autosummary::
      :toctree: modules/generated

      toolbox.conceptual.GlobalConceptualDFT
      toolbox.conceptual.LocalConceptualDFT
      toolbox.conceptual.CondensedConceptualDFT
      toolbox.densbased.DensityLocalTool
      toolbox.motbased.MOTBasedTool
      toolbox.interactions.NCI
      toolbox.interactions.ELF
      toolbox.interactions.LOL
      toolbox.kinetic.KED
      toolbox.orbsbased.OrbitalLocalTool
      denstools.densbased.DensTool
      denstools.densbased.DensGradTool
      denstools.densbased.DensGradLapTool
      denstools.densbased.DensGradLapKedTool
      conceptual.base.BaseGlobalTool
      conceptual.linear.LinearGlobalTool
      conceptual.quadratic.QuadraticGlobalTool
      conceptual.exponential.ExponentialGlobalTool
      conceptual.rational.RationalGlobalTool
      conceptual.cubic.CubicGlobalTool
      conceptual.general.GeneralGlobalTool
      conceptual.mixed.MixedGlobalTool
      conceptual.base.BaseLocalTool
      conceptual.linear.LinearLocalTool
      conceptual.quadratic.QuadraticLocalTool
      conceptual.mixed.MixedLocalTool
      conceptual.base.BaseCondensedTool
      conceptual.linear.LinearCondensedTool
      conceptual.quadratic.QuadraticCondensedTool
      conceptual.mixed.MixedCondensedTool
      topology.eigenvalues.EigenValueTool
      wrappers.molecule.Molecule
      wrappers.grid.MolecularGrid
      outputs.vmd.print_vmd_script_nci
      outputs.vmd.print_vmd_script_isosurface
      outputs.vmd.print_vmd_script_multiple_cube
      outputs.vmd.print_vmd_script_vector_field
      outputs.plot.plot_scatter
      utils.cube.UniformGrid
      utils.mesh.plane_mesh

