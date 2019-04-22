.. _api:
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

  * :class:`Noncovalent Interaction (NCI) <toolbox.nci.NCI>`
  * :class:`Kinetic Energy Density <toolbox.kinetic.KineticEnergyDensity>`

* Orbital Based Tools

  * :class:`Orbital Local Tool <toolbox.orbbased.OrbitalLocalTool>`


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

  * :class:`Density Local Tool <denstools.densitybased.DensityLocalTool>`


Wrappers Module
===============

* :class:`Molecule <wrappers.molecule.Molecule>`


Utility Module
==============

* :func:`doc_inherit <utils.utils.doc_inherit>`
* :class:`CubeGen <utils.cube.CubeGen>`


* Visualization Tools

  * VMD Scripts

    * :func:`print_vmd_script_nci <outputs.output_vmd.print_vmd_script_nci>`
    * :func:`print_vmd_script_isosurface <outputs.output_vmd.print_vmd_script_isosurface>`
    * :func:`print_vmd_script_multiple_cube <outputs.output_vmd.print_vmd_script_multiple_cube>`
    * :func:`print_vmd_script_vector_field <outputs.output_vmd.print_vmd_script_vector_field>`

.. Silent api generation
    .. autosummary::
      :toctree: modules/generated

      toolbox.conceptual.GlobalConceptualDFT
      toolbox.conceptual.LocalConceptualDFT
      toolbox.conceptual.CondensedConceptualDFT
      toolbox.nci.NCI
      toolbox.kinetic.KineticEnergyDensity
      toolbox.orbbased.OrbitalLocalTool
      denstools.densitybased.DensityLocalTool
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
      wrappers.molecule.Molecule
      outputs.output_vmd.print_vmd_script_nci
      outputs.output_vmd.print_vmd_script_isosurface
      outputs.output_vmd.print_vmd_script_multiple_cube
      outputs.output_vmd.print_vmd_script_vector_field
      utils.utils.doc_inherit
      utils.cube.CubeGen

