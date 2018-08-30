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
  * :class:`Condense Conceptual DFT Tools <toolbox.conceptual.CondensedConceptualDFT>`

* Density Based Tools

  * :class:`Noncovalent Interaction (NCI) <toolbox.nci.NCI>`

* Orbital Based Tools

  * :class:`Orbital Analysis <toolbox.orbitalbased.OrbitalAnalysis>`

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

  * :class:`Mixed Condensed Tool <conceptual.mixed.MixedCondensedTool>`


Density-Based Module
====================

* Density-Based Tools

  * :class:`Density Local Tool <toolbox.densitybased.DensityLocalTool>`

Orbital-Based Module
====================

* Orbital-Based Tools

  * :class:`Orbital Local Tool <toolbox.orbitalbased.OrbitalLocalTool>`

Utility Module
==============

* :func:`doc_inherit <utils.utils.doc_inherit>`
* :class:`BaseMolecule <utils.molecule.BaseMolecule>`
* :class:`HortonMolecule <utils.wrappers.HortonMolecule>`
* :class:`CubeGen <utils.cube.CubeGen>`

* Visualization Tools

  * VMD Scripts

    * :func:`print_vmd_script_nci <utils.output.print_vmd_script_nci>`
    * :func:`print_vmd_script_isosurface <utils.output.print_vmd_script_isosurface>`
    * :func:`print_vmd_script_multiple_cube <utils.output.print_vmd_script_multiple_cube>`
    * :func:`print_vmd_script_vector_field <utils.output.print_vmd_script_vector_field>`

.. Silent api generation
    .. autosummary::
      :toctree: modules/generated

      toolbox.conceptual.GlobalConceptualDFT
      toolbox.conceptual.LocalConceptualDFT
      toolbox.nci.NCI
      toolbox.orbitalbased.OrbitalAnalysis
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
      conceptual.mixed.MixedCondensedTool
      denstools.densitybased.DensityLocalTool
      orbtools.orbitalbased.OrbitalLocalTool
      utils.utils.doc_inherit
      utils.cube.CubeGen
      utils.molecule.BaseMolecule
      utils.wrappers.HortonMolecule
      utils.output.print_vmd_script_nci
      utils.output.print_vmd_script_isosurface
      utils.output.print_vmd_script_multiple_cube
      utils.output.print_vmd_script_vector_field

