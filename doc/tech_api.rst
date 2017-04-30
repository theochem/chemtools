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

Analysis Module
===============

* Conceptual DFT Tools

  * :class:`Global Conceptual DFT Tools <analysis.conceptual.GlobalConceptualDFT>`
  * :class:`Local Conceptual DFT Tools <analysis.conceptual.LocalConceptualDFT>`
  * :class:`Condense Conceptual DFT Tools <analysis.conceptual.CondensedConceptualDFT>`

* Density Based Tools

  * :class:`Noncovalent Interaction (NCI) <analysis.densitybased.NCI>`
  * :func:`_compute_hessian <analysis.densitybased._compute_hessian>`

* Orbital Based Tools

  * :class:`Orbital Analysis <analysis.orbitalbased.OrbitalAnalysis>`

ToolBox Module
==============

* Global Conceptual DFT Tools

  * :class:`Base Global Tool <toolbox.conceptualglobal.BaseGlobalTool>`
  * :class:`Linear Global Tool <toolbox.conceptualglobal.LinearGlobalTool>`
  * :class:`Quadratic Global Tool <toolbox.conceptualglobal.QuadraticGlobalTool>`
  * :class:`Exponential Global Tool <toolbox.conceptualglobal.ExponentialGlobalTool>`
  * :class:`Rational Global Tool <toolbox.conceptualglobal.RationalGlobalTool>`
  * :class:`General Global Tool <toolbox.conceptualglobal.GeneralGlobalTool>`

* Local Conceptual DFT Tools

  * :class:`Base Local Tool <toolbox.conceptuallocal.BaseLocalTool>`
  * :class:`Linear Local Tool <toolbox.conceptuallocal.LinearLocalTool>`
  * :class:`Quadratic Local Tool <toolbox.conceptuallocal.QuadraticLocalTool>`

* Condensed Conceptual DFT Tools

  * :class:`Condensed Conceptual DFT Tool <toolbox.conceptualcondense.CondensedTool>`

* Density Based Tools

  * :class:`Density Local Tool <toolbox.densitybased.DensityLocalTool>`

* Orbital Based Tools

  * :class:`Orbital Local Tool <toolbox.orbitalbased.OrbitalLocalTool>`

Utility Tools
=============

* :func:`doc_inherit <utils.utils.doc_inherit>`
* :class:`Molecule <utils.molecule.Molecule>`
* :class:`WaveFunctions <utils.molecule.WaveFunction>`
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

      analysis.conceptual.GlobalConceptualDFT
      analysis.conceptual.LocalConceptualDFT
      analysis.densitybased.NCI
      analysis.densitybased._compute_hessian
      analysis.orbitalbased.OrbitalAnalysis
      toolbox.conceptualglobal.BaseGlobalTool
      toolbox.conceptualglobal.LinearGlobalTool
      toolbox.conceptualglobal.QuadraticGlobalTool
      toolbox.conceptualglobal.ExponentialGlobalTool
      toolbox.conceptualglobal.RationalGlobalTool
      toolbox.conceptualglobal.GeneralGlobalTool
      toolbox.conceptuallocal.BaseLocalTool
      toolbox.conceptuallocal.LinearLocalTool
      toolbox.conceptuallocal.QuadraticLocalTool
      toolbox.conceptualcondense.CondensedTool
      toolbox.densitybased.DensityLocalTool
      toolbox.orbitalbased.OrbitalLocalTool
      utils.utils.doc_inherit
      utils.cube.CubeGen
      utils.molecule.Molecule
      utils.molecule.WaveFunction
      utils.output.print_vmd_script_nci
      utils.output.print_vmd_script_isosurface
      utils.output.print_vmd_script_multiple_cube
      utils.output.print_vmd_script_vector_field

