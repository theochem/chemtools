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

Analysis Tools
==============

* Conceptual DFT Tools

  * :class:`Global Conceptual DFT Tools <analysis.conceptual.GlobalConceptualDFT>`
  * :class:`Local Conceptual DFT Tools <analysis.conceptual.LocalConceptualDFT>`

* Other Tools

  * :class:`Noncovalent Interaction <analysis.others.NCI>`
  * :func:`_compute_hessian <analysis.others._compute_hessian>`

* Visualization Tools

  * VMD Output

    * :func:`_print_vmd_script_nci <analysis.output._print_vmd_script_nci>`
    * :func:`print_vmd_script_isosurface <analysis.output.print_vmd_script_isosurface>`

Conceptual Tools
================

* Global Conceptual DFT Tools

  * :class:`Base Global Tool <tool.globaltool.BaseGlobalTool>`
  * :class:`Linear Global Tool <tool.globaltool.LinearGlobalTool>`
  * :class:`Quadratic Global Tool <tool.globaltool.QuadraticGlobalTool>`
  * :class:`Exponential Global Tool <tool.globaltool.ExponentialGlobalTool>`
  * :class:`Rational Global Tool <tool.globaltool.RationalGlobalTool>`
  * :class:`General Global Tool <tool.globaltool.GeneralGlobalTool>`

* Local Conceptual DFT Tools

  * :class:`Base Local Tool <tool.localtool.BaseLocalTool>`
  * :class:`Linear Local Tool <tool.localtool.LinearLocalTool>`
  * :class:`Quadratic Local Tool <tool.localtool.QuadraticLocalTool>`
  * :class:`Density Local Tool <tool.densitytool.DensityLocalTool>`
  * :class:`Orbital Local Tool <tool.orbitaltool.OrbitalLocalTool>`

* :class:`Condensed Conceptual DFT Tool <tool.condensedtool.CondensedTool>`
 
Utility Tools
=============
* :func:`doc_inherit <utils.doc_inherit>`
* :class:`CubeGen <utils.CubeGen>`


.. Silent api generation
    .. autosummary::
      :toctree: modules/generated

      analysis.conceptual.GlobalConceptualDFT
      analysis.conceptual.LocalConceptualDFT
      analysis.others.NCI
      analysis.others._compute_hessian
      analysis.output._print_vmd_script_nci
      analysis.output.print_vmd_script_isosurface
      tool.globaltool.BaseGlobalTool
      tool.globaltool.LinearGlobalTool
      tool.globaltool.QuadraticGlobalTool
      tool.globaltool.ExponentialGlobalTool
      tool.globaltool.RationalGlobalTool
      tool.globaltool.GeneralGlobalTool
      tool.localtool.BaseLocalTool
      tool.localtool.LinearLocalTool
      tool.localtool.QuadraticLocalTool
      tool.densitytool.DensityLocalTool
      tool.orbitaltool.OrbitalLocalTool
      tool.condensedtool.CondensedTool
      utils.doc_inherit
      utils.CubeGen

