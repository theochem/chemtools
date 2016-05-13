# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2014-2015 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
'''Condensed Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


import numpy as np
from chemtools.utils import doc_inherit


class CondensedTool(object):
    '''
    Class of condensed conceptual DFT reactivity descriptors.

    So far only the Fragment of Molecular Response is used,
    where the weights do not depend on the number of electrons.
    '''
    def __init__(self, dens_part):
        '''
        Parameters
        ----------
        dens_part:
            A WPartClass object obtained from partitioning the molecule
        '''
        self._dens_part = dens_part

    def condense_atoms(self, local_property):
        r'''
        Return condensed values of the local descriptor :math:`p_{\text{local}}\left(\mathbf{r}\right)`
        into atomic contribution :math:`P_A` defined as:

        .. math::

           P_A = \int \omega_A\left(\mathbf{r}\right) p_{\text{local}}\left(\mathbf{r}\right) d\mathbf{r}

        Parameters
        ----------
        local_property : np.ndarray
            Local descriptor evaluated on grid.
        '''
        natom = self._dens_part.natom
        local_condensed = np.zeros(natom)
        for index in xrange(natom):
            at_grid = self._dens_part.get_grid(index)
            at_weight = self._dens_part.cache.load('at_weights',index)
            wcor = self._dens_part.get_wcor(index)
            local_prop = self._dens_part.to_atomic_grid(index, local_property)
            local_condensed[index] = at_grid.integrate(at_weight, local_prop, wcor)
        return local_condensed

    def condese_pairs(self, response):
        r'''
        Return condensed values of the response function :math:`f_{\text{response}}\left(\mathbf{r}, \mathbf{r'}\right)`
        into atomic pair contribution :math:`P_{AB}` defined as:

        .. math::

           P_{AB} = \int \omega_A\left(\mathbf{r}\right) \omega_B\left(\mathbf{r'}\right)
                       f_{\text{response}}\left(\mathbf{r}, \mathbf{r'}\right) d\mathbf{r} d\mathbf{r'}

        Parameters
        ----------
        response : np.ndarray
            Response evaluated on grid.
        '''
        pass
