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


class BaseCondensedTool(object):
    '''
    Base class of condensed conceptual DFT reactivity descriptors.
    '''
    def __init__(self, weights, grid):
        '''
        Parameters
        ----------
        weights :
            List of atomic weights
        grid :
            Molecular grid
        '''
        self._weights = weights
        self._grid = grid

    @property
    def weights(self):
        '''
        '''
        return self._weights

    @property
    def grid(self):
        '''
        '''
        return self._grid

    def condense_atoms(self, local_property):
        r'''
        Condense local descriptor :math:`p_{\text{local}}\left(\mathbf{r}\right)` into
        atomic contribution :math:`P_A` defined as:

        .. math::

           P_A = \int \omega_A\left(\mathbf{r}\right) p_{\text{local}}\left(\mathbf{r}\right) d\mathbf{r}

        Parameters
        ----------
        local_property : np.ndarray
            Local descriptor evaluated on grid.
        '''
        pass

    def condese_pairs(self, response):
        r'''
        Condense response function :math:`f_{\text{response}}\left(\mathbf{r}, \mathbf{r'}\right)` into
        atomic pair contribution :math:`P_{AB}` defined as:

        .. math::

           P_{AB} = \int \omega_A\left(\mathbf{r}\right) \omega_B\left(\mathbf{r'}\right)
                       f_{\text{response}}\left(\mathbf{r}, \mathbf{r'}\right) d\mathbf{r} d\mathbf{r'}

        Parameters
        ----------
        response : np.ndarray
            Response evaluated on grid.
        '''
        pass
