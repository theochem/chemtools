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
'''Density-Based Local Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


import numpy as np


class DensityLocalTool(object):
    '''
    Class of desnity-based local descriptive tools.
    '''
    def __init__(self, density, gradient=None, laplacian=None):
        '''
        Parameters
        ----------
        density : np.ndarray
            Density of the system evaluated on a grid
        gradient : np.ndarray, default=None
            Gradient of the density evaluated on a grid
        laplacian : np.ndarray, default=None
            Laplacian density evaluated on a grid
        '''
        self._density = density
        self._gradient = gradient
        self._laplacian = laplacian

    @property
    def density(self):
        '''
        '''
        return self._density

    @property
    def gradient(self):
        '''
        '''
        return self._gradient

    @property
    def laplacian(self):
        '''
        '''
        return self._laplacian

    @property
    def shanon_information(self):
        r'''
        Shanon information defined as :math:`\rho(r) \ln \rho(r)`.
        '''
        # masking might be needed
        value = self._density * np.log(self._density)
        return value
