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
'''Local Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


from horton import doc_inherit


class BaseLocalTool(object):
    '''
    Base Class of Local Conceptual DFT Reactivity Descriptors.
    '''
    def __init__(self, ff_plus=None, ff_minus=None, global_instance=None):
        '''
        Parameters
        ----------
        ff_plus : np.ndarray
            Positive Fukui Function.
        ff_minus : np.ndarray
            Negative Fukui Function.
        global_instance :
            Instance of ``GlobalConceptualTool`` class.
        '''
        self._ff_plus = ff_plus
        self._ff_minus = ff_minus
        self._global_instance = global_instance

        # Define zero Fukui Functional as the average of FF+ and FF-
        if (self._ff_plus is not None) and (self._ff_minus is not None):
            self._ff_zero = 0.5 * (self._ff_plus + self._ff_minus)
        else:
            self._ff_zero = None

    @property
    def ff_plus(self):
        '''Fukui Function from above (positive Fukui Function).'''
        return self._ff_plus

    @property
    def ff_minus(self):
        '''Fukui Function from below (negative Fukui Function).'''
        return self._ff_minus

    @property
    def ff_zero(self):
        '''Fukui Function from center.'''
        return self._ff_zero

    @property
    def dual_descriptor(self):
        '''Dual Descriptor.'''
        if (self._ff_plus is not None) and (self._ff_minus is not None):
            value = self._ff_plus - self._ff_minus
            return value

    def __getattr__(self, attr):
        '''
        '''
        # Identify the global property and the type of Fukui Function
        global_prop, ff_type = attr.rsplit('_', 1)

        # Check for availability of GlobalConceptualTool instance
        if self._global_instance is None:
            raise ValueError('The argument global_instance is None!')

        # Check for availability of global property
        if global_prop not in dir(self._global_instance):
            raise ValueError('The global property={0} is not known.'.format(global_prop))

        # Check for availability of the type of Fukui Function
        if ff_type not in ['plus', 'minus', 'zero']:
            raise ValueError('The attribute ff_{0} is not known.'.format(ff_type))

        # Get the global property & the Fukui Function
        global_descriptor = getattr(self._global_instance, global_prop)
        fukui_function = getattr(self, 'ff_' + ff_type)
        if fukui_function is None:
            raise ValueError('The ff_{0} is None!'.format(ff_type))

        # Compute the local property
        local_descriptor = global_descriptor * fukui_function

        return local_descriptor


class LinearLocalTool(BaseLocalTool):
    '''
    Class of Local Conceptual DFT Reactivity Descriptors based on the Linear Energy Model.
    '''

    @doc_inherit(BaseLocalTool)
    def __init__(self, ff_plus=None, ff_minus=None, global_instance=None):
        super(self.__class__, self).__init__(ff_plus, ff_minus, global_instance)

    @property
    def mu_plus(self):
        '''The local plus chemical potential in the linear model is the positive Fukui function'''
        return self._ff_plus

    @property
    def mu_minus(self):
        '''The local minus chemical potential in the linear model is the negative Fukui function'''
        return self._ff_minus

    @property
    def mu_zero(self):
        '''The local zero chemical potential in the linear model is the zero Fukui function'''
        return self._ff_zero


class QuadraticLocalTool(BaseLocalTool):
    '''
    Class of Local Conceptual DFT Reactivity Descriptors based on the Quadratic Energy Model.
    '''

    @doc_inherit(BaseLocalTool)
    def __init__(self, ff_plus=None, ff_minus=None, global_instance=None):
        super(self.__class__, self).__init__(ff_plus, ff_minus, global_instance)
