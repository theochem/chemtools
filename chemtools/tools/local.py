#!/usr/bin/env python
'''Local Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


class QuadraticLocalTool(object):
    '''
    Class of Local Conceptual DFT Reactivity Descriptors based on the Quadratic Energy Model.
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

        # Define central Fukui Functional as the average of FF+ and FF-
        if (self._ff_plus is not None) and (self._ff_minus is not None):
            self._ff_central = 0.5 * (self._ff_plus + self._ff_minus)
        else:
            self._ff_central = None

    @property
    def ff_plus(self):
        '''Fukui Function from above (positive Fukui Function).'''
        return self._ff_plus

    @property
    def ff_minus(self):
        '''Fukui Function from below (negative Fukui Function).'''
        return self._ff_minus

    @property
    def ff_central(self):
        '''Fukui Function from center.'''
        return self._ff_central

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
        if ff_type not in ['plus', 'minus', 'central']:
            raise ValueError('The attribute ff_{0} is not known.'.format(ff_type))

        # Get the global property & the Fukui Function
        global_descriptor = getattr(self._global_instance, global_prop)
        fukui_funcion = getattr(self, 'ff_' + ff_type)
        if fukui_funcion is None:
            raise ValueError('The ff_{0} is None!'.format(ff_type))

        # Compute the local property
        local_descriptor = global_descriptor * fukui_funcion

        return local_descriptor
