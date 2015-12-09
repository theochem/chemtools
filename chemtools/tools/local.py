#!/usr/bin/env python
'''Local Conceptual Density Functional Theory (DFT) Reactivity Tools.'''


class QuadraticLocalTool(object):
    '''
    Class of Local Conceptual DFT Reactivity Descriptors based on the Quadratic Energy Model.
    '''
    def __init__(self, ff_above=None, ff_below=None, global_instance=None):
        '''
        Parameters
        ----------
        ff_above : np.ndarray
            Positive Fukui Function.
        ff_below : np.ndarray
            Negative Fukui Function.
        global_instance :
            Instance of ``GlobalConceptualTool`` class.
        '''
        self._ff_above = ff_above
        self._ff_below = ff_below
        self._global_instance = global_instance

        # Define central Fukui Functional as the average of FF+ and FF-
        if (self._ff_above is not None) and (self._ff_below is not None):
            self._ff_central = 0.5 * (self._ff_above + self._ff_below)
        else:
            self._ff_central = None

    @property
    def ff_above(self):
        '''Fukui Function from above (positive Fukui Function).'''
        return self._ff_above

    @property
    def ff_below(self):
        '''Fukui Function from below (negative Fukui Function).'''
        return self._ff_below

    @property
    def ff_central(self):
        '''Fukui Function from center.'''
        return self._ff_central

    @property
    def dual_descriptor(self):
        '''Dual Descriptor.'''
        if (self._ff_above is not None) and (self._ff_below is not None):
            value = self._ff_above - self._ff_below
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
