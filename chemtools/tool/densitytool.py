
'''Density-Based Local Conceptual Density Functional Theory (DFT) Reactivity Tools.

'''

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
