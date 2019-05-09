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
r"""Contains class responsible for calculating various descriptors based on eigenvalues."""


import warnings
import numpy as np


__all__ = ["EigenDescriptor"]


class EigenDescriptor(object):
    r"""Class of descriptive tools based on eigenvalues."""

    def __init__(self, eigenvals, zero_eps=1e-15):
        r"""Initialize class.

        Parameters
        ----------
        eigenvals : np.ndarray
            A two-dimensional array holding the eigenvalues separately for each point.
        zero_eps : float, optional
            The error bound for being a zero eigenvalue.

        """
        if not isinstance(eigenvals, np.ndarray):
            raise TypeError("Eigenvalues should be a numpy array.")
        if not eigenvals.ndim == 2:
            raise TypeError("Eigenvalues should be a two dimensional array.")
        self._eigenvals = eigenvals
        self._zero_eps = zero_eps

    @property
    def eigenvals(self):
        r"""Set of two-dimensional eigenvalues."""
        return self._eigenvals

    @property
    def ellipticity(self):
        r"""Ellipticity of a bond critical point.

        .. math::
           \frac{\lambda_\text{max}}{\lambda_\text{max-1}} - 1

        Returns
        -------
        float :
            Returns the ratio of the largest two eigenvalues minus one. If second largest eigenvalue
            is zero, then infinity is returned.
        """
        # get the two largest eigenvalues
        index = np.argsort(self.eigenvals, axis=1)[:, -2:]
        max1 = self.eigenvals[np.arange(self.eigenvals.shape[0]), index[:, 1]]
        max2 = self.eigenvals[np.arange(self.eigenvals.shape[0]), index[:, 0]]
        # if np.abs(eigen2) < self._zero_eps:
        #    warnings.warn("Second largest eigenvalue is zero.")
        #     return np.inf
        return (max1 / max2) - 1.

    @property
    def bond_descriptor(self):
        r"""Bond descriptor defined as the ratio of average of positive and negative eigenvalues.

        .. math::
           \frac{\left(\frac{\sum_{\lambda_k > 0} \lambda_k}{\sum_{\lambda_k > 0} 1}\right)}
                {\left(\frac{\sum_{\lambda_k < 0} \lambda_k}{\sum_{\lambda_k < 0} 1}\right)}

        Returns
        -------
        float :
            Ratio of average positive eigenvalues to average negative eigenvalues. If there are no
            negative eigenvalues, then None is returned.
        """
        # compute numerator
        pos_mask = (self._eigenvals > self._zero_eps).astype(int)
        result = np.sum(self._eigenvals * pos_mask, axis=1)
        result /= np.sum(pos_mask, axis=1)
        # compute denominator
        neg_mask = (self._eigenvals < -self._zero_eps).astype(int)
        result /= (np.sum(self._eigenvals * neg_mask, axis=1) / np.sum(neg_mask, axis=1))
        return result

    @property
    def eccentricity(self):
        r"""Eccentricity (essentially the condition number) of the set of eigenvalues.

        .. math ::
            \sqrt{\frac{\lambda_\text{max}}{\lambda_\text{min}}}

        Returns
        -------
        float :
            The condition number, the square root of largest eigenval divided by minimum eigenval.
            If one of maximima or mininum is negative, then none is returned.
        """
        ratio = np.amax(self._eigenvals, axis=1) / np.amin(self._eigenvals, axis=1)
        # set negative values to None
        ratio[ratio < 0.] = np.nan
        return np.sqrt(ratio)

    @property
    def index(self):
        r"""Index of critical point which is the number of negative-curvature directions.

        .. math::
           \sum_{\lambda_k < 0} 1

        Returns
        -------
        int :
            Number of negative eigenvalues.
        """
        return np.sum(self._eigenvals < -self._zero_eps, axis=1)

    @property
    def rank(self):
        r"""Rank of the critical point.

        The rank of a critical point is the number of positive eigenvalues.
        This is used to classify critical points on it's stability (trajectories going in or out).

        .. math::
            \sum_{\lambda_i > 0} 1

        Returns
        -------
        int :
            The number of non-zero eigenvalues.
        """
        return np.sum(np.abs(self._eigenvals) > self._zero_eps, axis=1)

    @property
    def signature(self):
        r"""Signature of the critical point.

        This is used to classify saddle critical points (ie certain directions flow outwards and
        certain directions flow inwards).

        .. math::
            \sum_{\lambda_k > 0.} 1 - \sum_{\lambda_k < 0.} 1

        Returns
        -------
        int :
            The number of positive eigenvalues minus the number of negative eigenvalues.
        """
        result = np.sum(self._eigenvals > self._zero_eps, axis=1)
        result -= np.sum(self._eigenvals < -self._zero_eps, axis=1)
        return result

    @property
    def morse(self):
        r"""Rank and signature of the critical point.

        .. math::
            \left(\sum_{\lambda_k > 0} 1, \sum_{\lambda_k > 0.}1 - \sum_{\lambda_k < 0.} 1\right)

        A system is degenerate if it has a zero eigenvalue and consequently, it's critical point
        is said to be "catastrophe". It returns a warning in this case.

        Returns
        -------
        (int, int) :
            Returns the rank and signature of the critical point.
        """
        if np.any(np.abs(self._eigenvals) < self._zero_eps):
            warnings.warn("Near catastrophic eigenvalue (close to zero) been found.")
        return np.hstack([self.rank[:, np.newaxis], self.signature[:, np.newaxis]])
