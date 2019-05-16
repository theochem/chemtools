# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
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


__all__ = ["EigenValueTool"]


class EigenValueTool(object):
    r"""Class of descriptive tools based on eigenvalues."""

    def __init__(self, eigenvalues, eps=1e-15):
        r"""Initialize class.

        Parameters
        ----------
        eigenvalues : np.ndarray
            A two-dimensional array recording the eigenvalues of each point in a row.
        eps : float, optional
            The error bound for being a zero eigenvalue.

        """
        if not isinstance(eigenvalues, np.ndarray):
            raise TypeError("Eigenvalues should be a numpy array.")
        if not eigenvalues.ndim == 2:
            raise TypeError("Eigenvalues should be a two dimensional array.")
        self._eigenvalues = eigenvalues
        self._eps = eps

    @property
    def eigenvalues(self):
        r"""Eigenvalues of points."""
        return self._eigenvalues

    @property
    def ellipticity(self):
        r"""Ellipticity of electron density.

        .. math::
           \frac{\lambda_\text{min}}{\lambda_\text{min-1}} - 1
        """
        # get the two smallest eigenvalues
        index = np.argsort(self.eigenvalues, axis=1)[:, :2]
        min1 = self.eigenvalues[np.arange(self.eigenvalues.shape[0]), index[:, 0]]
        min2 = self.eigenvalues[np.arange(self.eigenvalues.shape[0]), index[:, 1]]
        # if np.abs(eigen2) < self._zero_eps:
        #    warnings.warn("Second largest eigenvalue is zero.")
        #     return np.inf
        return (min1 / min2) - 1.

    @property
    def bond_descriptor(self):
        r"""Bond descriptor defined as the ratio of average of positive and negative eigenvalues.

        .. math::
           \frac{\left(\frac{\sum_{\lambda_k > 0} \lambda_k}{\sum_{\lambda_k > 0} 1}\right)}
                {\left(\frac{\sum_{\lambda_k < 0} \lambda_k}{\sum_{\lambda_k < 0} 1}\right)}
        """
        # compute numerator
        pos_mask = (self.eigenvalues > self._eps).astype(int)
        result = np.sum(self.eigenvalues * pos_mask, axis=1) / np.sum(pos_mask, axis=1)
        # compute denominator
        neg_mask = (self.eigenvalues < -self._eps).astype(int)
        result /= (np.sum(self.eigenvalues * neg_mask, axis=1) / np.sum(neg_mask, axis=1))
        return result

    @property
    def eccentricity(self):
        r"""Eccentricity (essentially the condition number) of the set of eigenvalues.

        .. math ::
            \sqrt{\frac{\lambda_\text{max}}{\lambda_\text{min}}}
        """
        ratio = np.amax(self.eigenvalues, axis=1) / np.amin(self.eigenvalues, axis=1)
        # set negative values to None
        ratio[ratio < 0.] = np.nan
        return np.sqrt(ratio)

    @property
    def index(self):
        r"""Index of critical point which is the number of negative-curvature directions.

        .. math::
           \sum_{\lambda_k < 0} 1
        """
        return np.sum(self._eigenvalues < -self._eps, axis=1)

    @property
    def rank(self):
        r"""Rank of the critical point.

        The rank of a critical point is the number of positive eigenvalues.
        This is used to classify critical points on it's stability (trajectories going in or out).

        .. math::
            \sum_{\lambda_i > 0} 1
        """
        return np.sum(np.abs(self._eigenvalues) > self._eps, axis=1)

    @property
    def signature(self):
        r"""Signature of point(s) which is the difference of number of positive & negative values.

        .. math::
            \sum_{\lambda_k > 0.} 1 - \sum_{\lambda_k < 0.} 1
        """
        result = np.sum(self.eigenvalues > self._eps, axis=1)
        result -= np.sum(self.eigenvalues < -self._eps, axis=1)
        return result

    @property
    def morse(self):
        r"""Rank and signature of the critical point.

        .. math::
            \left(\sum_{\lambda_k > 0} 1, \sum_{\lambda_k > 0.}1 - \sum_{\lambda_k < 0.} 1\right)

        A system is degenerate if it has a zero eigenvalue and consequently, it's critical point
        is said to be "catastrophe". It returns a warning in this case.

        """
        if np.any(np.abs(self.eigenvalues) < self._eps):
            warnings.warn("Near catastrophic eigenvalue (close to zero) been found.")
        return np.hstack([self.rank[:, np.newaxis], self.signature[:, np.newaxis]])
