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
    r"""
    Stores eigenvalues with various descriptors based on the ratio between eigenvalues.

    Primarily used for Quantum Theory of Atoms in Molecules (QTAIM).
    """

    def __init__(self, eigenvals, zero_eps=1e-15):
        r"""
        Construct an EigenDescriptor class, only requiring an array of eigenvalues.

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

    def ellipticity(self, index):
        r"""
        Ellipticity of a bond critical point.

        Parameters
        ----------
        index : int
            The index of the eigenvalues.

        Returns
        -------
        float :
            Returns the ratio of the largest two eigenvalues minus one. If second largest eigenvalue
            is zero, then infinity is returned.
        """
        def _largest_two_numbers(arr):
            r"""Get the first two largest values of the array."""
            num1 = -np.inf
            num2 = -np.inf
            for eigen in arr:
                if eigen > num1:
                    num2 = num1
                    num1 = eigen
                elif eigen >= num2:
                    num2 = eigen
            return num1, num2

        eigen1, eigen2 = _largest_two_numbers(self.eigenvals[index])
        if np.abs(eigen2) < self._zero_eps:
            warnings.warn("Second largest eigenvalue is zero.")
            return np.inf

        return (eigen1 / eigen2) - 1.

    def bond_descriptor(self, index):
        r"""
        Bond descriptor is defined as the ratio of average of positive and negative eigenvalues.

        Parameters
        ----------
        index : int
            The index of the eigenvalues.

        Returns
        -------
        float :
            Ratio of average positive eigenvalues to average negative eigenvalues. If there are no
            negative eigenvalues, then None is returned.
        """
        # compute numerator
        pos_eigen = self._eigenvals[index][self._eigenvals[index] > self._zero_eps]
        result = np.sum(pos_eigen)
        if len(pos_eigen) != 0:
            result /= len(pos_eigen)
        # compute denominator
        neg_eigen = self._eigenvals[index][self._eigenvals[index] < -self._zero_eps]
        if len(neg_eigen) != 0:
            result /= (np.sum(neg_eigen) / len(neg_eigen))
        else:
            result = None
        return result

    def eccentricity(self, index):
        r"""
        Eccentricity of the set of eigenvalues.

        This is defined as
        .. math ::
            \sqrt{\frac{\lambda_{max}}{\lambda_{min}}}

        Parameters
        ----------
        index : int
            The index of the eigenvalues.

        Returns
        -------
        float :
            The condition number, the square root of largest eigenval divided by minimum eigenval.
            If one of maximima or mininum is negative, then none is returned.
        """
        max_eigen = np.max(self._eigenvals[index])
        min_eigen = np.min(self._eigenvals[index])
        ratio = max_eigen / min_eigen
        if ratio < 0.:
            return None
        return np.sqrt(max_eigen / min_eigen)

    def index_critical_pt(self, index):
        r"""
        Index of the critical point is the number of negative-curvature directions.

        Parameters
        ----------
        index : int
            The index of the eigenvalues, not to be confused with index of critical point.

        Returns
        -------
        int :
            Number of negative eigenvalues.
        """
        return np.sum(self._eigenvals[index] < -self._zero_eps, axis=0)

    def rank(self, index):
        r"""
        Rank of the critical point.

        The rank of a critical point is the number of positive eigenvalues.
        This is used to classify critical points on it's stability (trajectories going in or out).

        .. math ::
            \sum_{\lambda_i > 0} 1

        Parameters
        ----------
        index : int
            The index of the eigenvalue.

        Returns
        -------
        int :
            The number of non-zero eigenvalues.
        """
        nonzero_eigen = [x for x in self._eigenvals[index] if np.abs(x) > self._zero_eps]
        return len(nonzero_eigen)

    def signature(self, index):
        r"""
        Signature of the critical point.

        This is used to classify saddle critical points (ie certain directions flow outwards and
        certain directions flow inwards).

        .. math ::
            \sum_{\lambda_{i} > 0.}1 - \sum_{\lambda_{i} < 0.}1

        Parameters
        ----------
        index : int
            The index of the eigenvalues.

        Returns
        -------
        int :
            The number of positive eigenvalues minus the number of negative eigenvalues.
        """
        pos_eigen, neg_eigen = self._positive_negative_eigenvalues(index)
        return len(pos_eigen) - len(neg_eigen)

    def morse_critical_pt(self, index):
        r"""
        Return both the rank and signature of the critical point.

        .. math ::
            (\sum_{\lambda_i > 0} 1, \sum_{\lambda_{i} > 0.}1 - \sum_{\lambda_{i} < 0.}1)

        A system is degenerate if it has a zero eigenvalue and consequently, it's critical point
        is said to be "catastrophe". It returns a warning in this case.

        Returns
        -------
        (int, int) :
            Returns the rank and signature of the critical point.
        """
        if np.any(np.abs(self._eigenvals) < self._zero_eps):
            warnings.warn("Near catastrophic eigenvalue (close to zero) been found.")
        return self.rank(index), self.signature(index)

    def _positive_negative_eigenvalues(self, index):
        r"""Return the positive and negative eigenvalues."""
        pos_eigen = []
        neg_eigen = []
        for eigen in self._eigenvals[index]:
            if eigen > self._zero_eps:
                pos_eigen.append(eigen)
            elif eigen < -self._zero_eps:
                neg_eigen.append(eigen)
        return pos_eigen, neg_eigen

    @property
    def eigenvals(self):
        r"""Return two-dimensional set of eigenvalues."""
        return self._eigenvals
