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
"""Test critical point finder."""


from unittest import TestCase

from chemtools.topology.critical import Topology, CriticalPoint

import numpy as np
from numpy.testing import assert_allclose


class TestCriticalPoints(TestCase):
    # Test critical point finder class.
    def gauss_func(self, coors, centers=np.zeros(3), alphas=1):
        # Generate gaussian function value.
        return np.prod(np.exp((-alphas * (coors - centers) ** 2)), axis=-1)

    def gauss_deriv(self, coors, centers=np.zeros(3), alphas=1):
        # Generate 1st order derivative of gaussian function.
        if coors.ndim > 1:
            func_v = self.gauss_func(coors, centers, alphas)[:, None]
        else:
            func_v = self.gauss_func(coors, centers, alphas)
        return -2 * alphas * (coors - centers) * func_v

    def gauss_deriv2(self, coors, centers=np.zeros(3), alphas=1):
        # Generate 2nd order derivative of gaussian function.
        diag = (
            -2 * alphas + 4 * alphas ** 2 * (coors - centers) ** 2
        ) * self.gauss_func(coors, centers, alphas)
        hm = np.diag(diag)
        for i in range(3):
            for j in range(i + 1, 3):
                hm[i][j] = (
                    4
                    * (coors[i] - centers[i])
                    * (coors[j] - centers[j])
                    * self.gauss_func(coors, centers, alphas)
                )
                hm[j][i] = hm[i][j]
        return hm

    def test_instantiation_of_topology_class(self):
        # Test properly initiate topo instance.
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        pts = np.random.rand(4, 3)
        topo = Topology(self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts, coords=coors)
        assert isinstance(topo, Topology)

    def test_center_of_gaussian_is_a_nuclear_critical_point(self):
        r"""Test the center of a single Gaussian is a nuclear critical point."""
        coords = np.array([[0., 0., 0.]])
        pts = np.random.rand(4, 3)
        topo = Topology(self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts, coords)
        topo.find_critical_points()
        # Test that only one critical point is found.
        assert len(topo.cps) == 1
        # Test critical point is center
        assert_allclose(np.array([0., 0., 0.]), topo.cps[0].coordinate)
        # Test it is a nuclear critical point.
        assert len(topo.nna) == 1
        assert_allclose(np.array([0., 0., 0.]), topo.nna[0].coordinate)

    def test_find_critical_pts_of_sum_of_two_gaussians(self):
        r"""
        Test finding critical points of a sum of two Gaussian functions.

        There should be three critical points, two of them are the center of the Gaussians.
        The third should be in-between them.

        """
        # The center of the two Gaussians and the alpha parameters
        atoms = np.array([[-2, -2, -2], [2, 2, 2]])
        alf = 0.3

        def fun_v(coors):
            # Generate value function.
            return self.gauss_func(coors, atoms[0], alphas=alf) + self.gauss_func(
                coors, atoms[1], alphas=alf
            )

        def fun_d(coors):
            # Generate 1st order deriv function.
            return self.gauss_deriv(coors, atoms[0], alphas=alf) + self.gauss_deriv(
                coors, atoms[1], alphas=alf
            )

        def fun_d2(coors):
            # Generate 2nd order deriv function.
            return self.gauss_deriv2(coors, atoms[0], alphas=alf) + self.gauss_deriv2(
                coors, atoms[1], alphas=alf
            )

        pts = np.random.rand(20, 3) * 0.5
        pts = np.vstack(
            (pts, np.array([0.05, 0.05, 0.05]), np.array([-0.05, -0.05, -0.05]))
        )
        tp_ins = Topology(fun_v, fun_d, fun_d2, pts, atoms)
        tp_ins.find_critical_points()
        # Test that three critical points are found
        assert len(tp_ins.cps) == 3
        # Test that two of them are the centers
        assert len(tp_ins.nna) == 2
        assert_allclose(tp_ins.nna[0].coordinate, atoms[0, :], rtol=1e-5)
        assert_allclose(tp_ins.nna[1].coordinate, atoms[1, :], rtol=1e-5)
        assert_allclose(tp_ins.cps[-2].coordinate, atoms[0, :], rtol=1e-5)
        assert_allclose(tp_ins.cps[-1].coordinate, atoms[1, :], rtol=1e-5)
        # Test the bond critical-point inbetween the two centers
        assert len(tp_ins.bcp) == 1
        assert_allclose(tp_ins.bcp[0].coordinate, np.array([0., 0., 0.]), atol=1e-4)
        # Poincare-Hopf is satisfied
        assert tp_ins.poincare_hopf_equation
        # Test all other types of critical points weren't found
        assert len(tp_ins.ccp) == 0
        assert len(tp_ins.rcp) == 0

    def test_find_critical_points_of_three_gaussians_in_triangle_format(self):
        # Test three Gaussians centered at "atoms" with hyper-parameter alpha
        #   The three centers of the Gaussians should form a equaliteral triangle.
        #   Since it is a equaliteral triangle, there should 9 critical points,
        #   First three are the centers, the next three are inbetween the the tree Gaussians
        #   the final critical point is the center of the equilateral triangle.
        atoms = np.array([[-2, -2, 0], [2, -2, 0], [0, 2.0 * (np.sqrt(3.0) - 1.0), 0]])
        alf = 1

        def fun_v(coors):
            return (
                self.gauss_func(coors, atoms[0], alphas=alf)
                + self.gauss_func(coors, atoms[1], alphas=alf)
                + self.gauss_func(coors, atoms[2], alphas=alf)
            )

        def fun_d(coors):
            return (
                self.gauss_deriv(coors, atoms[0], alphas=alf)
                + self.gauss_deriv(coors, atoms[1], alphas=alf)
                + self.gauss_deriv(coors, atoms[2], alphas=alf)
            )

        def fun_d2(coors):
            return (
                self.gauss_deriv2(coors, atoms[0], alphas=alf)
                + self.gauss_deriv2(coors, atoms[1], alphas=alf)
                + self.gauss_deriv2(coors, atoms[2], alphas=alf)
            )
        # Construct a Three-Dimensional Grid
        one_dgrid = np.arange(-2.0, 2.0, 0.1)
        one_dgridz = np.arange(-0.5, 0.5, 0.1)  # Reduce computation time.
        pts = np.vstack(np.meshgrid(one_dgrid, one_dgrid, one_dgridz)).reshape(3,-1).T
        tp_ins = Topology(fun_v, fun_d, fun_d2, pts, atoms)
        tp_ins.find_critical_points()
        # Test that there are 7 critical points (c.p.) in total
        assert len(tp_ins.cps) == 7
        # Test that there are 3 nuclear c.p., 3 bond c.p., and 1 ring cp.p at the center.
        assert len(tp_ins.nna) == 3
        assert len(tp_ins.bcp) == 3
        assert len(tp_ins.rcp) == 1
        assert len(tp_ins.ccp) == 0
        # Test the Nuclear c.p. is the same as the center of Gaussians.
        assert_allclose(tp_ins.nna[0].coordinate, atoms[0], rtol=1e-6)
        assert_allclose(tp_ins.nna[1].coordinate, atoms[1], rtol=1e-6)
        assert_allclose(tp_ins.nna[2].coordinate, atoms[2], rtol=1e-6)
        # Test the bond c.p. is the same as the center between the center of Gaussians.
        assert_allclose(tp_ins.bcp[0].coordinate, (atoms[0] + atoms[1]) / 2.0, atol=1e-3)
        assert_allclose(tp_ins.bcp[1].coordinate, (atoms[1] + atoms[2]) / 2.0, atol=1e-3)
        assert_allclose(tp_ins.bcp[2].coordinate, (atoms[2] + atoms[0]) / 2.0, atol=1e-3)
        # Test the ring c.p. is the as the center of the equilateral triangle.
        #   Calculated using centroid of equilateral triangle online calculator.
        assert_allclose(tp_ins.rcp[0].coordinate, np.array([0 , -0.845, 0]), atol=1e-3)
