"""Test critical point finder."""

from unittest import TestCase

from chemtools.topology.critical_pts import Topology, CriticalPoint

import numpy as np
from numpy.testing import assert_allclose


class TestCriticalPoints(TestCase):
    """Test critical point finder class."""

    def gauss_func(self, coors, centers=np.zeros(3), alphas=1):
        # """Generate gaussian function value.
        return np.prod(np.exp((-alphas * (coors - centers) ** 2)), axis=-1)

    def gauss_deriv(self, coors, centers=np.zeros(3), alphas=1):
        # """Generate 1st order derivative of gaussian function."""
        if coors.ndim > 1:
            func_v = self.gauss_func(coors, centers, alphas)[:, None]
        else:
            func_v = self.gauss_func(coors, centers, alphas)
        return -2 * alphas * (coors - centers) * func_v

    def gauss_deriv2(self, coors, centers=np.zeros(3), alphas=1):
        # """Generate 2nd order derivative of gaussian function."""
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

    def test_topo(self):
        # """Test properly initiate topo instance."""
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        pts = np.random.rand(4, 3)
        topo = Topology(coors, self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts)
        assert isinstance(topo, Topology)

    def test_default_cube(self):
        # """Test default cube for points."""
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        topo = Topology(coors, self.gauss_func, self.gauss_deriv, self.gauss_deriv2)
        assert topo._kdtree.data.shape == (60 * 60 * 60, 3)

    def test_construct_cage(self):
        # """Test construct cage among target points."""
        pts = Topology._construct_cage(np.array([0, 0, 0]), 1)
        assert len(pts) == 4
        dis = np.linalg.norm(pts[:] - np.array([0, 0, 0]), axis=-1)
        assert_allclose(dis, np.ones(4) * 4.89898, rtol=1e-5)

    def test_get_gradietn(self):
        # """Test get proper gradient."""
        # initiate topo obj.
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        pts = np.random.rand(4, 3)
        topo = Topology(coors, self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts)
        # get cage points
        pts = Topology._construct_cage(np.array([0.5, 0.5, 0.5]), 0.1)
        g_pts = topo.get_gradient(pts)
        # assert len(g_pts) == 4
        ref_g = self.gauss_deriv(pts)
        assert_allclose(ref_g, g_pts)

    def test_add_critical_points(self):
        # """Test add critical points to class."""
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        pts = np.random.rand(4, 3)
        topo = Topology(coors, self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts)
        pt = CriticalPoint(np.random.rand(3), None, None)
        # bond
        ct_type = -1
        topo._add_critical_point(pt, ct_type)
        assert_allclose(topo._crit_bond[0].point, pt.point, atol=1e-10)
        # maxima
        ct_type = -3
        topo._add_critical_point(pt, ct_type)
        assert_allclose(topo._crit_max[0].point, pt.point, atol=1e-10)
        # ring
        ct_type = 1
        topo._add_critical_point(pt, ct_type)
        assert_allclose(topo._crit_ring[0].point, pt.point, atol=1e-10)
        # cage
        ct_type = 3
        topo._add_critical_point(pt, ct_type)
        assert_allclose(topo._crit_cage[0].point, pt.point, atol=1e-10)

    def test_is_coors_pt(self):
        # """Test same pts as nuclear position."""
        coors = np.array([[1, 1, 1], [-1, -1, -1]])
        pts = np.random.rand(4, 3)
        topo = Topology(coors, self.gauss_func, self.gauss_deriv, self.gauss_deriv2, pts)
        for i in coors:
            result = topo._is_coors_pt(i)
            assert result
        pt = (np.random.rand(4, 3) + 0.1) / 2
        for i in pt:
            result = topo._is_coors_pt(i)
            assert not result

    def test_find_one_critical_pt(self):
        # """Test two atoms critical pts."""
        atoms = np.array([[-2, -2, -2], [2, 2, 2]])
        alf = 0.3

        def fun_v(coors):
            """Generate value function."""
            return self.gauss_func(coors, atoms[0], alphas=alf) + self.gauss_func(
                coors, atoms[1], alphas=alf
            )

        def fun_d(coors):
            """Generate 1st order deriv function."""
            return self.gauss_deriv(coors, atoms[0], alphas=alf) + self.gauss_deriv(
                coors, atoms[1], alphas=alf
            )

        def fun_d2(coors):
            """Generate 2nd order deriv function."""
            return self.gauss_deriv2(coors, atoms[0], alphas=alf) + self.gauss_deriv2(
                coors, atoms[1], alphas=alf
            )

        pts = np.random.rand(20, 3) * 0.5
        pts = np.vstack(
            (pts, np.array([0.05, 0.05, 0.05]), np.array([-0.05, -0.05, -0.05]))
        )
        tp_ins = Topology(atoms, fun_v, fun_d, fun_d2, pts)
        tp_ins.find_critical_pts()
        # one critical pt
        assert len(tp_ins._found_ct) == 1
        assert len(tp_ins._crit_bond) == 1
        # one critical type bond
        assert tp_ins._found_ct_type[0] == -1
        # one critical point at origin
        assert_allclose(tp_ins._found_ct[0], [0, 0, 0], atol=1e-10)

    def test_find_triangle_critical_pt(self):
        # """Test three atom ring critical pts."""
        atoms = np.array([[-2, -2, 0], [2, -2, 0], [0, 1, 0]])
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

        tp_ins = Topology(atoms, fun_v, fun_d, fun_d2, extra=1)
        tp_ins.find_critical_pts()
        assert len(tp_ins._crit_bond) == 3
        assert len(tp_ins._crit_ring) == 1
        assert len(tp_ins._crit_max) == 0
        assert len(tp_ins._crit_cage) == 0
