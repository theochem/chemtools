# -*- coding: utf-8 -*-
# QTAIM is an atoms-in-molecules partitioning package based on
# Richard Bader's Quantum Theory of Atoms in Molecules.
#
# Copyright (C) 2014-2015 The QTAIM Development Team
#
# This file is part of QTAIM.
#
# QTAIM is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# QTAIM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
r"""This file contains functions for finding the critical points."""
import numpy as np
from scipy.optimize import root
from scipy.spatial import KDTree

import warnings

__all__ = ["find_critical_pts", "poincare_hopf_relation"]


class CriticalPoint(object):
    def __init__(self, point, eigenvalues, eigenvectors):
        self.point = point
        self.eigenvalues = eigenvalues
        self.eigenvectors = eigenvectors

    def __repr__(self):
        return "{}".format(self.point)


class Topo(object):
    def __init__(self,
                 molecule,
                 value_func,
                 gradian_func,
                 hess_func,
                 points=None):
        self.molecule = molecule
        self.v_f = value_func
        self.g_f = gradian_func
        self.h_f = hess_func
        # num of the maximum equals to num of atoms
        self._kdtree = KDTree(points)
        self._crit_max = []
        self._crit_bond = []
        self._crit_ring = []
        self._cirt_cage = []

    def add_points(self, points):
        new_points = np.vstack((self._kdtree.data, points))
        self._kdtree = KDTree(new_points)

    def _root_find(self, init_guess):
        sol = root(
            self.g_f,
            x0=init_guess,
            jac=self.h_f,
            tol=1e-15,
            method="hybr",
            options={'maxfev': 1000})
        return sol.success, sol.x

    def get_gradient(self, *points):
        g_list = np.array([self.g_f(i) for i in points])
        return np.linalg.norm(g_list, axis=1)

    @staticmethod
    def _construct_cage(point, length, n_points=4):
        if n_points == 4:
            # 489898 = sqrt(24), 0.942809 = sqrt(8) / 3
            # 0.471405 = sqrt(2) / 3, 0.816497 = sqrt(2 / 3)
            constant = 4.89898 * length
            p1 = point + constant * np.array([0.942809, 0., -0.333333])
            p2 = point + constant * np.array([-0.471405, 0.816497, -0.333333])
            p3 = point + constant * np.array([-0.471405, -0.816497, -0.333333])
            p4 = point + constant * np.array([0., 0., 1.])
        else:
            raise NotADirectoryError(
                "Given args n_point={} is not valid".format(n_points))
        return p1, p2, p3, p4

    def find_critical_pts(self):
        for init_point in self._kdtree.data:
            length, _ = self._kdtree.query(init_point, 4)
            tetrahedral = self._construct_cage(init_point, length[-1])
            g_values = self.get_gradient(tetrahedral)
            central_g = self.g_f(init_point)
            # initial guess points
            if np.linalg.norm(central_g) < min(g_values):
                converge, point = self._root_find(init_point)
                if converge:
                    self._classify_critical_pt(point)

    def _satisfy_poincare_hopf(self):
        pre_hopf = len(self._crit_max) - len(self._crit_bond) + \
                            len(self._crit_ring) - len(self._cirt_cage)
        post_hopf = pre_hopf + len(self.molecule.numbers)
        return post_hopf == 1

    def _classify_critical_pt(self, point, eigen_cutoff=1e-4):
        hess_crit = self.h_f(point)
        eigenvals, eigenvecs = np.linalg.eigh(hess_crit)
        signature = np.sum(np.sign(eigenvals))
        signature_dict = {
            -3: self._crit_max,
            3: self._cirt_cage,
            -1: self._crit_bond,
            1: self._crit_ring,
        }
        # Zero eigenvalues occur with points that are far away.
        # If eigenvalues too small, neglect this critical point
        if np.max(np.abs(eigenvals)) > eigen_cutoff:
            # Add this critical point to it's list
            crit_pt = CriticalPoint(point, eigenvals, eigenvecs)
            signature_dict[signature].append(crit_pt)
        self._satisfy_poincare_hopf()


class TopologyInfo(object):
    __slots__ = [
        "eigenvals", 'eigenvecs', "type", "poincare_hopf", "critical_pts"
    ]

    def __init__(self, numb_maxima):
        self.eigenvals = []
        self.eigenvecs = []
        self.type = []
        self.poincare_hopf = numb_maxima
        self.critical_pts = []

    def update_info(self, eval, evec, classify, phopf, cp):
        self.eigenvals.append(eval)
        self.eigenvecs.append(evec)
        self.type.append(classify)
        self.poincare_hopf += phopf
        self.critical_pts.append(cp)


def poincare_hopf_relation(num_maxima, num_bond, num_ring, num_cage):
    return num_maxima - num_bond + num_ring - num_cage


def _satisfy_poincare_hopf(num_maxima, num_bond, num_ring, num_cage):
    return poincare_hopf_relation(num_maxima, num_bond, num_ring,
                                  num_cage) == 1


def _remove_small_dens_vals(wat_pts, grid_pt_cont, tol):
    return wat_pts[np.where(grid_pt_cont.dens_arr > tol)[0]]


def _is_not_a_maxima_pt(pt, maximas, atom_eps):
    return np.all(np.linalg.norm(pt["x"] - maximas, axis=1) >= atom_eps)


def _density_value_not_small(grid_cont, index, dens_cutoff):
    return grid_cont.dens_arr[index] > dens_cutoff


def _is_not_already_found(critical_pts, root, eps):
    for y in critical_pts:
        if np.linalg.norm(y - root) <= eps:
            return False
    return True


def _passes_critical_pt_conditions(cart_pt, grid_cont, grad_func, dens_cutoff):
    grad_vals = grid_cont.func_vals_on_neighbours(cart_pt, grad_func)
    norm_grad = np.linalg.norm(grad_vals, axis=1)
    if np.all(np.linalg.norm(grad_func(cart_pt)) <= norm_grad):
        index_pt = grid_cont.cartesian_to_index(cart_pt)
        if _density_value_not_small(grid_cont, index_pt, dens_cutoff):
            return True
    return False


def _classify_critical_pt(critical_pt, hessian, eigen_cutoff=1e-4):
    hess_crit = hessian(critical_pt)
    eigenvals, eigenvecs = np.linalg.eigh(hess_crit)
    signature = np.sum(np.sign(eigenvals))

    classification = "Not classified"
    phopf_value = 0.
    # Zero eigenvalues occur with points that are far away.
    if np.max(np.abs(eigenvals)) < eigen_cutoff:
        return eigenvals, eigenvecs, classification, phopf_value
    elif signature == -3:
        classification = "Local Maxima"
        phopf_value = 1.
    elif signature == 3:
        classification = "Cage Critical"
        phopf_value = -1.
    elif signature == -1:
        classification = "Bond Critical"
        phopf_value = -1.
    elif signature == 1:
        classification = "Ring Critical"
        phopf_value = 1.
    return eigenvals, eigenvecs.T, classification, phopf_value


def find_critical_pts(watershed_pts,
                      maxima,
                      grid_cont,
                      grad_func,
                      hessian,
                      atom_eps=1e-3,
                      dens_cutoff=1e-4,
                      eigen_cuttoff=1e-5):
    # TODO: Add vander wall radius
    # TODO: Why do I have the gradient is less than gradient of neighbours.
    options = {"maxfev": 10000, "factor": 0.5, "xtol": 1e-15}
    topo = TopologyInfo(len(maxima))

    for _, wat_pt in enumerate(watershed_pts):
        # obtain near by points value
        if _passes_critical_pt_conditions(wat_pt, grid_cont, grad_func,
                                          dens_cutoff):
            print("Cartesian pt ", wat_pt)
            # solve nonlinear equation with scipy root
            rt = root(
                grad_func,
                wat_pt,
                jac=hessian,
                tol=1e-15,
                method="hybr",
                options=options)

            if rt["success"] and _is_not_a_maxima_pt(rt, maxima, atom_eps):
                if _is_not_already_found(topo.critical_pts, rt["x"], atom_eps):
                    topo_cp = _classify_critical_pt(rt["x"], hessian,
                                                    eigen_cuttoff)
                    if topo_cp[2] != "Not classified":
                        topo.update_info(topo_cp[0], topo_cp[1], topo_cp[2],
                                         topo_cp[3], rt["x"])

    if topo.poincare_hopf != 1:
        warnings.warn("Poincare Hopf value is not one! You may need to change"
                      "the different settings on the critical point finder or"
                      "try again with a more refine grid.")
    return {
        "critical_pts": np.array(topo.critical_pts),
        "eigenvals": topo.eigenvals,
        "eigenvecs": topo.eigenvecs,
        "type": topo.type,
        "poincare_hopf": topo.poincare_hopf
    }


def find_critical_pts_finite():
    pass
