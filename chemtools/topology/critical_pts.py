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
from scipy.linalg import eig

import warnings


__all__ = ["find_critical_pts", "poincare_hopf_relation"]


class _TopologyInfo(object):
    __slots__ = ["eigenvals", 'eigenvecs', "type", "poincare_hopf",
                 "critical_pts"]

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
    return poincare_hopf_relation(num_maxima, num_bond, num_ring, num_cage) == 1


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
    eigenvals, eigenvecs = eig(hess_crit)
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


def find_critical_pts(watershed_pts, maxima, grid_cont, grad_func, hessian,
                      atom_eps=1e-3, dens_cutoff=1e-4, eigen_cuttoff=1e-5):
    # TODO: Add vander wall radius
    # TODO: Why do I have the gradient is less than gradient of neighbours.
    options = {"maxfev": 10000, "factor": 0.5, "xtol": 1e-15}
    topo = _TopologyInfo(len(maxima))

    for i, wat_pt in enumerate(watershed_pts):
        if _passes_critical_pt_conditions(wat_pt, grid_cont, grad_func, dens_cutoff):
            print("Cartesian pt ", wat_pt)
            rt = root(grad_func, wat_pt, jac=hessian, tol=1e-15, method="hybr",
                      options=options)

            if rt["success"] and _is_not_a_maxima_pt(rt, maxima, atom_eps):
                if _is_not_already_found(topo.critical_pts, rt["x"], atom_eps):
                    topo_cp = _classify_critical_pt(rt["x"], hessian, eigen_cuttoff)
                    if topo_cp[2] != "Not classified":
                        topo.update_info(topo_cp[0], topo_cp[1], topo_cp[2],
                                         topo_cp[3], rt["x"])

    if topo.poincare_hopf != 1:
        warnings.warn("Poincare Hopf value is not one! You may need to change"
                      "the different settings on the critical point finder or"
                      "try again with a more refine grid.")
    return {"critical_pts": np.array(topo.critical_pts),
            "eigenvals": topo.eigenvals, "eigenvecs": topo.eigenvecs,
            "type": topo.type, "poincare_hopf": topo.poincare_hopf}


def find_critical_pts_finite():
    pass
