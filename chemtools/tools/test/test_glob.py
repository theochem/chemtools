#!/usr/bin/env python

from chemtools import *


def check_global_properties(global_instance, ip, ea):
    # chemical potential (mu) and hardness (eta)
    mu, eta = -0.5 * (ip + ea), (ip - ea)
    # check the attributes of the gloabl instance
    assert abs(global_instance.ip - ip) < 1.e-6
    assert abs(global_instance.ea - ea) < 1.e-6
    assert abs(global_instance.mu - mu) < 1.e-6
    assert abs(global_instance.eta - eta) < 1.e-6
    assert abs(global_instance.electronegativity + mu) < 1.e-6
    assert abs(global_instance.softness - 1.0 / eta) < 1.e-6
    assert abs(global_instance.n_max + mu / eta) < 1.e-6
    value = 0.5 * mu * mu / eta
    assert abs(global_instance.electrophilicity - value) < 1.e-6
    value = (ip - 3 * ea)**2 / (8 * (ip - ea))
    assert abs(global_instance.nucleofugality - value) < 1.e-6
    value = (3 * ip - ea)**2 / (8 * (ip - ea))
    assert abs(global_instance.electrofugality - value) < 1.e-6


def test_global_H():
    # H atom: IP=13.59843401, EA=0.754195
    g = QuadraticGlobalTool(13.59843401, 0.754195)
    check_global_properties(g, 13.59843401, 0.754195)


def test_global_Mg():
    # Mg atom: IP=7.646235, EA=0.0
    g = QuadraticGlobalTool(7.646235, 0.0)
    check_global_properties(g, 7.646235, 0.0)
