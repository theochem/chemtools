#!/usr/bin/env python

from chemtools import *
import numpy.testing

def check_global_properties(global_instance, energy_zero, energy_plus, energy_minus):
    ip = energy_minus - energy_zero
    ea = energy_zero - energy_plus
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
    g = QuadraticGlobalTool(0, -0.754195, 13.59843401)
    check_global_properties(g, 0, -0.754195, 13.59843401)


def test_global_Mg():
    # Mg atom: IP=7.646235, EA=0.0
    g = QuadraticGlobalTool(0, 0, -7.646235)
    check_global_properties(g,0, 0, -7.646235)

def check_global_linear_properties(global_instance, energy_zero, energy_plus, energy_minus):
    # mu_plus, mu_minus, mu_zero
    ip = energy_minus - energy_zero
    ea = energy_zero - energy_plus
    # chemical potential (mu) and hardness (eta)
    mu_plus, mu_minus, mu_zero = -ea, -ip, -0.5*(ea + ip)
    # check the attributes of the global instance
    numpy.testing.assert_almost_equal(global_instance.mu_plus, mu_plus, 6)
    numpy.testing.assert_almost_equal(global_instance.mu_minus, mu_minus, 6)
    numpy.testing.assert_almost_equal(global_instance.mu_zero, mu_zero, 6)

def test_global_linear_H():
    # H atom: IP=13.59843401, EA=0.754195
    g = LinearGlobalTool(0, -0.754195, 13.59843401)
    check_global_linear_properties(g, 0, -0.754195, 13.59843401)

def test_global_linear_Mg():
    # Mg atom: IP=7,646235, EA=0.0
    g = LinearGlobalTool(0, 0, -7.646235)
    check_global_linear_properties(g, 0, 0, -7.646235)
    g = QuadraticGlobalTool(0, 0, -7.646235)
    check_global_properties(g, 0, 0, -7.646235)
