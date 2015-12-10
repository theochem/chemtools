import sympy as sp
import numpy as np
from eqn import solve


n, n0, a, b, g = sp.symbols("n, n0, a, b, g")
expr = a * sp.exp(-g * (n - n0)) + b

n_eval = 5
n_energies = {4: 6., 5: 5., 6: 3.}
guess = {a: -1., b: 4., g: -np.log(3.)}

dexpr_computed = solve(expr, n_energies, n, n_eval, n0, guess)

answer = {a: -2., b: 7., g: -np.log(2.)}
dexpr_actual = expr.subs([ (param, value) for (param, value) in answer.iteritems() ])
dexpr_actual = dexpr_actual.subs(n0, n_eval)
dexpr_actual = dexpr_actual.diff(n)

vals_computed = [ dexpr_computed.subs(n, value) for value in range(0,10) ]
vals_actual = [ dexpr_actual.subs(n, value) for value in range(0,10) ]
for i in range(0,10):
    assert np.abs(vals_computed[i] - vals_actual[i]) < 1e-9

# vim: set textwidth=90 :
