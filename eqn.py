import numpy as np
import sympy as sp
from scipy.optimize import root


def solve(expr, n_energies, n_symbol, n_eval, n_eval_symbol=None, guess=None, ndiff=1, opts=None):
    """
    Doc string.

    Parameters
    ----------
    expr : sp.Expr
        The expression to differentiate.
    n_energies : dict
        int (electron-number) keys, float (energy) values.
    n_symbol : sp.Symbol
        The symbol in `expr` that represents electron-number.
    n_eval : int
        The electron number at which to evaluate `expr`.
    n_eval_symbol: sp.Symbol, optional
        The symbol in `expr` that represents electron-number at which to evaluate `expr`.
        If not specified, assume that it is already expressed numerically in `expr`.
    guess : dict, optional
        Initial guess at the values of `expr`'s parameters.  sp.Symbol keys, float
        values.
    ndiff : int, optional
        The order of the derivative to take.  Defaults to 1.
    opts : dict, optional
        Additional options to pass to the internal SciPy solver.

    Returns
    -------
    dexpr : sp.Expr
        The expression equivalent to d(`expr`)/d(`n_symbol`).

    """

    # Handle `opts` argument
    if not opts:
        opts = {}

    # Parse the non-electron-number parameters of `expr`
    if n_eval_symbol:
        expr = expr.subs(n_eval_symbol, n_eval)
    params = [symbol for symbol in expr.atoms() if isinstance(symbol, sp.Symbol) and
              (symbol is not n_symbol) and (symbol is not n_eval_symbol)]

    # Fill in missing non- energy/electron-number symbols from `guess`
    if guess:
        guess.update({symbol: 0. for symbol in params if symbol not in guess})

    # Create initial guess at `params`' values if `guess` is not given
    else:
        guess = {symbol: 0. for symbol in params}

    # Construct system of equations to _solve for `params`
    assert len(params) <= len(n_energies), \
        "The system is underdetermined.  Inlcude more (electron-number, energy) pairs."
    system_eqns = []
    for nelec, energy in n_energies.iteritems():
        eqn = sp.lambdify((params,), expr.subs(n_symbol, nelec) - energy, "numpy")
        system_eqns.append(eqn)

    # Solve the system of equations for `params`
    def objective(args):
        return np.array([fun(args) for fun in system_eqns])

    root_guess = np.array([guess[symbol] for symbol in params])
    roots = root(objective, root_guess, **opts)
    assert roots.success, \
        "The system of equations could not be solved."

    # Substitute in the values of `params`
    expr = expr.subs([(params[i], roots.x[i]) for i in range(len(params))])

    # Return the derivative of `expr` wrt `n_symbol`
    for i in range(ndiff):
        expr = expr.diff(n_symbol)
    return expr

# vim: set textwidth=90 :
