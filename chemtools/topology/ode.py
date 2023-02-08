import numpy as np
from scipy.spatial.distance import cdist

from scipy.integrate import solve_ivp

__all__ = ["steepest_ascent_rk45", "gradient_path"]


def _get_normalized_gradient_func(grad_func):
    r"""Returns the normalized version of the function."""
    def norm_grad_func(x):
        grad = grad_func(x)
        return grad / np.linalg.norm(grad, axis=1)[:, None]
    return norm_grad_func


def _RK45_step(pts, grad_func, step_size):
    r"""
    Runge-Kutta fourth and five-order step for the following ode system:

    .. math::
        \frac{d(r(s))}{ds} = \nabla \rho (r(s)),

    where :math:`\nabla \rho(x)` is the gradient of a function.

    Parameters
    ----------
    pts : ndarray(M, 3)
        Points to take a step in.
    grad_func: callable(ndarray(M,3), ndarray(M,3))
        Callable function that takes in points and obtains the gradient.
    step_size: float
        Stepsize for the step

    Returns
    -------
    (ndarray(M,3), ndarray(M,3)):
        Returns fourth-order and fifth-order Runge-Kutta step, respectively.

    """
    # Fourth and Fifth-Order Runge-Kutta
    k1 = step_size * grad_func(pts)
    k2 = step_size * grad_func(pts + 0.4 * k1)
    k3 = step_size * grad_func(pts + (3.0 / 32) * k1 + (9.0 / 32.0) * k2)
    k4 = step_size * grad_func(pts + (1932 / 2197) * k1 - (7200 / 2197) * k2 + (7296 / 2197) * k3)
    k5 = step_size * grad_func(
        pts + (439 / 216) * k1 - 8.0 * k2 + (3680 / 513) * k3 - (845 / 4104) * k4
    )
    k6 = step_size * grad_func(
        pts
        - (8.0 / 27.0) * k1
        + 2.0 * k2
        - (3544 / 2565) * k3
        + (1859 / 4104) * k4
        - (11 / 40) * k5
    )

    # Get the fourth and five-order approximation
    y_four = pts + (25.0 / 216.0) * k1 + (1408 / 2565) * k3 + (2197 / 4101) * k4 - k5 / 5.0
    y_five = (
        pts +
        (16.0 / 135.0) * k1
        + (6656 / 12825) * k3
        + (28561 / 56430) * k4
        - (9.0 / 50.0) * k5
        + (2.0 / 55.0) * k6
    )
    return y_four, y_five


def steepest_ascent_rk45(
    initial_pts, dens_func, grad_func, beta_spheres, maximas, ss_0=1e-7,
    tol=1e-7, max_ss=0.25, terminate_if_other_basin_found=False
):
    r"""
    Solves the following problem ODE using Runge-Kutta of order 4(5) with adaptive step-size

    .. math::
        \frac{d(r(t))}{dt} = \frac{\nabla \rho(r(t))}{|| \rho(r() ||}

    over a set of points.

    Parameters
    ----------
    initial_pts: ndarray(N, 3)
        Initial points to solve for steepest-ascent/backtracing.
    dens_func: callable(ndarray(N,3), ndarray(N))
        The electron density function.
    grad_func: callable(ndarray(N,3), ndarray(N,3))
        The gradient of the electron density.
    beta_spheres: ndarray(M,)
        The beta-sphere/trust-region radius of each atom. These are spheres
        centered at each maxima that reduces the convergence of each point.
    maximas: ndarray(M,3)
        The position of each atoms in the molecule.
    ss_0: float, optional
        The initial step-size of the ODE (RK45) solver.
    max_ss: float, optional
        Maximum step-size of the ODE (RK45) solver.
    tol: float, optional
        Tolerance for the adaptive step-size.
    terminate_if_other_basin_found : bool
        If true, then if multiple basin values were found, then the ODE solver will exit.
        If false, then the ODE solver will run until all points enter one of the
        beta-sphere/trust-region.

    Returns
    -------
    ndarray(N,):
        Integer array that assigns each point to a basin/maxima/atom of the molecule.
        If value is negative one, then the point wasn't assigned to a basin.
    
    """
    norm_grad_func = _get_normalized_gradient_func(grad_func)

    numb_pts = initial_pts.shape[0]
    if isinstance(ss_0, float):
        ss = np.ones((numb_pts, 1)) * ss_0
    elif isinstance(ss_0, np.ndarray):
        if not ss_0.shape[1] == 1:
            raise ValueError(f"Steps-size {ss_0.shape} should have shape of the form (N, 1).")
        ss = ss_0.copy()
    else:
        raise TypeError(f"Step-size ss_0 {type(ss_0)} should have type float or array.")

    pts = initial_pts.copy()
    dens_vals0 = dens_func(initial_pts)
    # print("Intial density values ", dens_vals0)

    assigned_basins = (-1) * np.ones((numb_pts,), dtype=np.int)
    not_found_indices = np.arange(numb_pts)
    first_basin = -1  # First basin value that was found
    print("START STEEPEST-ASCENT")
    import time
    while len(not_found_indices) != 0:
        #start = time.time()
        y_four, y_five = _RK45_step(pts, norm_grad_func, ss)
        # final = time.time()
        # print("RK Step ", final - start)

        # Update step-size
        # print("Step size used", ss)
        ss = (tol * ss / (2.0 * np.linalg.norm(y_five - y_four, axis=1)[:, None]))**0.25
        ss[ss > max_ss] = max_ss

        # Get density values and if the density-values decreased, reduce step-size
        dens_vals1 = dens_func(y_five)
        indices = np.where(dens_vals1 <= dens_vals0)[0]
        # print("New density values ", dens_vals1)
        # print("Indices that density decreaed ", indices)#, dens_vals1[indices], dens_vals0[indices])
        if len(indices) != 0:
            # print("Density Decreased")
            # print("Gradients here", grad_func(pts[indices, :]))
            y_five[indices, :] = pts[indices, :]
            ss[indices] *= 0.25
        # TODO: Check here if the density is equal to the isosurface value, then stop.

        # Check any points are within the beta-spheres and remove them if they converged.
        dist_maxima = cdist(y_five, maximas)
        # print("Distance maxima", dist_maxima)
        beta_sph = dist_maxima <= beta_spheres
        which_basins = np.where(beta_sph)
        # print("beta_sphereS ", beta_spheres)
        # print("which pts are within basin based on beta-sphere", which_basins)
        # print("Conv to bet asphere", which_converged)
        if len(which_basins[0]) != 0:
            assigned_basins[not_found_indices[which_basins[0]]] = which_basins[1]
            not_found_indices = np.delete(not_found_indices, which_basins[0])
            y_five = np.delete(y_five, which_basins[0], axis=0)
            ss = np.delete(ss, which_basins[0])[:, None]
            dens_vals1 = np.delete(dens_vals1, which_basins[0])

            # Terminate early if multiple basins were found
            if terminate_if_other_basin_found:
                unique_basins = np.unique(which_basins[1])
                # If the number of basins is greater than one then exit
                if len(unique_basins) > 1:
                    return assigned_basins
                else:
                    # If a new basin is found then exit
                    if first_basin is None:
                        first_basin = unique_basins[0]
                    elif first_basin != unique_basins[0]:
                        return assigned_basins

        # print("New indices to still do: ", not_found_indices)

        # Update next iteration
        pts = y_five.copy()
        # print(pts)
        dens_vals0 = dens_vals1

        # input("Next step")
    # print("Final basins ", assigned_basins)
    return assigned_basins


def gradient_path(pt, grad_func, use_norm=False, maximas=None, t_span=(0, 1000), method="LSODA", max_step=100,
                  t_inc=400, max_tries=10, first_step=1e-3, beta_spheres=-np.inf):
    r"""
    Solves the final point of steepest-ascent starting at a single pt.

    If maximas is not provided, then termination occurs due to convergence.

    Parameters
    ----------
    pt : ndarray(3,)
        The initial point of the ODE.
    grad_func : callable(ndarray(N, 3) -> ndarray(M, 3))
        The gradient function for the steepest-ascent
    use_norm: bool
        If true, then the grad_func output is normalized. This isn't useful for optimizing
        the centers to obtain the local maximas.
    maximas: (ndarray(M, 3), None)
        The maximas of the density. If it isn't provided then termiantion occurs due to
        convergence between two consequent points of the ODE.
    t_span: (float, float)
        The lower-bound and upper-bound of solving the ODE time-step.

    Returns
    -------
    ndarray(3,)
        The final point of the steepest-ascent path.

    """
    is_converged = False
    y0 = pt.copy()
    numb_times = 0

    norm_grad_func = grad_func
    if use_norm:
        norm_grad_func = _get_normalized_gradient_func(grad_func)

    def grad(t, x):
        return norm_grad_func(np.array([x]))[0]

    while not is_converged and numb_times < max_tries:
        sol = solve_ivp(
            grad,
            y0=y0,
            t_span=t_span,
            method=method,
            max_step=max_step,
            first_step=first_step,
        )
        # print(sol)
        assert sol["success"], "ODE was not successful."
        # If it is close to a maxima or within any of the beta-spheres, then stop.
        if maximas is not None:
            last_y_val = sol["y"][:, -1]
            dist_maxima = np.linalg.norm(last_y_val - maximas, axis=1)
            if np.any(dist_maxima < 0.1) or np.any(dist_maxima <= beta_spheres):
                return sol["y"][:, -1]
        # if maximas not specified, then just look at if it converged.
        else:
            convergence = np.linalg.norm(sol["y"][:, -2] - sol["y"][:, -1])
            if convergence < 1e-1:
                return sol["y"][:, -1]
        # No convergence occured, so increaes t-span.
        # print(sol["y"][:, -1], t_span, "YE")
        t_span = (t_span[1], t_span[1] + t_inc)
        y0 = sol["y"][:, -1]
        numb_times += 1

    if numb_times == max_tries:
        raise RuntimeError(f"No convergence in normalized_gradient path pt {pt},"
                           f" solution {sol['y'][:, -1]}, t_span {t_span}")
