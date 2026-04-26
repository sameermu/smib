"""Implicit trapezoidal integrator with fixed-point corrector.

This is the solver PSSE uses (with minor variants). A-stable, second-order
accurate. Step size h is fixed; pick it based on fastest time constant.

For a general DAE:
    dx/dt = f(x, V)           (differential part: all model states)
    0     = g(x, V)           (algebraic part:   network equations)

the partitioned scheme we use here does:
  1) at time t_k: (x_k, V_k) consistent
  2) predict:  x_pred = x_k + h * f(x_k, V_k)
  3) solve network at x_pred to get V_pred
  4) corrector fixed-point:
        x_{k+1} <- x_k + 0.5*h*( f(x_k, V_k) + f(x_{k+1}, V_{k+1}) )
     with V_{k+1} re-solved each time
  5) when |x_{k+1} - x_prev| < tol, accept

Simultaneous Newton (xdot + network jointly) is more robust for stiff DAEs
but we don't need it for SMIB. If you see convergence failures later, swap.
"""
import numpy as np


def trapezoidal_step(x0, f_fn, solve_network_fn, h,
                     tol=1e-8, max_iter=15):
    """One trapezoidal step.

    Parameters
    ----------
    x0 : np.ndarray
        State vector at time t_k.
    f_fn : callable
        f(x, V) -> dx/dt. V comes from solve_network_fn.
    solve_network_fn : callable
        V = solve_network_fn(x). Returns current bus voltages (complex array).
    h : float
        Step size in seconds.
    tol : float
        Absolute tolerance on corrector fixed-point.
    max_iter : int
        Max corrector iterations.

    Returns
    -------
    x1 : np.ndarray
    info : dict
        {"iters": int, "residual": float}
    """
    V0 = solve_network_fn(x0)
    f0 = f_fn(x0, V0)
    # predictor (explicit Euler)
    x = x0 + h * f0
    for k in range(max_iter):
        V = solve_network_fn(x)
        f = f_fn(x, V)
        x_new = x0 + 0.5 * h * (f0 + f)
        res = float(np.max(np.abs(x_new - x)))
        x = x_new
        if res < tol:
            return x, {"iters": k + 1, "residual": res}
    raise RuntimeError(
        f"Trapezoidal corrector did not converge: residual={res:.3e}"
    )
