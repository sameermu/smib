"""Time-loop runner: GENCLS + Network + scenarios -> trace dict.

Designed for SMIB with a single dynamic model (Phase 1 — GENCLS).  When
we add IBR/AVR/PSS/Gov in later phases this module will grow to
orchestrate multiple models, but the pattern is the same:
  - solve the augmented network for V given current model states
  - call model.derivatives() to get dx/dt
  - hand (x, dx/dt) to the trapezoidal integrator
  - log every algebraic_output() each step

The augmented Ybus substitution
-------------------------------
A GENCLS injects current  I_inj = (E' - V) / (j*X'd), which depends on
V.  Substituting into the network equation
  Y00*V0 + Y01*Vslack = I_inj
yields the explicit form
  V0 = (E' / (j*X'd) - Y01*Vslack) / (Y00 + 1/(j*X'd))
which is exactly "include the machine's Norton admittance in Y00".  No
fixed-point iteration on the algebraic loop is needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np

from .models.gencls import GENCLS
from .network import Network
from .solver import trapezoidal_step


def _solve_network_with_gencls(state_arr, gencls: GENCLS, network: Network) -> complex:
    """V_terminal at the machine bus, with the Norton admittance folded
    into Y00 so the algebraic loop closes in one shot."""
    # Sync state into the model so derivatives use the right delta.
    gencls.unflatten(state_arr)
    Y_mach = gencls.norton_admittance()
    I_src = gencls.norton_source()
    Y = network.ybus()
    V0 = (I_src - Y[0, 1] * network.V_slack) / (Y[0, 0] + Y_mach)
    return V0


@dataclass
class SimResult:
    """Container for everything a notebook wants from a run."""
    t: np.ndarray
    traces: dict
    final_state: dict
    info: dict


def run_smib_gencls(gencls: GENCLS, network: Network,
                    t_end: float = 5.0, h: float = 5e-3,
                    scenarios: Iterable[Callable] = (),
                    init_V: complex | None = None,
                    init_S: complex | None = None,
                    require_steady_state: bool = True) -> SimResult:
    """Run a single-machine SMIB scenario.

    Parameters
    ----------
    gencls    : initialised or uninitialised GENCLS instance.  If
                init_V and init_S are passed, this function calls
                gencls.initialise(init_V, init_S) before stepping.
    network   : Network instance (already configured with R, X, V_slack).
    t_end     : simulation horizon in seconds.
    h         : fixed integration step.  5 ms is a safe default for
                GENCLS swing dynamics; reduce to 1 ms if you see the
                trapezoidal residual climb.
    scenarios : iterable of `apply(t_now, dt, network)` callables that
                may mutate the network mid-run (faults, V steps, ...).
    init_V, init_S : if provided, initialise the machine before the run.

    Returns
    -------
    SimResult with time vector, all logged traces, and convergence info.
    """
    if init_V is not None and init_S is not None:
        gencls.initialise(init_V, init_S)

    # Verify flat-line floor at t=0 (suppressible for perturbation tests).
    V0 = _solve_network_with_gencls(gencls.flatten(), gencls, network)
    gencls.inputs["V_terminal"] = V0
    d0 = gencls.derivatives()
    worst = max(abs(v) for v in d0.values())
    if require_steady_state and worst > 1e-6:
        raise RuntimeError(
            f"GENCLS not at steady state at t=0: max|dx/dt| = {worst:.3e}. "
            "Initialisation is inconsistent — check V, S inputs, or pass "
            "require_steady_state=False if you intentionally perturbed the state."
        )

    # Time grid.
    n_steps = int(round(t_end / h))
    t = np.linspace(0.0, t_end, n_steps + 1)

    # Pre-allocate trace arrays from the algebraic output dict.
    out0 = gencls.algebraic_output()
    traces = {k: np.zeros(n_steps + 1) for k in out0.keys()}
    for k, v in out0.items():
        traces[k][0] = v

    iters_log = np.zeros(n_steps + 1, dtype=int)

    # Closures for the integrator.
    def f_fn(x, V):
        gencls.unflatten(x)
        gencls.inputs["V_terminal"] = V
        d = gencls.derivatives()
        return np.array([d[k] for k in gencls.state_keys])

    def solve_net(x):
        return _solve_network_with_gencls(x, gencls, network)

    x = gencls.flatten()
    for k in range(1, n_steps + 1):
        # Apply any scheduled scenario events.
        for sc in scenarios:
            sc(t[k], h, network)

        x, info = trapezoidal_step(x, f_fn, solve_net, h)
        iters_log[k] = info["iters"]

        # Sync state and log algebraic outputs.
        gencls.unflatten(x)
        gencls.inputs["V_terminal"] = solve_net(x)
        gencls.derivatives()  # refresh _last_Pe, _last_I
        out = gencls.algebraic_output()
        for kk, vv in out.items():
            traces[kk][k] = vv

    return SimResult(
        t=t,
        traces=traces,
        final_state=dict(gencls.state),
        info={"mean_iters": float(iters_log.mean()), "max_iters": int(iters_log.max())},
    )
