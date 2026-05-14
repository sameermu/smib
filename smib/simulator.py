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

import math
from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np

from .models.gencls import GENCLS
from .models.genrou import GENROU
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


# =====================================================================
# GENROU (Phase 2.0) runner with 2x2 saliency-aware network solve
# =====================================================================

def _solve_network_with_genrou(state_arr, genrou: GENROU, network: Network) -> complex:
    """V_terminal at the machine bus for a saliency-aware GENROU.

    GENROU has X'_d != X'_q so a single Norton admittance scalar is only
    approximate.  We instead solve the 2x2 algebraic system directly:

        I_inject_global(V) = Y_bus[0,0] * V + Y_bus[0,1] * V_slack

    where I_inject_global is LINEAR in V (given fixed delta, Eqp, Edp):

        Re(I_inject) = M[0,0]*v_re + M[0,1]*v_im + c[0]
        Im(I_inject) = M[1,0]*v_re + M[1,1]*v_im + c[1]

    The 2x2 matrix M and constants c come from the model's stator
    current equations expressed in (re, im) of the global DQ frame.
    See genrou.py for the underlying algebra.

    Returns the complex V_terminal that satisfies KCL at bus 0.
    """
    genrou.unflatten(state_arr)
    delta = genrou.state["delta"]
    Eqp = genrou.state["Eqp"]
    Edp = genrou.state["Edp"]
    Xdp = genrou.params["Xdp"]
    Xqp = genrou.params["Xqp"]
    c, s = math.cos(delta), math.sin(delta)

    # Linear coefficients of the model's I_inject(V) in re/im of V_global.
    M = np.array([
        [s * c * (1 / Xqp - 1 / Xdp), -c * c / Xqp - s * s / Xdp],
        [s * s / Xqp + c * c / Xdp,    s * c * (1 / Xdp - 1 / Xqp)],
    ])
    c_vec = np.array([
        -Edp * c / Xqp + Eqp * s / Xdp,
        -Edp * s / Xqp - Eqp * c / Xdp,
    ])

    # Network linear coefficients of I_inject(V) at bus 0 (KCL).
    Y = network.ybus()
    Y00 = Y[0, 0]
    Y01 = Y[0, 1]
    N = np.array([
        [Y00.real, -Y00.imag],
        [Y00.imag,  Y00.real],
    ])
    rhs_const = Y01 * network.V_slack
    n_const = np.array([rhs_const.real, rhs_const.imag])

    # Solve (M - N) v = n_const - c_vec.
    A = M - N
    b = n_const - c_vec
    v = np.linalg.solve(A, b)
    return complex(v[0], v[1])


def run_smib_genrou(genrou: GENROU, network: Network,
                    t_end: float = 5.0, h: float = 5e-3,
                    scenarios: Iterable[Callable] = (),
                    init_V: complex | None = None,
                    init_S: complex | None = None,
                    require_steady_state: bool = True) -> SimResult:
    """Run a single-machine SMIB scenario with the 4-state GENROU model.

    Same interface as ``run_smib_gencls`` but uses the saliency-aware
    2x2 algebraic network solve.

    Parameters
    ----------
    genrou : initialised or uninitialised GENROU instance.  If init_V
             and init_S are passed, this function calls
             ``genrou.initialise(init_V, init_S)`` before stepping.
    network : Network instance (already configured with R, X, V_slack).
    t_end : simulation horizon [s].
    h : fixed integration step [s].  2 ms is a safe default for the
        transient time scale; reduce to 0.5 ms if you see corrector
        residual climb.
    scenarios : iterable of ``apply(t_now, dt, network)`` callables.
    init_V, init_S : if provided, initialise the machine before the run.
    require_steady_state : if True (default), enforce that all 4
        derivatives are < 1e-6 at t=0.  Disable for perturbation tests.

    Returns
    -------
    SimResult with time vector, traces, final state, and convergence
    info.
    """
    if init_V is not None and init_S is not None:
        genrou.initialise(init_V, init_S)

    # Verify flat-line floor at t=0.
    V0 = _solve_network_with_genrou(genrou.flatten(), genrou, network)
    genrou.inputs["V_terminal"] = V0
    d0 = genrou.derivatives()
    worst = max(abs(v) for v in d0.values())
    if require_steady_state and worst > 1e-6:
        raise RuntimeError(
            f"GENROU not at steady state at t=0: max|dx/dt| = {worst:.3e}. "
            "Initialisation is inconsistent — check V, S inputs, or pass "
            "require_steady_state=False if you intentionally perturbed."
        )

    n_steps = int(round(t_end / h))
    t = np.linspace(0.0, t_end, n_steps + 1)

    out0 = genrou.algebraic_output()
    traces = {k: np.zeros(n_steps + 1) for k in out0.keys()}
    for k, v in out0.items():
        traces[k][0] = v
    iters_log = np.zeros(n_steps + 1, dtype=int)

    def f_fn(x, V):
        genrou.unflatten(x)
        genrou.inputs["V_terminal"] = V
        d = genrou.derivatives()
        return np.array([d[k] for k in genrou.state_keys])

    def solve_net(x):
        return _solve_network_with_genrou(x, genrou, network)

    x = genrou.flatten()
    for k in range(1, n_steps + 1):
        for sc in scenarios:
            sc(t[k], h, network)

        x, info = trapezoidal_step(x, f_fn, solve_net, h)
        iters_log[k] = info["iters"]

        genrou.unflatten(x)
        genrou.inputs["V_terminal"] = solve_net(x)
        genrou.derivatives()
        out = genrou.algebraic_output()
        for kk, vv in out.items():
            traces[kk][k] = vv

    return SimResult(
        t=t,
        traces=traces,
        final_state=dict(genrou.state),
        info={"mean_iters": float(iters_log.mean()),
              "max_iters": int(iters_log.max())},
    )


# =====================================================================
# GENROU + ST1A (Phase 2.1) runner — multi-model with signal flow
# =====================================================================

def run_smib_genrou_avr(genrou: "GENROU", avr, network: Network,
                        t_end: float = 5.0, h: float = 5e-3,
                        scenarios=(),
                        init_V: complex | None = None,
                        init_S: complex | None = None,
                        require_steady_state: bool = True) -> SimResult:
    """Run a SMIB scenario with GENROU + ST1A coupled together.

    Signal flow each timestep:
        |V_terminal|  →  AVR.inputs.V_terminal_mag
        AVR.algebraic_output.Efd  →  GENROU.params.Efd

    The state vector concatenates GENROU's 4 states with the AVR's 2
    states:

        x = [delta, omega, Eqp, Edp,    Vc, x_LL]
                  GENROU             |    ST1A
    """
    from .models.st1a import ST1A
    from .models.genrou import GENROU
    if not isinstance(genrou, GENROU):
        raise TypeError("genrou must be a GENROU instance")
    if not isinstance(avr, ST1A):
        raise TypeError("avr must be an ST1A instance")

    # Initialise both models in the right order.
    if init_V is not None and init_S is not None:
        genrou.initialise(init_V, init_S)
        avr.initialise(init_V, init_S, Efd_init=genrou.params["Efd"])

    # Compose the combined state (model_order: GENROU first, then ST1A).
    n_g = len(genrou.state_keys)
    n_a = len(avr.state_keys)

    def get_x():
        return np.concatenate([genrou.flatten(), avr.flatten()])

    def set_x(x):
        genrou.unflatten(x[:n_g])
        avr.unflatten(x[n_g:])

    def f_fn(x, V):
        set_x(x)
        # Signal flow: |V| → AVR
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.derivatives()  # to refresh _last_Efd cache
        # AVR Efd → GENROU
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        d_g = genrou.derivatives()
        d_a = avr.derivatives()
        return np.concatenate([
            np.array([d_g[k] for k in genrou.state_keys]),
            np.array([d_a[k] for k in avr.state_keys]),
        ])

    def solve_net(x):
        set_x(x)
        return _solve_network_with_genrou(x[:n_g], genrou, network)

    # Steady-state check.
    V0 = solve_net(get_x())
    avr.inputs["V_terminal_mag"] = abs(V0)
    avr.derivatives()
    genrou.params["Efd"] = avr.algebraic_output()["Efd"]
    genrou.inputs["V_terminal"] = V0
    d0_g = genrou.derivatives()
    d0_a = avr.derivatives()
    worst = max(max(abs(v) for v in d0_g.values()),
                max(abs(v) for v in d0_a.values()))
    if require_steady_state and worst > 1e-5:
        raise RuntimeError(
            f"GENROU+ST1A not at steady state: max|dx/dt| = {worst:.3e}"
        )

    n_steps = int(round(t_end / h))
    t = np.linspace(0.0, t_end, n_steps + 1)

    # Pre-allocate trace arrays — combine outputs from both models.
    #
    # Key-collision policy: when both models report the same key (notably
    # "Efd"), GENROU's value is the *commanded* field voltage that fed
    # into its derivatives this step, and the AVR's value is the *output*
    # of the regulator.  At every steady state and every converged step
    # they agree.  We keep one canonical "Efd" trace driven by the AVR
    # (the live regulator output) and stash GENROU's view under
    # "genrou_Efd" for diagnostic plotting.
    out0_g = genrou.algebraic_output()
    out0_a = avr.algebraic_output()
    traces = {}
    for k, v in out0_g.items():
        new_k = f"genrou_{k}" if k in out0_a else k
        traces[new_k] = np.zeros(n_steps + 1); traces[new_k][0] = v
    for k, v in out0_a.items():
        # AVR keys take precedence on collision.
        traces[k] = np.zeros(n_steps + 1); traces[k][0] = v
    iters_log = np.zeros(n_steps + 1, dtype=int)

    x = get_x()
    for k in range(1, n_steps + 1):
        for sc in scenarios:
            sc(t[k], h, network)

        x, info = trapezoidal_step(x, f_fn, solve_net, h)
        iters_log[k] = info["iters"]

        set_x(x)
        V = solve_net(x)
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()

        out_g = genrou.algebraic_output()
        out_a = avr.algebraic_output()
        for kk, vv in out_g.items():
            target = f"genrou_{kk}" if kk in out_a else kk
            traces[target][k] = vv
        for kk, vv in out_a.items():
            traces[kk][k] = vv

    return SimResult(
        t=t,
        traces=traces,
        final_state={**genrou.state, **avr.state},
        info={"mean_iters": float(iters_log.mean()),
              "max_iters": int(iters_log.max())},
    )


def run_smib_genrou_avr_pss(genrou: "GENROU", avr, pss, network: Network,
                            t_end: float = 5.0, h: float = 5e-3,
                            scenarios=(),
                            init_V: complex | None = None,
                            init_S: complex | None = None,
                            require_steady_state: bool = True) -> SimResult:
    """Run a SMIB scenario with GENROU + ST1A + PSS1A coupled together.

    Signal flow each timestep:

        |V_terminal|        →  AVR.inputs.V_terminal_mag
        ω̄ (rotor slip)      →  PSS.inputs.Delta_omega
        PSS.Vpss            →  AVR.inputs.Vpss
        AVR.Efd             →  GENROU.params.Efd

    The state vector concatenates GENROU's 4 states, the AVR's 2
    states, and the PSS's 3 states:

        x = [delta, omega, Eqp, Edp,    Vc, x_LL,    x_w, x_LL1, x_LL2]
              GENROU                  | ST1A      | PSS1A

    The PSS contributes no current to the network — it's a pure
    measurement-and-feedback block — so the algebraic Y_bus solve is
    structurally identical to Phase 2.1's run_smib_genrou_avr.  The
    PSS adds three more differential states for the trapezoidal
    integrator to track.
    """
    from .models.st1a import ST1A
    from .models.pss1a import PSS1A
    from .models.genrou import GENROU
    if not isinstance(genrou, GENROU):
        raise TypeError("genrou must be a GENROU instance")
    if not isinstance(avr, ST1A):
        raise TypeError("avr must be an ST1A instance")
    if not isinstance(pss, PSS1A):
        raise TypeError("pss must be a PSS1A instance")

    # Init order: GENROU first (knows Efd_init), then AVR (consumes Efd_init
    # and back-solves Vref), then PSS (all zeros at steady state).
    if init_V is not None and init_S is not None:
        genrou.initialise(init_V, init_S)
        avr.initialise(init_V, init_S, Efd_init=genrou.params["Efd"])
        pss.initialise()

    n_g = len(genrou.state_keys)
    n_a = len(avr.state_keys)
    n_p = len(pss.state_keys)

    def get_x():
        return np.concatenate([genrou.flatten(), avr.flatten(), pss.flatten()])

    def set_x(x):
        genrou.unflatten(x[:n_g])
        avr.unflatten(x[n_g:n_g + n_a])
        pss.unflatten(x[n_g + n_a:])

    def f_fn(x, V):
        set_x(x)
        # PSS sees rotor slip ω̄ (stored as genrou.state["omega"]).
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        # Evaluate PSS to update Vpss output.
        pss.derivatives()
        # AVR sees terminal voltage magnitude AND Vpss summed in.
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = pss.algebraic_output()["Vpss"]
        avr.derivatives()  # refresh _last_Efd cache
        # AVR Efd → GENROU's field input.
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        d_g = genrou.derivatives()
        d_a = avr.derivatives()
        d_p = pss.derivatives()
        return np.concatenate([
            np.array([d_g[k] for k in genrou.state_keys]),
            np.array([d_a[k] for k in avr.state_keys]),
            np.array([d_p[k] for k in pss.state_keys]),
        ])

    def solve_net(x):
        set_x(x)
        return _solve_network_with_genrou(x[:n_g], genrou, network)

    # Steady-state DAE-consistency check.
    V0 = solve_net(get_x())
    pss.inputs["Delta_omega"] = genrou.state["omega"]
    pss.derivatives()
    avr.inputs["V_terminal_mag"] = abs(V0)
    avr.inputs["Vpss"] = pss.algebraic_output()["Vpss"]
    avr.derivatives()
    genrou.params["Efd"] = avr.algebraic_output()["Efd"]
    genrou.inputs["V_terminal"] = V0
    d0_g = genrou.derivatives()
    d0_a = avr.derivatives()
    d0_p = pss.derivatives()
    worst = max(max(abs(v) for v in d0_g.values()),
                max(abs(v) for v in d0_a.values()),
                max(abs(v) for v in d0_p.values()))
    if require_steady_state and worst > 1e-5:
        raise RuntimeError(
            f"GENROU+ST1A+PSS1A not at steady state: max|dx/dt| = {worst:.3e}"
        )

    n_steps = int(round(t_end / h))
    t = np.linspace(0.0, t_end, n_steps + 1)

    out0_g = genrou.algebraic_output()
    out0_a = avr.algebraic_output()
    out0_p = pss.algebraic_output()
    traces = {}
    # Build trace dict with collision-aware naming.  Precedence (highest →
    # lowest) is PSS, AVR, GENROU — so a shared key like "Efd" lands as
    # the AVR's live regulator output, and the PSS's outputs (Vpss, y_w,
    # …) keep their own slot.
    for k, v in out0_g.items():
        new_k = k
        if k in out0_a or k in out0_p:
            new_k = f"genrou_{k}"
        traces[new_k] = np.zeros(n_steps + 1); traces[new_k][0] = v
    for k, v in out0_a.items():
        new_k = k
        if k in out0_p:
            new_k = f"avr_{k}"
        traces[new_k] = np.zeros(n_steps + 1); traces[new_k][0] = v
    for k, v in out0_p.items():
        traces[k] = np.zeros(n_steps + 1); traces[k][0] = v
    iters_log = np.zeros(n_steps + 1, dtype=int)

    x = get_x()
    for k in range(1, n_steps + 1):
        for sc in scenarios:
            sc(t[k], h, network)

        x, info = trapezoidal_step(x, f_fn, solve_net, h)
        iters_log[k] = info["iters"]

        set_x(x)
        V = solve_net(x)
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        pss.derivatives()
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = pss.algebraic_output()["Vpss"]
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()

        out_g = genrou.algebraic_output()
        out_a = avr.algebraic_output()
        out_p = pss.algebraic_output()
        for kk, vv in out_g.items():
            target = f"genrou_{kk}" if (kk in out_a or kk in out_p) else kk
            traces[target][k] = vv
        for kk, vv in out_a.items():
            target = f"avr_{kk}" if kk in out_p else kk
            traces[target][k] = vv
        for kk, vv in out_p.items():
            traces[kk][k] = vv

    return SimResult(
        t=t,
        traces=traces,
        final_state={**genrou.state, **avr.state, **pss.state},
        info={"mean_iters": float(iters_log.mean()),
              "max_iters": int(iters_log.max())},
    )
