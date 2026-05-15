"""Frequency-domain analysis tools for smib.

Three layers:

1. **Numerical linearisation** — finite-difference Jacobians of the
   nonlinear DAE around any operating point.  Works for any combination
   of models in the smib stack (GENROU alone, GENROU + ST1A, GENROU +
   ST1A + PSS1A).  Returns continuous-time state-space matrices
   ``(A, B, C, D)`` ready for scipy.signal.

2. **Transfer-function utilities** — build SISO transfer functions
   from MIMO state-space, multiply TFs (series composition), compose
   open-loop ``L(s) = C(s) · P(s)`` for AVR / PSS design, evaluate
   Bode magnitude/phase and Nyquist locus.

3. **Stability margins + eigenvalue analysis** — compute phase margin,
   gain margin, and crossover frequencies from a Bode sweep.
   ``eigenmodes(A)`` returns the closed-loop poles with their damping
   ratio, natural frequency, and dominant participating state.

All ground-truth checks use scipy.signal for the Laplace inverse
work; this module is mostly bookkeeping + smib-specific glue.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np
import scipy.signal as sig


# =====================================================================
# Numerical linearisation
# =====================================================================

def numerical_jacobian(f: Callable[[np.ndarray], np.ndarray],
                       x0: np.ndarray,
                       eps: float = 1e-6) -> np.ndarray:
    """Central-difference Jacobian ∂f/∂x at x0.

    For a function ``f : R^n -> R^m``, returns the (m × n) Jacobian
    matrix ``J[i,j] = ∂f_i/∂x_j``.
    """
    x0 = np.asarray(x0, dtype=float)
    f0 = np.asarray(f(x0), dtype=float)
    n = x0.size
    m = f0.size
    J = np.zeros((m, n))
    for j in range(n):
        e = np.zeros(n); e[j] = eps
        f_plus = np.asarray(f(x0 + e), dtype=float)
        f_minus = np.asarray(f(x0 - e), dtype=float)
        J[:, j] = (f_plus - f_minus) / (2 * eps)
    return J


# =====================================================================
# Closed-loop linearisation: GENROU alone
# =====================================================================

@dataclass
class StateSpace:
    """Continuous-time state-space (A, B, C, D) with named axes."""
    A: np.ndarray
    B: np.ndarray
    C: np.ndarray
    D: np.ndarray
    state_names: tuple
    input_names: tuple
    output_names: tuple


def linearise_genrou(genrou, network, V_op: complex, S_op: complex,
                     eps: float = 1e-6) -> StateSpace:
    """Linearise bare GENROU + network around (V_op, S_op).

    State x = (delta, omega, Eqp, Edp).
    Inputs u = (Efd, Pm).  Both are *parameters* of GENROU in smib,
    so we treat them as small-signal exogenous inputs here.
    Outputs y = (|V|, P, Q, Id, Iq, Te, delta, omega, Eqp).

    Returns a StateSpace with A (4×4), B (4×2), C (9×4), D (9×2).
    Uses scipy-compatible sign convention: dx/dt = A·x + B·u,
    y = C·x + D·u, both linearised around the operating point so
    `x` and `u` represent *deviations* from the OP.
    """
    from .simulator import _solve_network_with_genrou
    # Re-init at the operating point to make sure state is consistent.
    genrou.initialise(V_op, S_op)
    x0 = genrou.flatten()
    Efd0 = genrou.params["Efd"]
    Pm0 = genrou.params["Pm"]

    def _set_and_eval_derivs(x: np.ndarray) -> np.ndarray:
        genrou.unflatten(x)
        V = _solve_network_with_genrou(x, genrou, network)
        genrou.inputs["V_terminal"] = V
        d = genrou.derivatives()
        return np.array([d[k] for k in genrou.state_keys])

    def _set_and_eval_outputs(x: np.ndarray) -> np.ndarray:
        genrou.unflatten(x)
        V = _solve_network_with_genrou(x, genrou, network)
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()
        out = genrou.algebraic_output()
        return np.array([out["|V|"], out["P"], out["Q"],
                         out["Id"], out["Iq"], out["Pe"],
                         genrou.state["delta"], genrou.state["omega"],
                         genrou.state["Eqp"]])

    def _set_and_eval_derivs_u(u: np.ndarray) -> np.ndarray:
        genrou.params["Efd"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        return _set_and_eval_derivs(x0)

    def _set_and_eval_outputs_u(u: np.ndarray) -> np.ndarray:
        genrou.params["Efd"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        return _set_and_eval_outputs(x0)

    # Save and restore Efd/Pm so jacobians don't permanently shift them.
    A = numerical_jacobian(_set_and_eval_derivs, x0, eps)
    C = numerical_jacobian(_set_and_eval_outputs, x0, eps)
    # Reset to nominal so the u-Jacobians use the OP starting point.
    genrou.params["Efd"] = Efd0
    genrou.params["Pm"] = Pm0
    u0 = np.array([Efd0, Pm0])
    B = numerical_jacobian(_set_and_eval_derivs_u, u0, eps)
    D = numerical_jacobian(_set_and_eval_outputs_u, u0, eps)
    # Restore once more after the u-Jacobian's last call.
    genrou.params["Efd"] = Efd0
    genrou.params["Pm"] = Pm0
    genrou.unflatten(x0)

    return StateSpace(
        A=A, B=B, C=C, D=D,
        state_names=("delta", "omega", "Eqp", "Edp"),
        input_names=("Efd", "Pm"),
        output_names=("|V|", "P", "Q", "Id", "Iq", "Te",
                      "delta", "omega", "Eqp"),
    )


# =====================================================================
# Closed-loop linearisation: GENROU + ST1A
# =====================================================================

def linearise_genrou_avr(genrou, avr, network, V_op: complex, S_op: complex,
                         eps: float = 1e-6) -> StateSpace:
    """Linearise GENROU + ST1A closed-loop around (V_op, S_op).

    State x = (delta, omega, Eqp, Edp, Vc, x_LL).  6 states.
    Inputs u = (Vref, Pm, Vpss).  3 inputs.
    Outputs as in :func:`linearise_genrou`, plus Efd.

    Vpss is treated as an external input so the analysis can be reused
    for PSS design (you compute the transfer function from Vpss to Te
    and design the PSS to phase-compensate it).
    """
    from .simulator import _solve_network_with_genrou
    genrou.initialise(V_op, S_op)
    avr.initialise(V_op, S_op, Efd_init=genrou.params["Efd"])

    n_g = len(genrou.state_keys)
    n_a = len(avr.state_keys)

    def _combine_state():
        return np.concatenate([genrou.flatten(), avr.flatten()])

    def _set_combined(x):
        genrou.unflatten(x[:n_g])
        avr.unflatten(x[n_g:])

    x0 = _combine_state()
    Vref0 = avr.params["Vref"]
    Pm0 = genrou.params["Pm"]
    Vpss0 = avr.inputs.get("Vpss", 0.0)

    def _eval_derivs(x):
        _set_combined(x)
        V = _solve_network_with_genrou(x[:n_g], genrou, network)
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        d_g = genrou.derivatives()
        d_a = avr.derivatives()
        return np.concatenate([
            np.array([d_g[k] for k in genrou.state_keys]),
            np.array([d_a[k] for k in avr.state_keys]),
        ])

    def _eval_outputs(x):
        _set_combined(x)
        V = _solve_network_with_genrou(x[:n_g], genrou, network)
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()
        out = genrou.algebraic_output()
        avr_out = avr.algebraic_output()
        return np.array([out["|V|"], out["P"], out["Q"],
                         out["Id"], out["Iq"], out["Pe"],
                         genrou.state["delta"], genrou.state["omega"],
                         genrou.state["Eqp"], avr_out["Efd"]])

    def _eval_derivs_u(u):
        avr.params["Vref"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        avr.inputs["Vpss"] = float(u[2])
        return _eval_derivs(x0)

    def _eval_outputs_u(u):
        avr.params["Vref"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        avr.inputs["Vpss"] = float(u[2])
        return _eval_outputs(x0)

    A = numerical_jacobian(_eval_derivs, x0, eps)
    C = numerical_jacobian(_eval_outputs, x0, eps)
    # Restore before u-Jacobians.
    avr.params["Vref"] = Vref0
    genrou.params["Pm"] = Pm0
    avr.inputs["Vpss"] = Vpss0
    u0 = np.array([Vref0, Pm0, Vpss0])
    B = numerical_jacobian(_eval_derivs_u, u0, eps)
    D = numerical_jacobian(_eval_outputs_u, u0, eps)
    # Restore once more.
    avr.params["Vref"] = Vref0
    genrou.params["Pm"] = Pm0
    avr.inputs["Vpss"] = Vpss0
    _set_combined(x0)

    return StateSpace(
        A=A, B=B, C=C, D=D,
        state_names=("delta", "omega", "Eqp", "Edp", "Vc", "x_LL"),
        input_names=("Vref", "Pm", "Vpss"),
        output_names=("|V|", "P", "Q", "Id", "Iq", "Te",
                      "delta", "omega", "Eqp", "Efd"),
    )


# =====================================================================
# Closed-loop linearisation: GENROU + ST1A + PSS1A
# =====================================================================

def linearise_genrou_avr_pss(genrou, avr, pss, network,
                             V_op: complex, S_op: complex,
                             eps: float = 1e-6) -> StateSpace:
    """Linearise GENROU + ST1A + PSS1A around (V_op, S_op).

    State x = (delta, omega, Eqp, Edp, Vc, x_LL, x_w, x_LL1, x_LL2).
    Inputs u = (Vref, Pm, Vpss_ext).  Vpss_ext is an *external* Vpss
    injection on top of whatever the PSS itself produces, useful for
    perturbing the PSS path.
    Outputs as before plus Vpss.
    """
    from .simulator import _solve_network_with_genrou
    genrou.initialise(V_op, S_op)
    avr.initialise(V_op, S_op, Efd_init=genrou.params["Efd"])
    pss.initialise()

    n_g = len(genrou.state_keys)
    n_a = len(avr.state_keys)
    n_p = len(pss.state_keys)

    def _combine():
        return np.concatenate([genrou.flatten(), avr.flatten(), pss.flatten()])

    def _set_combined(x):
        genrou.unflatten(x[:n_g])
        avr.unflatten(x[n_g:n_g + n_a])
        pss.unflatten(x[n_g + n_a:])

    x0 = _combine()
    Vref0 = avr.params["Vref"]
    Pm0 = genrou.params["Pm"]
    Vpss_ext0 = 0.0  # external probe

    def _eval_derivs(x):
        _set_combined(x)
        V = _solve_network_with_genrou(x[:n_g], genrou, network)
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        pss.derivatives()
        Vpss_total = pss.algebraic_output()["Vpss"] + Vpss_ext0
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = Vpss_total
        avr.derivatives()
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

    def _eval_outputs(x):
        _set_combined(x)
        V = _solve_network_with_genrou(x[:n_g], genrou, network)
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        pss.derivatives()
        Vpss_total = pss.algebraic_output()["Vpss"] + Vpss_ext0
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = Vpss_total
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()
        out = genrou.algebraic_output()
        avr_out = avr.algebraic_output()
        pss_out = pss.algebraic_output()
        return np.array([out["|V|"], out["P"], out["Q"],
                         out["Id"], out["Iq"], out["Pe"],
                         genrou.state["delta"], genrou.state["omega"],
                         genrou.state["Eqp"], avr_out["Efd"],
                         pss_out["Vpss"]])

    # Use a non-local for Vpss_ext so the closures pick up changes.
    container = {"Vpss_ext": 0.0}
    def _eval_derivs_with_u(u):
        avr.params["Vref"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        container["Vpss_ext"] = float(u[2])
        # re-eval but with the external Vpss added
        _set_combined(x0)
        V = _solve_network_with_genrou(x0[:n_g], genrou, network)
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        pss.derivatives()
        Vpss_total = pss.algebraic_output()["Vpss"] + container["Vpss_ext"]
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = Vpss_total
        avr.derivatives()
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

    def _eval_outputs_with_u(u):
        avr.params["Vref"] = float(u[0])
        genrou.params["Pm"] = float(u[1])
        container["Vpss_ext"] = float(u[2])
        _set_combined(x0)
        V = _solve_network_with_genrou(x0[:n_g], genrou, network)
        pss.inputs["Delta_omega"] = genrou.state["omega"]
        pss.derivatives()
        Vpss_total = pss.algebraic_output()["Vpss"] + container["Vpss_ext"]
        avr.inputs["V_terminal_mag"] = abs(V)
        avr.inputs["Vpss"] = Vpss_total
        avr.derivatives()
        genrou.params["Efd"] = avr.algebraic_output()["Efd"]
        genrou.inputs["V_terminal"] = V
        genrou.derivatives()
        out = genrou.algebraic_output()
        avr_out = avr.algebraic_output()
        pss_out = pss.algebraic_output()
        return np.array([out["|V|"], out["P"], out["Q"],
                         out["Id"], out["Iq"], out["Pe"],
                         genrou.state["delta"], genrou.state["omega"],
                         genrou.state["Eqp"], avr_out["Efd"],
                         pss_out["Vpss"]])

    A = numerical_jacobian(_eval_derivs, x0, eps)
    C = numerical_jacobian(_eval_outputs, x0, eps)
    avr.params["Vref"] = Vref0
    genrou.params["Pm"] = Pm0
    container["Vpss_ext"] = 0.0
    u0 = np.array([Vref0, Pm0, 0.0])
    B = numerical_jacobian(_eval_derivs_with_u, u0, eps)
    D = numerical_jacobian(_eval_outputs_with_u, u0, eps)
    # Restore
    avr.params["Vref"] = Vref0
    genrou.params["Pm"] = Pm0
    _set_combined(x0)

    return StateSpace(
        A=A, B=B, C=C, D=D,
        state_names=("delta", "omega", "Eqp", "Edp", "Vc", "x_LL",
                     "x_w", "x_LL1", "x_LL2"),
        input_names=("Vref", "Pm", "Vpss_ext"),
        output_names=("|V|", "P", "Q", "Id", "Iq", "Te",
                      "delta", "omega", "Eqp", "Efd", "Vpss"),
    )


# =====================================================================
# Transfer-function extraction
# =====================================================================

def siso_tf(ss: StateSpace, input_name: str, output_name: str):
    """Pick a SISO transfer function from a MIMO state-space.

    Returns (num, den) polynomial coefficients (highest power first),
    suitable for scipy.signal.TransferFunction / Bode / Nyquist.
    """
    j = ss.input_names.index(input_name)
    i = ss.output_names.index(output_name)
    Bi = ss.B[:, j:j+1]
    Ci = ss.C[i:i+1, :]
    Di = ss.D[i:i+1, j:j+1]
    sys = sig.StateSpace(ss.A, Bi, Ci, Di)
    tf = sys.to_tf()
    return np.atleast_1d(np.asarray(tf.num).flatten()), \
           np.atleast_1d(np.asarray(tf.den).flatten())


def series(tf1, tf2):
    """Series composition of two TFs:  num = num1·num2,  den = den1·den2."""
    n1, d1 = tf1
    n2, d2 = tf2
    return np.polymul(n1, n2), np.polymul(d1, d2)


# =====================================================================
# Controller transfer functions
# =====================================================================

def avr_st1a_tf(Ka: float, Tr: float, Tb: float, Tc: float):
    """ST1A regulator TF (simplified, no Ta or rate feedback):

        C_AVR(s) = K_a · (1 + s·T_c) / [ (1 + s·T_b) · (1 + s·T_r) ]
    """
    num = np.array([Ka * Tc, Ka])
    den = np.polymul([Tb, 1.0], [Tr, 1.0])
    return num, den


def pss1a_tf(Ks: float, Tw: float, T1: float, T2: float, T3: float, T4: float):
    """PSS1A TF:

        C_PSS(s) = K_s · (s·T_w / (1+s·T_w))
                       · ((1+s·T_1) / (1+s·T_2))
                       · ((1+s·T_3) / (1+s·T_4))
    """
    num_w = np.array([Tw, 0.0])
    den_w = np.array([Tw, 1.0])
    num_l1 = np.array([T1, 1.0])
    den_l1 = np.array([T2, 1.0])
    num_l2 = np.array([T3, 1.0])
    den_l2 = np.array([T4, 1.0])
    num = Ks * np.polymul(np.polymul(num_w, num_l1), num_l2)
    den = np.polymul(np.polymul(den_w, den_l1), den_l2)
    return num, den


# =====================================================================
# Bode, Nyquist, margins
# =====================================================================

def bode(num, den, w=None):
    """Bode plot magnitude (dB) and phase (deg) for L(s) = num(s)/den(s).

    If ``w`` is not given, uses a log-spaced grid 0.01..1000 rad/s.
    """
    if w is None:
        w = np.logspace(-2, 3, 2000)
    sys = sig.TransferFunction(num, den)
    _, mag_db, phase_deg = sig.bode(sys, w=w)
    return w, mag_db, phase_deg


def nyquist(num, den, w=None):
    """Nyquist locus L(jω) on the complex plane.

    Returns (w, L_complex) where L_complex is an array of complex
    numbers.  Plot ``L_complex.real`` vs ``L_complex.imag`` to see
    the locus; mark −1+0j to apply the encirclement criterion.
    """
    if w is None:
        w = np.logspace(-2, 3, 4000)
    s = 1j * w
    L = np.polyval(num, s) / np.polyval(den, s)
    return w, L


def stability_margins(num, den, w=None) -> dict:
    """Phase and gain margins, plus crossover frequencies.

    Returns a dict with keys ``pm_deg``, ``gm_db``, ``w_gc``,
    ``w_pc``.  Any of them may be ``None`` if the corresponding
    crossover doesn't exist in the swept range.

    Phase margin (PM):  at the gain-crossover frequency ω_gc where
      |L(jω_gc)| = 1, PM = 180° + ∠L(jω_gc).  Positive PM means
      stable, larger means more damped.

    Gain margin (GM):  at the phase-crossover frequency ω_pc where
      ∠L(jω_pc) = -180°, GM = -|L(jω_pc)|_dB.  Positive GM in dB
      means stable.
    """
    if w is None:
        w = np.logspace(-3, 4, 8000)
    sys = sig.TransferFunction(num, den)
    _, mag_db, phase_deg = sig.bode(sys, w=w)
    mag_lin = 10 ** (mag_db / 20.0)
    phase_rad = np.unwrap(np.deg2rad(phase_deg))

    out = {"pm_deg": None, "gm_db": None, "w_gc": None, "w_pc": None}

    # Gain crossover: |L| crosses 1 from above.
    diff_mag = mag_lin - 1.0
    sgn = np.sign(diff_mag)
    cross_g = np.where(np.diff(sgn) != 0)[0]
    for i in cross_g:
        # Linear interpolate.
        a, b = diff_mag[i], diff_mag[i + 1]
        if a == b:
            continue
        f = a / (a - b)
        w_gc = w[i] + f * (w[i + 1] - w[i])
        ph = phase_rad[i] + f * (phase_rad[i + 1] - phase_rad[i])
        pm = 180.0 + np.rad2deg(ph) % 360.0 - 360.0
        # Normalise to (-180, 180]
        while pm > 180.0:
            pm -= 360.0
        while pm <= -180.0:
            pm += 360.0
        out["w_gc"], out["pm_deg"] = float(w_gc), float(pm)
        break  # first crossing only

    # Phase crossover: phase crosses -π (from above-towards or from-below).
    diff_ph = phase_rad - (-np.pi)
    cross_p = np.where(np.diff(np.sign(diff_ph)) != 0)[0]
    for i in cross_p:
        a, b = diff_ph[i], diff_ph[i + 1]
        if a == b:
            continue
        f = a / (a - b)
        w_pc = w[i] + f * (w[i + 1] - w[i])
        mag_at = mag_db[i] + f * (mag_db[i + 1] - mag_db[i])
        out["w_pc"], out["gm_db"] = float(w_pc), float(-mag_at)
        break

    return out


# =====================================================================
# Eigenmode analysis
# =====================================================================

@dataclass
class Mode:
    eigenvalue: complex
    natural_frequency_hz: float
    damping_ratio: float
    dominant_state: str


def eigenmodes(ss: StateSpace) -> list[Mode]:
    """Compute eigenvalues of A and tag each with damping ratio,
    natural frequency (Hz), and the state that participates most
    strongly (modulus of the corresponding right eigenvector entry).
    """
    eigs, V = np.linalg.eig(ss.A)
    modes: list[Mode] = []
    for k, lam in enumerate(eigs):
        wn = abs(lam)
        if wn > 1e-12:
            zeta = float(-lam.real / wn)
        else:
            zeta = float("nan")
        # Dominant state = argmax of |right eigenvector component|
        v = V[:, k]
        j = int(np.argmax(np.abs(v)))
        modes.append(Mode(
            eigenvalue=complex(lam),
            natural_frequency_hz=float(wn / (2 * np.pi)),
            damping_ratio=zeta,
            dominant_state=ss.state_names[j],
        ))
    return modes


__all__ = [
    "StateSpace", "Mode",
    "numerical_jacobian",
    "linearise_genrou", "linearise_genrou_avr", "linearise_genrou_avr_pss",
    "siso_tf", "series",
    "avr_st1a_tf", "pss1a_tf",
    "bode", "nyquist", "stability_margins",
    "eigenmodes",
]
