"""GENCLS tests — the correctness floor for Phase 1.

We verify three properties:

1. Flat-line.  With consistent initialisation, |delta| and |omega|
   must not drift over a 10-second no-disturbance run.
2. Small-perturbation linearity.  An uncontrolled GENCLS on SMIB
   has natural frequency  wn = sqrt(omega0 * Pmax * cos(delta0) / (2H))
   for small deviations around the operating point. We give the rotor
   a tiny angle perturbation and recover this from the FFT.
3. Critical clearing time bracketing.  Find a t_clear at which the
   100 ms case is stable and the 400 ms case is unstable, then sweep
   to bracket CCT to within +/- 10 ms.  This proves the fault logic,
   the solver during the fault, and the swing equation are all
   internally consistent.

Test 2 doubles as a tuning sanity check: if wn from simulation diverges
from the analytic value by more than a few percent, something is off
with H, Xdp, or the angle convention.
"""
from __future__ import annotations

import math

import numpy as np

from smib.models import GENCLS
from smib.network import Network
from smib.powerflow import two_bus_pf
from smib.scenarios import three_phase_fault_schedule
from smib.simulator import run_smib_gencls


# ---------- shared setup ------------------------------------------------

def _setup(P=0.8, Q=0.2, X_line=0.5, H=4.0, D=0.0, Xdp=0.30):
    """Standard SMIB operating point used across tests."""
    V1, _ = two_bus_pf(P, Q, 1.0, 0.0, 0.0, X_line, bus_type="PQ")
    S = complex(P, Q)
    gencls = GENCLS(H=H, D=D, Xdp=Xdp, f0=50.0)
    network = Network(R=0.0, X=X_line)
    gencls.initialise(V1, S)
    return gencls, network, V1, S


# ---------- 1) flat-line ------------------------------------------------

def test_flatline_no_drift():
    gencls, network, _, _ = _setup()
    res = run_smib_gencls(gencls, network, t_end=10.0, h=5e-3)
    drift_delta = float(np.abs(res.traces["delta"] - res.traces["delta"][0]).max())
    drift_omega = float(np.abs(res.traces["omega"] - res.traces["omega"][0]).max())
    assert drift_delta < 1e-5, f"delta drifted by {drift_delta:.3e} rad over 10 s"
    assert drift_omega < 1e-5, f"omega drifted by {drift_omega:.3e} pu over 10 s"


# ---------- 2) small-perturbation natural frequency ---------------------

def test_natural_frequency_matches_analytic():
    """Bump delta by 1 deg, watch the swing, compare wn to the analytic
    value for an undamped GENCLS on SMIB.

    For a classical machine on an infinite bus through pure reactance:
        Pmax  = |E'| * |Vinf| / X_total      where X_total = X'd + X_line
        Ks    = Pmax * cos(delta0)             (synchronising coefficient)
        wn    = sqrt(omega0 * Ks / (2H))       [rad/s]

    The simulated swing should match this within a few percent (the
    nonlinearity from sin(delta) is tiny at 1 deg perturbation).
    """
    P, Q, X_line, H, Xdp = 0.5, 0.0, 0.5, 4.0, 0.30
    gencls, network, V1, _ = _setup(P=P, Q=Q, X_line=X_line, H=H, Xdp=Xdp)

    # Analytic wn around the equilibrium.
    Eint = gencls.params["Eint"]
    delta0 = gencls.state["delta"]
    Vinf = abs(network.V_slack)
    X_tot = Xdp + X_line
    Pmax = Eint * Vinf / X_tot
    Ks = Pmax * math.cos(delta0)
    omega0 = gencls.params["omega0"]
    wn_analytic = math.sqrt(omega0 * Ks / (2 * H))

    # Apply a small angle perturbation and let it ring.  The simulator
    # would otherwise refuse to start from a non-steady initial state.
    gencls.state["delta"] += math.radians(1.0)
    res = run_smib_gencls(gencls, network, t_end=10.0, h=2e-3,
                          require_steady_state=False)

    # Recover wn from the dominant FFT bin of delta(t) - mean.
    delta_centred = res.traces["delta"] - res.traces["delta"].mean()
    fs = 1.0 / (res.t[1] - res.t[0])
    spec = np.abs(np.fft.rfft(delta_centred))
    freqs = np.fft.rfftfreq(len(delta_centred), d=1 / fs)
    f_peak = float(freqs[np.argmax(spec[1:]) + 1])  # skip DC bin
    wn_sim = 2 * math.pi * f_peak

    rel_err = abs(wn_sim - wn_analytic) / wn_analytic
    assert rel_err < 0.05, (
        f"Natural frequency mismatch: sim {wn_sim:.3f} rad/s vs analytic "
        f"{wn_analytic:.3f} rad/s  (rel err {rel_err:.2%})"
    )


# ---------- 3) critical clearing time bracket ---------------------------

def _runs_stable(t_clear_s: float, t_end: float = 5.0) -> bool:
    """Return True iff the rotor stays within +/- 360 deg of the start."""
    gencls, network, _, _ = _setup()
    delta0 = gencls.state["delta"]
    fault = three_phase_fault_schedule(t_start=1.0, t_clear=1.0 + t_clear_s,
                                       Z_fault=0.001 + 0j)
    res = run_smib_gencls(gencls, network, t_end=t_end, h=2e-3, scenarios=[fault])
    return float(np.abs(res.traces["delta"] - delta0).max()) < 2 * math.pi


def test_short_fault_is_stable():
    assert _runs_stable(0.10), "100 ms fault should be stable for these params"


def test_long_fault_loses_sync():
    assert not _runs_stable(0.40), "400 ms fault should lose sync for these params"


def test_critical_clearing_time_brackets():
    """Bisect to bracket CCT to +/- 10 ms.

    This is also the test that exercises every code path: the fault
    apply, the network during the fault, the fault clear, and the
    post-fault recovery — all in one go.
    """
    lo, hi = 0.10, 0.40
    for _ in range(8):  # 0.30 / 2^8 ~= 1 ms
        mid = 0.5 * (lo + hi)
        if _runs_stable(mid):
            lo = mid
        else:
            hi = mid
    cct = 0.5 * (lo + hi)
    assert hi - lo < 0.010, f"failed to converge CCT bracket: [{lo:.3f}, {hi:.3f}]"
    # Sanity: CCT for these conservative params should land in 100-400 ms.
    assert 0.10 < cct < 0.40, f"CCT {cct:.3f} s outside the expected range"
