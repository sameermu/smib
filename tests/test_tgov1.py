"""TGOV1 governor tests — the correctness floor for Phase 2.3.

Five tests, mirroring the structure of test_st1a.py and
test_pss1a.py:

1. **Flat-line drift** — with the full GENROU + ST1A + PSS1A + TGOV1
   stack wired together, no disturbance, all 11 states stay at machine
   epsilon over 10 s.  Catches sign bugs in the ω̄ → governor → Pm →
   swing-equation coupling and trace-plumbing bugs in
   ``run_smib_genrou_avr_pss_gov``.

2. **Init self-consistency** — the governor's steady-state Pm output
   must round-trip GENROU's demanded Pm exactly, and the combined
   11-state derivative vector at t=0 must be < 1e-6.  This pins the
   master-plan pitfall: Pref must equal Pm0 directly — the 1/R droop
   gain acts only on Δω, never on Pref (wiring Pref through 1/R puts
   everything out by a factor of 20).

3. **Droop steady state** — sustained grid under-frequency of
   Δf = -0.004 pu (-0.2 Hz on 50 Hz), applied via the infinite-bus
   angle ramp.  After settling, the rotor slip must equal the grid
   deviation and the mechanical power rise must satisfy
   ΔPm = -Δω/R within 2 %.  This is the defining primary-frequency-
   response check — if it fails, the droop gain, the summer sign, or
   the valve chain is wrong.

4. **Governor locked reproduces Phase 2.2** — with droop R huge
   (1/R ≈ 0) the governor holds Pm = Pref regardless of slip, which
   is exactly the constant-Pm boundary condition of the PSS stack.
   CCT on the deep inductive fault must agree with
   run_smib_genrou_avr_pss within a few ms.  Any divergence means
   the governor wiring perturbs the machine even when it should be
   passive.

5. **Valve limit + anti-windup** — a deep under-frequency event that
   demands more than Vmax must clamp the valve at Vmax (Pm ceiling)
   without the state winding up beyond the limit; after the frequency
   recovers, Pm must come back down promptly (no slow unwind from an
   over-charged integrator).
"""
from __future__ import annotations

import math

import numpy as np

from smib.models.genrou import GENROU
from smib.models.pss1a import PSS1A
from smib.models.st1a import ST1A
from smib.models.tgov1 import TGOV1
from smib.network import Network
from smib.powerflow import two_bus_pf
from smib.scenarios import (grid_frequency_ramp_schedule,
                            three_phase_fault_schedule)
from smib.simulator import (run_smib_genrou_avr_pss,
                            run_smib_genrou_avr_pss_gov)


# ---------- shared setup ----------------------------------------------

def _setup(P=0.8, Q=0.2, X_line=0.5):
    V1, _ = two_bus_pf(P, Q, 1.0, 0.0, 0.0, X_line, bus_type="PQ")
    S = complex(P, Q)
    return V1, S


def _fresh_stack(**gov_kwargs):
    return GENROU(D=3.0), ST1A(), PSS1A(), TGOV1(**gov_kwargs)


# ---------- 1) flat-line ----------------------------------------------

def test_tgov1_flatline_no_drift():
    """All eleven combined states must hold at machine epsilon over
    10 s with no disturbance."""
    V1, S = _setup()
    g, a, p, gv = _fresh_stack()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    res = run_smib_genrou_avr_pss_gov(g, a, p, gv, n, t_end=10.0, h=2e-3,
                                      init_V=V1, init_S=S)
    drifts = {}
    for k in ("delta", "omega", "Eqp", "Edp", "Vc", "Efd", "Vpss",
              "Pm", "gov_x_v", "|V|"):
        arr = res.traces[k]
        drifts[k] = float(np.abs(arr - arr[0]).max())
    worst = max(drifts.values())
    assert worst < 1e-5, f"Full-stack drift over 10 s: {drifts}"


# ---------- 2) init self-consistency ----------------------------------

def test_tgov1_init_self_consistent():
    """The governor's SS Pm output must equal GENROU's demanded Pm,
    and the combined 11-state derivative vector must be ~0 at t=0."""
    V1, S = _setup()
    g = GENROU(D=3.0)
    g.initialise(V1, S)
    Pm_demanded = g.params["Pm"]

    gv = TGOV1()
    gv.initialise(Pm_init=Pm_demanded)
    gv.derivatives()
    Pm_supplied = gv.algebraic_output()["Pm"]
    assert abs(Pm_supplied - Pm_demanded) < 1e-12, (
        f"Governor init mismatch: GENROU asked for Pm={Pm_demanded:.6f}, "
        f"TGOV1 returned Pm={Pm_supplied:.6f}"
    )
    # Pref must be Pm0 itself (NOT Pm0·R or Pm0/R — the factor-of-20 trap).
    assert abs(gv.params["Pref"] - Pm_demanded) < 1e-12

    # Combined derivative vector at machine epsilon.
    a = ST1A(); p = PSS1A()
    a.initialise(V1, S, Efd_init=g.params["Efd"])
    p.initialise()
    a.inputs["V_terminal_mag"] = abs(V1)
    g.params["Efd"] = a.algebraic_output()["Efd"]
    g.inputs["V_terminal"] = V1
    worst = max(
        max(abs(v) for v in g.derivatives().values()),
        max(abs(v) for v in a.derivatives().values()),
        max(abs(v) for v in p.derivatives().values()),
        max(abs(v) for v in gv.derivatives().values()),
    )
    assert worst < 1e-6, f"Combined stack not at steady state: {worst:.3e}"


# ---------- 3) droop steady state --------------------------------------

def test_tgov1_droop_steady_state():
    """Sustained -0.004 pu grid frequency deviation → slip locks to the
    grid and ΔPm = -Δω/R within 2 %."""
    V1, S = _setup()
    g, a, p, gv = _fresh_stack()
    R = gv.params["R"]
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    dw = -0.004
    ev = grid_frequency_ramp_schedule(t_start=1.0, delta_f_pu=dw,
                                      ramp_time=2.0)
    res = run_smib_genrou_avr_pss_gov(g, a, p, gv, n, t_end=45.0, h=2e-3,
                                      scenarios=[ev], init_V=V1, init_S=S)
    slip_final = res.traces["omega"][-1]
    dPm = res.traces["Pm"][-1] - res.traces["Pm"][0]
    dPm_expected = -dw / R  # = 0.08 pu for R = 5 %

    assert abs(slip_final - dw) < 1e-4, (
        f"Rotor slip did not lock to grid deviation: {slip_final:+.5f} "
        f"vs {dw:+.5f}"
    )
    assert abs(dPm - dPm_expected) / dPm_expected < 0.02, (
        f"Droop response off: dPm = {dPm:.4f} pu, expected "
        f"-dw/R = {dPm_expected:.4f} pu"
    )


# ---------- 4) governor locked reproduces Phase 2.2 -------------------

def test_tgov1_locked_matches_pss_stack():
    """With R huge the governor is inert (Pm pinned at Pref) — CCT on
    the deep inductive fault must match the Phase 2.2 stack within a
    few ms."""
    V1, S = _setup()

    def is_stable_pss(t_clear):
        g = GENROU(D=3.0); a = ST1A(); p = PSS1A()
        n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
        f = three_phase_fault_schedule(1.0, 1.0 + t_clear, 0 + 0.10j)
        r = run_smib_genrou_avr_pss(g, a, p, n, t_end=5.0, h=2e-3,
                                    scenarios=[f], init_V=V1, init_S=S)
        d0 = r.traces["delta"][0]
        return abs(r.traces["delta"] - d0).max() < 2 * math.pi

    def is_stable_locked_gov(t_clear):
        g, a, p, gv = _fresh_stack(R=1e6)
        n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
        f = three_phase_fault_schedule(1.0, 1.0 + t_clear, 0 + 0.10j)
        r = run_smib_genrou_avr_pss_gov(g, a, p, gv, n, t_end=5.0, h=2e-3,
                                        scenarios=[f], init_V=V1, init_S=S)
        d0 = r.traces["delta"][0]
        return abs(r.traces["delta"] - d0).max() < 2 * math.pi

    def bisect(is_stable, lo=0.05, hi=0.60, n_iter=8):
        for _ in range(n_iter):
            mid = 0.5 * (lo + hi)
            if is_stable(mid):
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    cct_pss = bisect(is_stable_pss)
    cct_locked = bisect(is_stable_locked_gov)
    diff_ms = abs(cct_locked - cct_pss) * 1000
    assert diff_ms < 5.0, (
        f"Gov-locked CCT {cct_locked*1000:.0f} ms diverged from "
        f"Phase 2.2 CCT {cct_pss*1000:.0f} ms by {diff_ms:.1f} ms"
    )


# ---------- 5) valve limit + anti-windup -------------------------------

def test_tgov1_valve_limit_and_antiwindup():
    """A -2 % frequency event demands ΔPm = +0.4 pu — beyond the
    Vmax = 1.0 valve ceiling (Pref = 0.8).  The valve must clamp at
    Vmax, not wind up beyond it; on frequency recovery Pm must fall
    back towards Pref without a long unwind tail."""
    V1, S = _setup()
    g, a, p, gv = _fresh_stack()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    dw = -0.02   # demands dPm = 0.4 pu; headroom is only 0.2 pu

    # Down at t=1 (2 s ramp), recover at t=25 (2 s ramp back to 0).
    down = grid_frequency_ramp_schedule(t_start=1.0, delta_f_pu=dw,
                                        ramp_time=2.0)
    up = grid_frequency_ramp_schedule(t_start=25.0, delta_f_pu=-dw,
                                      ramp_time=2.0)
    res = run_smib_genrou_avr_pss_gov(g, a, p, gv, n, t_end=50.0, h=2e-3,
                                      scenarios=[down, up],
                                      init_V=V1, init_S=S)
    t = res.t
    x_v = res.traces["gov_x_v"]
    Pm = res.traces["Pm"]
    Vmax = gv.params["Vmax"]

    # (a) the valve saturates at Vmax.  The non-windup limit zeroes the
    # servo derivative at the boundary; the trapezoidal corrector can
    # overshoot by O(h·dx/dt) on the crossing step, so allow a hair —
    # what matters is that the state cannot WIND UP (drift further in).
    assert x_v.max() <= Vmax + 1e-3, (
        f"Valve wound up past Vmax: max x_v = {x_v.max():.6f}"
    )
    window = (t > 15.0) & (t < 25.0)
    assert x_v[window].min() > Vmax - 1e-6, (
        "Valve should be pinned at Vmax during the deep event"
    )
    # Pm ceiling: the reheat (T3 = 7.5 s) passes the clamped valve
    # through asymptotically — by the end of the 10 s pinned window
    # Pm must be within 1 % of the Vmax ceiling.
    assert abs(Pm[window][-1] - Vmax) < 1e-2

    # (b) anti-windup: after recovery the valve must be moving back
    # down promptly — within ~3 servo time constants of the recovery
    # completing, x_v must have left the limit by a finite margin.
    t_check = 27.0 + 3.0 * gv.params["T1"]
    idx = int(np.searchsorted(t, t_check))
    assert x_v[idx] < Vmax - 1e-3, (
        f"Valve stuck at limit after recovery: x_v({t_check:.1f} s) = "
        f"{x_v[idx]:.6f} — integrator wind-up suspected"
    )
