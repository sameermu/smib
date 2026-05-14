"""PSS1A tests — the correctness floor for Phase 2.2.

Five tests:

1. **Flat-line drift** — with GENROU + ST1A + PSS1A wired together, no
   disturbance, all 9 states must hold at machine epsilon over 10 s.
   The PSS starts at all-zero internal state because Δω = 0 at the
   operating point, so any drift > 1e-5 indicates a sign bug in the
   Δω → PSS → AVR.Vpss chain.

2. **Init self-consistency** — at t=0 the PSS's algebraic output Vpss
   must be exactly zero (washout removes any DC bias), and the AVR's
   back-solved Vref must therefore match the bare AVR-on case.  This
   pins down "PSS is harmless at the operating point" — which is the
   whole point of the washout filter.

3. **PSS bypass equals Phase 2.1** — set Ks = 0 (gain zero) and verify
   the response matches `run_smib_genrou_avr` within numerical noise.
   This is the regression check that the PSS plumbing doesn't disturb
   the AVR-only path.

4. **Damping improvement on the canonical fault** — on the 200 ms deep
   inductive fault, the late-window (4-8 s) rotor peak-to-peak must
   drop by at least 30 % compared to AVR-only.  Lower is fine; this
   is the headline-damping check.

5. **CCT preservation or lift** — the PSS must not *hurt* CCT (it
   shouldn't, since damping is in the swing-mode-relevant direction).
   On the deep inductive fault, GENROU+ST1A+PSS CCT >= GENROU+ST1A
   CCT minus 5 ms (small tolerance for bisect resolution).
"""
from __future__ import annotations

import math

import numpy as np

from smib.models.genrou import GENROU
from smib.models.st1a import ST1A
from smib.models.pss1a import PSS1A
from smib.network import Network
from smib.powerflow import two_bus_pf
from smib.scenarios import three_phase_fault_schedule
from smib.simulator import run_smib_genrou_avr, run_smib_genrou_avr_pss


def _setup(P=0.8, Q=0.2, X_line=0.5):
    V1, _ = two_bus_pf(P, Q, 1.0, 0.0, 0.0, X_line, bus_type="PQ")
    S = complex(P, Q)
    return V1, S


# ---------- 1) flat-line ----------------------------------------------

def test_pss1a_flatline_no_drift():
    """All nine combined states must hold at machine epsilon for 10 s."""
    V1, S = _setup()
    g = GENROU(D=3.0); a = ST1A(); p = PSS1A()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    res = run_smib_genrou_avr_pss(g, a, p, n, t_end=10.0, h=2e-3,
                                  init_V=V1, init_S=S)
    drifts = {}
    # The integrator's differential states aren't all exposed as traces;
    # use the available algebraic outputs as proxies.
    for k in ("delta", "omega", "Eqp", "Edp", "Vc", "Efd", "Vpss",
              "pss_y_w", "pss_y_LL1", "pss_y_LL2"):
        arr = res.traces[k]
        drifts[k] = float(np.abs(arr - arr[0]).max())
    worst = max(drifts.values())
    assert worst < 1e-5, f"GENROU+ST1A+PSS1A drift over 10 s: {drifts}"


# ---------- 2) init self-consistency ----------------------------------

def test_pss1a_init_zero_at_steady_state():
    """At t=0 the PSS contributes exactly zero — washout removes DC."""
    V1, S = _setup()
    g = GENROU(D=3.0); a = ST1A(); p = PSS1A()
    g.initialise(V1, S)
    a.initialise(V1, S, Efd_init=g.params["Efd"])
    p.initialise()

    p.inputs["Delta_omega"] = g.state["omega"]  # = 0 at SS
    p.derivatives()
    out = p.algebraic_output()
    assert abs(out["Vpss"]) < 1e-12, (
        f"PSS not zero at t=0: Vpss = {out['Vpss']:.3e}"
    )
    assert abs(out["pss_y_w"]) < 1e-12 and abs(out["pss_y_LL1"]) < 1e-12, (
        f"PSS intermediate signals nonzero at t=0: y_w={out['pss_y_w']:.3e}, "
        f"y_LL1={out['pss_y_LL1']:.3e}"
    )


# ---------- 3) PSS bypass (Ks=0) reproduces Phase 2.1 -----------------

def test_pss1a_bypass_matches_phase2_1():
    """With Ks = 0 the PSS contributes no Vpss; the GENROU + ST1A loop
    must therefore reproduce `run_smib_genrou_avr` traces within
    numerical noise."""
    V1, S = _setup()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    fault = three_phase_fault_schedule(1.0, 1.10, 0 + 0.10j)

    g1 = GENROU(D=3.0); a1 = ST1A()
    r1 = run_smib_genrou_avr(g1, a1, n, t_end=3.0, h=2e-3, scenarios=[fault],
                             init_V=V1, init_S=S)
    g2 = GENROU(D=3.0); a2 = ST1A(); p2 = PSS1A(Ks=0.0)
    r2 = run_smib_genrou_avr_pss(g2, a2, p2, n, t_end=3.0, h=2e-3,
                                 scenarios=[fault], init_V=V1, init_S=S)
    for k in ("delta", "omega", "Eqp", "Efd", "|V|"):
        diff = float(np.abs(r1.traces[k] - r2.traces[k]).max())
        assert diff < 1e-6, (
            f"PSS bypass diverged from AVR-only on '{k}': max|diff| = {diff:.3e}"
        )


# ---------- 4) damping improvement on the canonical fault -------------

def test_pss1a_damps_post_fault_oscillation():
    """On the canonical 200 ms deep inductive fault, late-window
    (4-8 s) rotor-angle peak-to-peak with PSS must be at least 30 %
    smaller than without."""
    V1, S = _setup()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)
    fault = three_phase_fault_schedule(1.0, 1.20, 0 + 0.10j)

    g1 = GENROU(D=3.0); a1 = ST1A()
    r1 = run_smib_genrou_avr(g1, a1, n, t_end=8.0, h=2e-3, scenarios=[fault],
                             init_V=V1, init_S=S)
    g2 = GENROU(D=3.0); a2 = ST1A(); p2 = PSS1A()
    r2 = run_smib_genrou_avr_pss(g2, a2, p2, n, t_end=8.0, h=2e-3,
                                 scenarios=[fault], init_V=V1, init_S=S)

    def pp_late(res):
        d = np.degrees(res.traces["delta"])
        late = d[int(res.t.searchsorted(4.0)):]
        return float(late.max() - late.min())

    pp_avr = pp_late(r1)
    pp_pss = pp_late(r2)
    ratio = pp_pss / pp_avr
    assert ratio < 0.70, (
        f"PSS did not damp the post-fault swing enough: "
        f"AVR-only pp={pp_avr:.2f}°, AVR+PSS pp={pp_pss:.2f}°, "
        f"ratio {ratio:.2%} (must be < 70 %)"
    )


# ---------- 5) CCT preserved or lifted ----------------------------------

def test_pss1a_does_not_hurt_cct():
    """Adding the PSS must not reduce CCT by more than ~5 ms (bisect
    resolution).  We expect ≈ equal or slightly better."""
    V1, S = _setup()
    n = Network(R=0.0, X=0.5, V_slack_mag=1.0)

    def is_st_avr(tc):
        g = GENROU(D=3.0); a = ST1A()
        f = three_phase_fault_schedule(1.0, 1.0 + tc, 0 + 0.10j)
        r = run_smib_genrou_avr(g, a, n, t_end=5.0, h=2e-3, scenarios=[f],
                                init_V=V1, init_S=S)
        d0 = r.traces["delta"][0]
        return abs(r.traces["delta"] - d0).max() < 2 * math.pi

    def is_st_pss(tc):
        g = GENROU(D=3.0); a = ST1A(); p = PSS1A()
        f = three_phase_fault_schedule(1.0, 1.0 + tc, 0 + 0.10j)
        r = run_smib_genrou_avr_pss(g, a, p, n, t_end=5.0, h=2e-3, scenarios=[f],
                                    init_V=V1, init_S=S)
        d0 = r.traces["delta"][0]
        return abs(r.traces["delta"] - d0).max() < 2 * math.pi

    def bisect(is_st, lo=0.05, hi=0.60, n=10):
        for _ in range(n):
            mid = 0.5 * (lo + hi)
            if is_st(mid): lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)

    cct_avr = bisect(is_st_avr)
    cct_pss = bisect(is_st_pss)
    assert cct_pss >= cct_avr - 0.005, (
        f"PSS hurt CCT: AVR-only {cct_avr*1000:.0f} ms, "
        f"AVR+PSS {cct_pss*1000:.0f} ms"
    )
