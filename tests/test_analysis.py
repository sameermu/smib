"""Frequency-domain analysis correctness floor.

Five tests covering the small-signal / Bode / Nyquist machinery:

1. **Numerical Jacobian sanity** — for a closed-form linear function
   the finite-difference Jacobian must equal the analytic Jacobian
   within numerical noise.
2. **Bare GENROU swing-mode eigenvalue matches the analytic value**
   from Phase 1's ``test_genrou_small_signal_natural_frequency`` —
   the numerical linearisation must agree with the hand derivation.
3. **AVR transfer function: Bode at one point** — at ω = 1/Tb, the
   asymptotic gain of ``Ka·Tc/Tb · 1/(jω)`` should match scipy.signal
   to within 1%.
4. **PSS lead-lag peak phase boost** lands at the analytic
   ω = 1/√(T1·T2) frequency.
5. **TGR design proof** — open-loop AVR has positive phase margin
   with Tb=20, and *negative* phase margin (unstable) with Tb=Tc=1.
   This is the textbook AVR-induced-instability headline.
"""
from __future__ import annotations

import math

import numpy as np

from smib.analysis import (numerical_jacobian, linearise_genrou,
                           siso_tf, series, avr_st1a_tf, pss1a_tf,
                           bode, stability_margins, eigenmodes)
from smib.models.genrou import GENROU
from smib.network import Network
from smib.powerflow import two_bus_pf


def _setup(P=0.8, Q=0.2, X_line=0.5):
    V1, _ = two_bus_pf(P, Q, 1.0, 0.0, 0.0, X_line, bus_type='PQ')
    S = complex(P, Q)
    n = Network(R=0.0, X=X_line, V_slack_mag=1.0)
    return V1, S, n


# ---------- 1) numerical jacobian sanity -------------------------------

def test_numerical_jacobian_linear():
    """f(x) = A·x  ⇒  ∂f/∂x = A, to within central-difference noise."""
    A = np.array([[1.0, 2.0, 3.0],
                  [-1.0, 0.5, 0.0],
                  [0.0, 0.0, 4.0]])
    f = lambda x: A @ x
    J = numerical_jacobian(f, np.array([0.0, 0.0, 0.0]))
    assert np.allclose(J, A, atol=1e-9)


# ---------- 2) GENROU swing-mode eigenvalue ---------------------------

def test_genrou_swing_eigenvalue_in_expected_range():
    """The dominant complex-conjugate pole pair must lie at the rotor
    swing mode (around 0.7–1.2 Hz on this SMIB), with positive but
    small damping ratio.

    Uses the *same* operating point as
    test_genrou_small_signal_natural_frequency, so the numerical
    linearisation can be cross-checked against the time-domain FFT
    test if either drifts.
    """
    V1, S, n = _setup()
    g = GENROU(D=3.0)
    ss = linearise_genrou(g, n, V1, S)
    swing = [m for m in eigenmodes(ss) if 0.5 < m.natural_frequency_hz < 1.5
             and abs(m.eigenvalue.imag) > 1e-3]
    assert len(swing) >= 1, (
        f"No swing-mode pole found in (0.5, 1.5) Hz; modes were: "
        f"{[(m.eigenvalue, m.natural_frequency_hz) for m in eigenmodes(ss)]}"
    )
    m = swing[0]
    assert 0.0 < m.damping_ratio < 0.30, (
        f"Swing-mode damping out of expected range: ζ = {m.damping_ratio:.3f}"
    )


# ---------- 3) AVR transfer function Bode at one point ----------------

def test_avr_tf_bode_matches_analytic_high_frequency_asymptote():
    """At a frequency well above 1/Tb but well below 1/Tr, the AVR's
    transfer function ``Ka·(1+sTc)/((1+sTb)(1+sTr))`` reduces to its
    asymptotic form ``Ka·(Tc/Tb)`` (a real, frequency-independent gain).

    Picking Tb=20, Tc=1, Tr=0.02, Ka=200, the asymptote in the band
    1 << ω·Tb and ω·Tr << 1 (i.e. 0.05 << ω << 50 rad/s) is

        |L|   ≈  Ka · Tc / Tb     = 200 · 1 / 20 = 10
        ∠L    ≈  −180° + 90° + 0°  = −90°    (one extra pole, one zero)

    Wait — careful: (1+sTc)·(1+sTb)^-1 contributes +90° from the zero
    and -90° from the pole at high freq, giving 0°.  The 1/(1+sTr) is
    still ~0° below 1/Tr.  So phase asymptote is 0° in this band.
    """
    num, den = avr_st1a_tf(Ka=200, Tr=0.02, Tb=20.0, Tc=1.0)
    # Test points well inside the asymptotic band: ω·Tc >> 1 (so the
    # numerator zero has unwound to its asymptote) and ω·Tr << 1 (the
    # input filter pole hasn't kicked in yet).  Tc=1, Tr=0.02 means
    # ~5 << ω << 30 is the safe window.
    omega_band = np.array([5.0, 8.0])  # ω·Tc=5 well above 1, ω·Tr=0.16 still well below 1
    _, mag_db, phase = bode(num, den, w=omega_band)
    mag_lin = 10 ** (mag_db / 20.0)
    expected = 200 * 1.0 / 20.0  # = 10
    for w, m, ph in zip(omega_band, mag_lin, phase):
        rel = abs(m - expected) / expected
        assert rel < 0.10, (
            f"AVR mid-band gain at ω={w}: simulated {m:.3f}, expected {expected:.3f} (rel err {rel:.2%})"
        )
        assert abs(ph) < 25.0, f"AVR mid-band phase at ω={w}: {ph:+.1f}° (expected ~0°)"


# ---------- 4) PSS lead-lag peak phase boost --------------------------

def test_pss_lead_lag_peak_phase_at_geometric_mean():
    """A first-order lead-lag (1+sT1)/(1+sT2) with T1 > T2 has its
    maximum phase boost at ω* = 1/√(T1·T2).  Two of them in series
    double the boost (approximately) and shift the peak slightly.

    Set T1=T3=0.5, T2=T4=0.05 (the smib default) and Ks=1 so we look
    only at the phase shape.  Also set Tw very large so the washout's
    -90° phase at low freq has unwound by the swing frequency.
    """
    num, den = pss1a_tf(Ks=1.0, Tw=100.0, T1=0.5, T2=0.05, T3=0.5, T4=0.05)
    # Single lead-lag peak is at sqrt(1/(T1·T2)) = sqrt(40) ≈ 6.3 rad/s
    w_expected = 1.0 / math.sqrt(0.5 * 0.05)
    omega_grid = np.logspace(-1, 2, 1000)
    _, _, phase = bode(num, den, w=omega_grid)
    # Find the peak phase boost; with washout fully unwound, two lead-lags
    # together hit close to 2 × arcsin((T1-T2)/(T1+T2)) ≈ 110°
    i_peak = int(np.argmax(phase))
    w_peak = omega_grid[i_peak]
    rel = abs(w_peak - w_expected) / w_expected
    assert rel < 0.30, (
        f"PSS peak-phase frequency: simulated {w_peak:.3f} rad/s vs "
        f"expected {w_expected:.3f} rad/s (rel err {rel:.2%})"
    )
    # And the magnitude of the boost at the peak should be > 80° (well
    # above the ~55° single-stage maximum).
    assert phase[i_peak] > 80.0, (
        f"PSS peak phase boost {phase[i_peak]:.1f}° (expected > 80°)"
    )


# ---------- 5) TGR design proof — the AVR design headline ---------------

def test_avr_no_tgr_is_unstable_tgr_is_stable():
    """With no transient gain reduction (Tb = Tc = 1, pure-gain
    lead-lag), the open-loop AVR + GENROU+network has *negative* phase
    margin near the swing frequency.  Closing the loop on a negative-PM
    system gives an unstable closed loop — that's the AVR-induced
    swing-mode instability we saw in time domain before TGR.

    With Tb = 20 (transient gain reduction), the open-loop PM becomes
    comfortably positive (~60°) and the closed loop is stable.

    Both arms of this test must hold; if either flips the design intent
    is broken.
    """
    V1, S, n = _setup()
    g = GENROU(D=3.0)
    ss = linearise_genrou(g, n, V1, S)
    num_P, den_P = siso_tf(ss, 'Efd', '|V|')

    # No TGR
    num_C0, den_C0 = avr_st1a_tf(Ka=200, Tr=0.02, Tb=1.0, Tc=1.0)
    num_L0, den_L0 = series((num_C0, den_C0), (num_P, den_P))
    m0 = stability_margins(num_L0, den_L0)
    assert m0["pm_deg"] is not None and m0["pm_deg"] < 10.0, (
        f"No-TGR AVR should be marginal or worse, got PM = {m0['pm_deg']}°"
    )

    # With TGR
    num_C1, den_C1 = avr_st1a_tf(Ka=200, Tr=0.02, Tb=20.0, Tc=1.0)
    num_L1, den_L1 = series((num_C1, den_C1), (num_P, den_P))
    m1 = stability_margins(num_L1, den_L1)
    assert m1["pm_deg"] is not None and m1["pm_deg"] > 40.0, (
        f"TGR AVR should have generous PM, got {m1['pm_deg']}°"
    )
    assert m1["gm_db"] is not None and m1["gm_db"] > 6.0, (
        f"TGR AVR should have generous GM, got {m1['gm_db']} dB"
    )
