"""Synchronous reference frame PLL (MHI PSCAD topology, phasor-domain).

Reference
---------
Manitoba Hydro International Ltd., "Phase Locked Loop (PLL) Component",
For PSCAD Version 5.0, initial release, 30 Jan 2020. Author: F. S. Wasylik.
Document lives at power_sys/phase_locked_loop_pll.pdf in this workspace.

The MHI block diagram (Figure 1 of the cited document) is:

    Va,Vb,Vc ─► Park(theta) ─► Vd, Vq ─► ATAN2(Vd, Vq) ─► phi_err
              (eq 1)                     (phase vector technique)

    phi_err ─► (subtract offset_angle) ─► PI(Gp, Gi/s) ─► Δω
           ─► (add ω0) ─► frequency limiter ─► 1/s ─► theta

Three details from the MHI document that differ from the textbook SRF-PLL:

1. The phase detector is atan2(Vd, Vq), not a normalised Vq.  atan2 is
   linear over the full ±π range, so large-angle step tests (MHI's
   Example 1 uses a 90° step) converge cleanly.  A sin-detector saturates
   near ±π/2 and its loop gain collapses there.

2. The axis convention in MHI's Park transform is SWAPPED from the usual
   IEEE/PSSE convention:

       Vq  ← cosine-aligned  (top row of eq 1)
       Vd  ← sine-aligned    (bottom row of eq 1)

   That means atan2(Vd, Vq) returns the angle of the complex number
   (Vq + j·Vd), and this angle goes to zero when the PLL is locked
   (Vd → 0, Vq → |V|).  Every comment in this module follows MHI's
   convention; any dq variable crossing the module boundary should be
   clearly labelled to avoid sign bugs.

3. The base frequency ω0 is added *after* the PI, not folded into it.
   So the PI computes Δω (a frequency deviation), and the integrator
   rate is ω0 + Δω.  A frequency limiter clamps |Δω| before the 1/s,
   which is also where anti-windup lives.

Phasor-domain collapse
----------------------
In an RMS/phasor simulator we don't carry abc time-domain signals; we
have a complex terminal phasor V_terminal in the global synchronous
frame.  The Park transform plus atan2 collapses to one line:

    let V_local = V_terminal * exp(-j * theta_pll)
        Vq = Re(V_local) = |V| * cos(Δ)
        Vd = Im(V_local) = |V| * sin(Δ)       # MHI axis convention
        atan2(Vd, Vq) = Δ = angle(V_terminal) - theta_pll

So in code we just compute `wrap(angle(V_terminal) - theta_pll)` and
skip the 3x3 matrix entirely.  The magnitude |V| cancels, which
automatically makes the detector robust at low voltage (no hand-rolled
`/|V|` guard needed).

Caveat — negative sequence
--------------------------
MHI Example 3 shows a 120 Hz ripple on tracked frequency during an
unbalanced fault, caused by negative-sequence content in the abc
signals.  Positive-sequence phasor simulation does not represent that
ripple.  The PLL will look cleaner in our runs than in a PSCAD EMT run
of the same case.  This is a limitation of RMS, not a modelling bug.

State and parameters
--------------------
States : theta_pll  [rad]   integrator output, wraps to (-π, π]
         x_I        [rad/s] PI integrator state (contribution to Δω)

Params : Kp                 [rad/s per rad of phi_err]
         Ki                 [rad/s² per rad of phi_err]
         omega0             [rad/s]  nominal grid frequency (2π·50 or 2π·60)
         omega_max_dev      [rad/s]  symmetric limiter on Δω
         offset_angle       [rad]   reference phase offset

Inputs : V_terminal         complex phasor at PLL bus (pu, global DQ frame)

Presets (see `preset()`)
------------------------
"bulk_AC"  — MHI Model 1 (Table 1), converted from deg-input to rad-input.
             Kp ≈ 5.24, Ki ≈ 17.45  →  ωn ≈ 4.2 rad/s (~0.67 Hz), ζ ≈ 0.63.
             Slow, well-damped. Matches a typical bulk-AC study PLL.

"ibr_fast" — Textbook IBR-vendor-like PLL, ωn = 2π·30 rad/s, ζ = 1/√2.
             Kp = 2ζωn ≈ 266, Ki = ωn² ≈ 35530.
             Will go unstable on weak grids (SCR ≲ 2-3), which is the
             pedagogical payoff of the Phase 3 SCR sweep.
"""
from __future__ import annotations

import math

import numpy as np

from .base import Model


# ---------- helpers ----------------------------------------------------

def _wrap_pi(x: float) -> float:
    """Wrap an angle to (-π, π]."""
    return (x + math.pi) % (2 * math.pi) - math.pi


def preset(name: str, f0: float = 60.0) -> dict:
    """Return a params dict for a named preset.

    f0 is the nominal grid frequency in Hz.  omega_max_dev defaults to
    2π·5 rad/s (±5 Hz), wide enough for all the scenarios we simulate.
    """
    omega0 = 2 * math.pi * f0
    common = {"omega0": omega0, "omega_max_dev": 2 * math.pi * 5, "offset_angle": 0.0}

    if name == "bulk_AC":
        # MHI Model 1 (Table 1): Kp=300, Ki=1000 in deg-input units.
        # Convert to rad-input: divide by (180/π).
        deg_to_rad = math.pi / 180.0
        return {**common, "Kp": 300.0 * deg_to_rad, "Ki": 1000.0 * deg_to_rad}

    if name == "ibr_fast":
        wn = 2 * math.pi * 30.0  # 30 Hz natural frequency
        zeta = 1.0 / math.sqrt(2.0)
        return {**common, "Kp": 2 * zeta * wn, "Ki": wn * wn}

    raise ValueError(f"Unknown PLL preset: {name!r}. Try 'bulk_AC' or 'ibr_fast'.")


# ---------- model ------------------------------------------------------

class PLL(Model):
    """Phasor-domain SRF-PLL, MHI PSCAD topology.

    Outputs for downstream consumers (IBR controllers, measurement blocks):
        theta_pll : tracked angle (rad, wrapped to (-π, π])
        omega_pll : tracked angular frequency (rad/s)
        phi_err   : phase error at the detector (rad)
    """

    name = "PLL"
    state_keys = ("theta_pll", "x_I")

    def __init__(self, name: str = "PLL", params: dict | str = "bulk_AC"):
        # Allow passing a preset name directly, for notebook readability.
        if isinstance(params, str):
            params = preset(params)
        super().__init__(name, params)
        # Inputs filled each step by the simulator.
        self.inputs = {"V_terminal": complex(1.0, 0.0)}
        # Cached last-computed outputs so plots can grab them without
        # re-running derivatives().
        self._last_phi_err: float = 0.0
        self._last_omega_pll: float = self.params["omega0"]

    # ----- Model interface --------------------------------------------

    def initialise(self, V: complex, S: complex = 0j, **kwargs) -> None:
        """Lock the PLL to the terminal voltage so dx/dt = 0 at t=0.

        The load-flow solution gives V at steady state.  Setting
        theta_pll = angle(V) and x_I = 0 makes phi_err = 0, Δω = 0,
        d(theta_pll)/dt = ω0, d(x_I)/dt = 0.  That is *almost* a
        flat-line — the angle itself rotates at ω0, which is expected
        because the global DQ frame also rotates at ω0.  The check
        `assert_initialised` allows a residual of ω0 on theta_pll;
        see `steady_state_residual()` below.
        """
        self.state["theta_pll"] = float(np.angle(V))
        self.state["x_I"] = 0.0
        self.inputs["V_terminal"] = V

    def derivatives(self) -> dict:
        """Compute state derivatives from current inputs and states."""
        V = self.inputs["V_terminal"]
        Kp = self.params["Kp"]
        Ki = self.params["Ki"]
        w0 = self.params["omega0"]
        w_max = self.params["omega_max_dev"]
        offset = self.params["offset_angle"]

        theta = self.state["theta_pll"]
        x_I = self.state["x_I"]

        # Phase detector (phasor collapse of MHI Park + atan2).
        phi_err = _wrap_pi(float(np.angle(V)) - theta - offset)

        # PI output = proportional + integral state = Δω (unlimited).
        dw_unlim = Kp * phi_err + x_I

        # Frequency limiter.  Total frequency = ω0 + clip(Δω).
        if dw_unlim > w_max:
            dw = w_max
            saturated_positive = True
            saturated_negative = False
        elif dw_unlim < -w_max:
            dw = -w_max
            saturated_positive = False
            saturated_negative = True
        else:
            dw = dw_unlim
            saturated_positive = False
            saturated_negative = False

        omega_pll = w0 + dw

        # Conditional-integration anti-windup: freeze x_I when the
        # limiter is active and phi_err would drive it further into
        # saturation.
        if (saturated_positive and phi_err > 0) or (saturated_negative and phi_err < 0):
            dxI_dt = 0.0
        else:
            dxI_dt = Ki * phi_err

        # Cache for algebraic_output / plotting.
        self._last_phi_err = phi_err
        self._last_omega_pll = omega_pll

        return {"theta_pll": omega_pll, "x_I": dxI_dt}

    def current_injection(self, V: complex) -> complex:
        """A PLL is a measurement block — it does not inject current."""
        return 0j

    def algebraic_output(self) -> dict:
        """Signals the rest of the model stack needs."""
        return {
            "theta_pll": self.state["theta_pll"],
            "omega_pll": self._last_omega_pll,
            "phi_err": self._last_phi_err,
        }

    # ----- diagnostics -------------------------------------------------

    def steady_state_residual(self) -> dict:
        """Residuals that should be zero at lock.

        theta_pll rotates at ω0 in steady state, which is the DQ-frame
        convention we use everywhere.  So the *angle drift* relative to
        ω0 is the thing that should be zero, not d(theta_pll)/dt itself.
        x_I should be flat.
        """
        d = self.derivatives()
        return {
            "theta_pll_vs_w0": d["theta_pll"] - self.params["omega0"],
            "x_I": d["x_I"],
        }

    def natural_frequency(self) -> float:
        """Closed-loop ωn (rad/s) — useful for tuning narratives."""
        return math.sqrt(self.params["Ki"])

    def damping_ratio(self) -> float:
        """Closed-loop ζ — useful for tuning narratives."""
        wn = self.natural_frequency()
        if wn == 0.0:
            return float("inf")
        return self.params["Kp"] / (2.0 * wn)
