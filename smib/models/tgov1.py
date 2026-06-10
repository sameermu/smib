"""TGOV1 — simplest steam turbine governor (PSSE / IEEE PES-TR1 §2.2).

The governor is the machine's *frequency* controller, completing the
classical plant stack:

    GENROU  — inertia + rotor flux          (Phase 2.0)
    ST1A    — terminal voltage regulation   (Phase 2.1)
    PSS1A   — swing-mode damping torque     (Phase 2.2)
    TGOV1   — primary frequency response    (Phase 2.3, this file)

When rotor speed drops below synchronous, TGOV1 opens the steam valve
and pushes more mechanical power into the shaft.  The defining
characteristic is the **droop** R: in steady state,

    ΔPm = -Δω / R        (pu power per pu speed deviation)

With R = 0.05 (5 %), a 1 % under-frequency event produces a 20 %
rise in mechanical power — the textbook primary-frequency-response
number every grid code quotes.

Topology (PSSE TGOV1, IEEE PES-TR1 "Dynamic Models for Turbine-
Governors in Power System Studies", 2013, §2.2)
------------------------------------------------------------------

                 Pref
                  |
                  v   +
    Δω ─[× 1/R]─>(−)──> [ 1/(1+sT1) ] ──> x_v  (valve, clamp Vmax/Vmin)
                                            |
                                            v
                               [ (1+sT2)/(1+sT3) ]  (reheat lead-lag)
                                            |
                                            v        −
                                          (sum) <──────[× Dt]── Δω
                                            |
                                            v
                                           Pm  →  GENROU shaft

States (2 total)
----------------
    x_v   [pu]   valve / gate servo position (output of the T1 lag,
                 clamped to Vmax/Vmin)
    x_r   [pu]   reheat lead-lag internal state (one-state realisation,
                 same form as ST1A / PSS1A lead-lags)

Differential equations
----------------------
::

    T1 · dx_v/dt = (Pref - Δω/R) - x_v     (valve servo; clamp + anti-windup)
    T3 · dx_r/dt = x_v - x_r               (reheat)

with the algebraic chain::

    y_reheat = (T2/T3) · x_v + (1 - T2/T3) · x_r
    Pm       = y_reheat - Dt · Δω

T1 (~0.5 s) is the speed-relay / servo lag — how fast the valve
physically moves.  T3 (~7.5 s) is the reheat time constant: steam that
has passed the high-pressure turbine must be re-heated before the
intermediate/low-pressure stages can use it, so most of the power
change arrives on a ~7.5 s first-order.  T2/T3 (~1/3) is the HP power
fraction: the high-pressure stage acts immediately on valve motion,
which the lead-lag zero captures.  Dt is a turbine damping term
(usually 0).

Valve limits and anti-windup
----------------------------
Vmax (rated valve opening, 1.0 pu) and Vmin (0.0 — steam valves do
not go negative) clamp x_v as a *non-windup* limit: when x_v sits on
a limit and the servo pushes further into it, dx_v/dt is forced to
zero so the state cannot wind up beyond the physical gate position.
This mirrors how PSSE implements the TGOV1 valve limit.

Inputs
------
    Delta_omega   [pu]   rotor speed deviation ω̄ (GENROU's "omega" state)

Outputs (algebraic_output)
--------------------------
    Pm        [pu]   mechanical power to the shaft (→ GENROU.params.Pm)
    gov_x_v   [pu]   valve position (for plotting the servo motion)
    gov_Pref  [pu]   the reference (constant after init)

Initialisation
--------------
At steady state Δω = 0, so the whole chain passes Pref through:

    x_v  = Pref
    x_r  = Pref
    Pm   = Pref

and DAE consistency with the machine requires Pref = Pm0, the
mechanical power GENROU's own init derived (Pm0 = Pe0 at steady
state).  This is the "governor Pref" pitfall from the master plan:
the droop gain 1/R acts ONLY on Δω, never on Pref — wiring Pref
through the 1/R gain puts the whole loop out by a factor of 20.

Default parameters (typical reheat steam unit, PES-TR1 §2.2)
------------------------------------------------------------
    R    = 0.05 pu    5 % droop (the near-universal grid-code value)
    T1   = 0.5  s     valve servo lag
    T2   = 2.5  s     reheat lead (HP fraction T2/T3 = 1/3)
    T3   = 7.5  s     reheat time constant
    Vmax = 1.0  pu    valve fully open at rated steam flow
    Vmin = 0.0  pu    valve closed
    Dt   = 0.0  pu    turbine damping (off)

References: IEEE PES-TR1 (2013) §2.2; PSSE Model Library "TGOV1";
Kundur, *Power System Stability and Control* (1994) §11.1 for droop.
"""
from __future__ import annotations

from .base import Model


class TGOV1(Model):
    """Steam turbine governor with droop, valve servo, and reheat.

    Two states: x_v (valve position), x_r (reheat lead-lag internal).
    Input: Delta_omega (rotor slip, pu).  Output: Pm to the shaft.
    """

    name = "TGOV1"
    state_keys = ("x_v", "x_r")

    def __init__(self,
                 name: str = "TGOV1",
                 R: float = 0.05,      # droop [pu speed / pu power]
                 T1: float = 0.5,      # valve servo time constant [s]
                 T2: float = 2.5,      # reheat lead (HP turbine fraction) [s]
                 T3: float = 7.5,      # reheat lag (reheat time constant) [s]
                 Vmax: float = 1.0,    # valve max opening [pu]
                 Vmin: float = 0.0,    # valve min opening [pu]
                 Dt: float = 0.0):     # turbine damping coefficient [pu]
        params = {
            "R": R, "T1": T1, "T2": T2, "T3": T3,
            "Vmax": Vmax, "Vmin": Vmin, "Dt": Dt,
            "Pref": 0.0,   # filled at init (= Pm0 from the machine)
        }
        super().__init__(name, params)
        self.inputs = {"Delta_omega": 0.0}
        self._last_Pm: float = 0.0

    # ----- helpers -----------------------------------------------------

    def _clamp_valve(self, x: float) -> float:
        return max(self.params["Vmin"], min(self.params["Vmax"], x))

    # ----- Model interface --------------------------------------------

    def initialise(self, V: complex = 1 + 0j, S: complex = 0j,
                   Pm_init: float = 0.0, **kwargs) -> None:
        """Back-calculate Pref and both states so dx/dt = 0 at t=0.

        ``Pm_init`` is the mechanical power demanded by the machine at
        steady state (GENROU's params['Pm'] after its initialise()).
        At SS the whole governor chain is a unity pass-through of Pref,
        so every internal quantity equals Pm_init.
        """
        Pref = float(Pm_init)
        if not (self.params["Vmin"] <= Pref <= self.params["Vmax"]):
            raise ValueError(
                f"TGOV1 init: Pref = {Pref:.4f} pu outside valve range "
                f"[{self.params['Vmin']}, {self.params['Vmax']}] — the "
                "operating point is not feasible for this governor."
            )
        self.params["Pref"] = Pref
        self.state["x_v"] = Pref
        self.state["x_r"] = Pref
        self.inputs["Delta_omega"] = 0.0
        self._last_Pm = Pref

    def derivatives(self) -> dict:
        R, T1 = self.params["R"], self.params["T1"]
        T2, T3 = self.params["T2"], self.params["T3"]
        Dt = self.params["Dt"]
        Pref = self.params["Pref"]
        Vmax, Vmin = self.params["Vmax"], self.params["Vmin"]

        dom = float(self.inputs["Delta_omega"])
        x_v = self.state["x_v"]
        x_r = self.state["x_r"]

        # 1. Droop summer and valve servo (non-windup limit).
        u_valve = Pref - dom / R
        dxv_dt = (u_valve - x_v) / T1
        # Anti-windup: freeze the servo when sitting on a limit and
        # being pushed further into it.
        if (x_v >= Vmax and dxv_dt > 0.0) or (x_v <= Vmin and dxv_dt < 0.0):
            dxv_dt = 0.0

        # 2. Reheat lead-lag (one-state realisation).
        x_v_lim = self._clamp_valve(x_v)
        y_reheat = (T2 / T3) * x_v_lim + (1.0 - T2 / T3) * x_r
        dxr_dt = (x_v_lim - x_r) / T3

        # 3. Mechanical power output.
        Pm = y_reheat - Dt * dom
        self._last_Pm = float(Pm)

        return {"x_v": dxv_dt, "x_r": dxr_dt}

    def current_injection(self, V: complex) -> complex:
        """The governor is a mechanical control block — no current."""
        return 0j

    def algebraic_output(self) -> dict:
        return {
            "Pm": float(self._last_Pm),
            "gov_x_v": float(self.state["x_v"]),
            "gov_Pref": float(self.params["Pref"]),
        }
