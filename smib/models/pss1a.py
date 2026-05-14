"""PSS1A — IEEE 421.5-2016 §6.1 single-input power system stabiliser.

The PSS sits in parallel to the AVR's voltage error loop.  It senses
the rotor speed deviation (Δω = ω − ω0 = ω0 · ω̄), processes it through
a washout filter and two stages of lead-lag phase compensation, scales
by a gain Ks, clamps, and injects the result as Vpss into the ST1A
summer.

Topology
--------

    Δω  →  [ s·Tw / (1 + s·Tw) ]  →  x_w     # washout, removes DC
        ↓
        [ (1 + s·T1) / (1 + s·T2) ]  →  x_LL1  # first lead-lag (phase)
        ↓
        [ (1 + s·T3) / (1 + s·T4) ]  →  x_LL2  # second lead-lag (phase)
        ↓
        [ × Ks ]                              # gain
        ↓
        [ clamp Vpssmin..Vpssmax ]            # signal-level limiter
        ↓
        Vpss  →  into the ST1A summer

States (3 total)
----------------

    x_w     [pu]   washout filter internal state
    x_LL1   [pu]   first lead-lag internal state
    x_LL2   [pu]   second lead-lag internal state

The washout has the canonical form

    Tw · dx_w/dt  =  Δω − x_w
    y_w           =  Δω − x_w

so at DC y_w → 0 and at high frequency y_w → Δω.  Tw is chosen large
enough (typically 1–10 s) that the swing mode at ~1 Hz passes through
essentially unattenuated.  In Laplace this realises sTw/(1+sTw).

Each lead-lag block uses the textbook one-state form

    T2 · dx/dt  =  u − x
    y           =  (T1/T2) · u + (1 − T1/T2) · x

with input u (the previous block's output) and internal state x.

Differential equations
----------------------

::

    Tw · dx_w/dt    =  Δω − x_w                        (washout state)
    T2 · dx_LL1/dt  =  y_w − x_LL1                     (first LL state)
    T4 · dx_LL2/dt  =  y_LL1 − x_LL2                   (second LL state)

with the algebraic chain

::

    y_w     =  Δω − x_w
    y_LL1   =  (T1/T2) · y_w   + (1 − T1/T2) · x_LL1
    y_LL2   =  (T3/T4) · y_LL1 + (1 − T3/T4) · x_LL2
    Vpss    =  clamp( Ks · y_LL2,  Vpssmin..Vpssmax )

Inputs
------

    Delta_omega   [pu]   per-unit rotor speed deviation (ω̄ in smib's
                          GENROU; the simulator harness passes this in)

Outputs
-------

    Vpss          [pu]   stabiliser output, fed into ST1A.inputs["Vpss"]
    y_w, y_LL1    [pu]   intermediate signals, exposed for plot debugging

Initialisation
--------------

At steady state Δω = 0, so every state and every internal signal is
zero:

    x_w = 0,  x_LL1 = 0,  x_LL2 = 0,  Vpss = 0.

This means the PSS has *zero contribution* at the operating point —
exactly what we want.  It only activates when the rotor deviates from
synchronous speed, which is the textbook design intent.

Tuning notes (defaults match this SMIB at the canonical operating point)
-----------------------------------------------------------------------

The PSS is meant to add damping torque at the rotor swing frequency.
For SMIB + GENROU + ST1A at the canonical operating point the swing
mode lives at ~0.7 Hz (ωn ≈ 4.4 rad/s).  We design the lead-lag pair
to provide about +90° of phase boost at that frequency, which
compensates the -90° phase lag through the AVR + rotor flux path so
that the resulting Vpss → Te contribution lands in phase with Δω
(pure damping torque).

Default values (typical small-signal-stable thermal unit on SMIB):

    Tw      = 5.0 s         washout
    T1 = T3 = 0.50 s        lead zeros
    T2 = T4 = 0.05 s        lag poles
    Ks      = 20            gain
    Vpssmax = +0.10 pu      signal limit (typical industrial choice)
    Vpssmin = -0.10 pu

Each lead-lag (T1/T2 = 10) contributes up to ~55° phase boost at
ω = 1/sqrt(T1·T2) ≈ 6.3 rad/s.  Two in series give ~95° at swing
frequency — about right for the AVR-driven -90° we need to offset.

Reference: IEEE 421.5-2016 §6.1, Kundur "Power System Stability and
Control" §17.4-17.5 for the phase-compensation derivation.
"""
from __future__ import annotations

from .base import Model


class PSS1A(Model):
    """Single-input PSS with washout + two lead-lag stages.

    Three states: (x_w, x_LL1, x_LL2).  Inputs: Delta_omega (rotor slip
    in pu).  Outputs: Vpss to the AVR summer.
    """

    name = "PSS1A"
    state_keys = ("x_w", "x_LL1", "x_LL2")

    def __init__(self,
                 name: str = "PSS1A",
                 Tw: float = 5.0,
                 T1: float = 0.50,
                 T2: float = 0.05,
                 T3: float = 0.50,
                 T4: float = 0.05,
                 Ks: float = 20.0,
                 Vpssmax: float = 0.10,
                 Vpssmin: float = -0.10):
        params = {
            "Tw": Tw, "T1": T1, "T2": T2, "T3": T3, "T4": T4,
            "Ks": Ks, "Vpssmax": Vpssmax, "Vpssmin": Vpssmin,
        }
        super().__init__(name, params)
        self.inputs = {"Delta_omega": 0.0}
        # Cached intermediate signals for plotting.
        self._y_w: float = 0.0
        self._y_LL1: float = 0.0
        self._y_LL2: float = 0.0
        self._Vpss: float = 0.0

    # ----- helpers -----------------------------------------------------

    def _clamp(self, x: float) -> float:
        return max(self.params["Vpssmin"], min(self.params["Vpssmax"], x))

    # ----- Model interface --------------------------------------------

    def initialise(self, V: complex = 1+0j, S: complex = 0j, **kwargs) -> None:
        """At steady state every state is zero (Δω = 0)."""
        self.state["x_w"] = 0.0
        self.state["x_LL1"] = 0.0
        self.state["x_LL2"] = 0.0
        self.inputs["Delta_omega"] = 0.0
        self._y_w = 0.0
        self._y_LL1 = 0.0
        self._y_LL2 = 0.0
        self._Vpss = 0.0

    def derivatives(self) -> dict:
        Tw = self.params["Tw"]
        T1, T2 = self.params["T1"], self.params["T2"]
        T3, T4 = self.params["T3"], self.params["T4"]
        Ks = self.params["Ks"]

        dom = float(self.inputs["Delta_omega"])
        x_w = self.state["x_w"]
        x_LL1 = self.state["x_LL1"]
        x_LL2 = self.state["x_LL2"]

        # Algebraic chain (each signal is the input to the next block).
        y_w = dom - x_w                                       # washout output
        y_LL1 = (T1 / T2) * y_w + (1.0 - T1 / T2) * x_LL1     # first lead-lag
        y_LL2 = (T3 / T4) * y_LL1 + (1.0 - T3 / T4) * x_LL2   # second lead-lag
        Vpss = self._clamp(Ks * y_LL2)

        # Differential states.
        dxw_dt = (dom - x_w) / Tw
        dxLL1_dt = (y_w - x_LL1) / T2
        dxLL2_dt = (y_LL1 - x_LL2) / T4

        # Cache outputs for plotting / downstream consumption.
        self._y_w = y_w
        self._y_LL1 = y_LL1
        self._y_LL2 = y_LL2
        self._Vpss = Vpss

        return {"x_w": dxw_dt, "x_LL1": dxLL1_dt, "x_LL2": dxLL2_dt}

    def current_injection(self, V: complex) -> complex:
        """The PSS is a measurement/control block — no current contribution."""
        return 0j

    def algebraic_output(self) -> dict:
        return {
            "Vpss": float(self._Vpss),
            "pss_y_w": float(self._y_w),
            "pss_y_LL1": float(self._y_LL1),
            "pss_y_LL2": float(self._y_LL2),
        }
