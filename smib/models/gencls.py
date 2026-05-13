"""GENCLS — classical synchronous machine.

Constant voltage E' behind transient reactance X'd, swing equation for
rotor dynamics, no field flux dynamics, no saliency, no AVR, no PSS, no
governor.  This is the simplest dynamic generator model and the right
starting point for transient stability — every concept (rotor angle,
inertia, equal-area criterion, critical clearing time) appears here in
its purest form.

Model
-----
Per-unit, common rotor base.  Two states:
    delta : rotor angle relative to the global synchronous DQ frame  [rad]
    omega : per-unit speed deviation (slip), so omega = 0 at synchronous

Swing equation:
    2H * d(omega)/dt = Pm - Pe - D*omega         (rotor dynamics)
    d(delta)/dt      = omega0 * omega            (kinematic relation)

Electrical:
    E' = |E'| * exp(j * delta)         (internal EMF, global DQ frame)
    I  = (E' - V_terminal) / (j * X'd) (current INTO the bus from the machine)
    Pe = Re(V_terminal * conj(I))      (electrical power output)

Initialisation
--------------
From PF: given V (complex terminal phasor) and S = P + jQ (terminal
power injection, generator convention),
    I0     = conj(S / V)
    E0     = V + j * X'd * I0
    delta0 = angle(E0)
    Eint   = |E0|
    Pm0    = Re(V * conj(I0))     # mechanical input balances Pe at steady state
    omega0_state = 0               # slip is zero at synchronous

After this, derivatives() returns dx/dt = 0 to within machine epsilon.
That is the flat-line floor.

Network coupling
----------------
GENCLS is not a pure current source — its injection depends on V.  The
clean way to handle this is to expose the Norton equivalent:
    Y_norton = 1 / (j * X'd)
    I_norton = E' * Y_norton
and let the network solver fold Y_norton into Ybus.  See
`smib.simulator.solve_network_with_machine` for the substitution.

Machine dq quantities (for the Id / Iq canonical traces)
--------------------------------------------------------
PSSE / Kundur convention: q-axis aligned with the rotor EMF E', d-axis
90° behind q (in the direction of rotation, q leads d).  In machine
dq the EMF is purely on the q-axis (E_d = 0, E_q = |E'|), and the
current splits into:

  Iq — torque-producing component, parallel to E'.  Positive at
       steady state for an exporting generator.
  Id — demagnetising/magnetising component, perpendicular to E'.
       Positive for lagging power factor (current vector trailing E').
       During a heavy LV event the rotor flux fights to stay constant
       (constant flux linkage theorem), so a large Id transient
       appears — this is the "natural fault current contribution" of
       a synchronous machine and is what makes generators the dominant
       short-circuit source on the network.

To go from a global-DQ phasor X to its machine dq components in this
convention we apply the rotation `j * exp(-j*delta)` (rotate by
-delta, then by +90°), which lines up the q-axis with the EMF and
puts q on the positive real axis, d on the positive imaginary axis.

Note on LV ride-through: with a physically realistic *inductive*
fault impedance, GENCLS satisfies the LV validation rule on its own
— the constant flux linkage theorem keeps |E'| pinned during the
fault, the resulting current is mostly reactive (q-axis), and at
the terminal Q rises while P drops.  This is the **inherent reactive
support** of synchronous machines.  Caveat: GENCLS holds |E'|
constant *forever* by assumption; in a real machine |E'| decays over
the field-winding time constant T'_d0 ~ 5-10 s, so the inherent
support drops off after ~1 s.  An AVR (Phase 2) sustains and
amplifies this support.

A purely *resistive* shunt fault is a pathological case: it
dissipates real power locally at the gen bus, pulls extra Pe through
the gen terminal, and can flip the sign of terminal Q.  Real fault
paths are predominantly inductive (transformer leakages, ground
impedances), so the deep fault scenario in the Phase 1 notebook uses
Z_f = j*0.10 pu.
"""
from __future__ import annotations

import math

import numpy as np

from .base import Model


class GENCLS(Model):
    """Classical machine: constant E' behind X'd, swing equation only."""

    name = "GENCLS"
    state_keys = ("delta", "omega")

    def __init__(self, name: str = "GENCLS",
                 H: float = 4.0,
                 D: float = 0.0,
                 Xdp: float = 0.30,
                 f0: float = 50.0):
        """
        Parameters
        ----------
        H   : inertia constant on machine MVA base [s]. Typical 3-6 for thermal,
              2-4 for hydro, 0.3-0.6 for IBR-rated synchronous condensers.
        D   : damping coefficient on slip [pu torque per pu slip]. Real machines
              rarely have natural damping > 1; 0 is the textbook GENCLS choice
              and also the worst-case for stability margin demonstrations.
        Xdp : transient reactance X'd [pu]. Typical 0.2-0.4.
        f0  : nominal grid frequency [Hz].
        """
        params = {"H": H, "D": D, "Xdp": Xdp, "omega0": 2 * math.pi * f0}
        # Eint and Pm are state-of-the-machine constants set at initialise().
        params["Eint"] = 1.0
        params["Pm"] = 0.0
        super().__init__(name, params)
        self.inputs = {"V_terminal": complex(1.0, 0.0)}
        # Cached for plotting.
        self._last_Pe: float = 0.0
        self._last_I: complex = 0j

    # ----- initialisation ---------------------------------------------

    def initialise(self, V: complex, S: complex, **kwargs) -> None:
        """Back-calculate delta, |E'|, and Pm so dx/dt = 0 at t=0.

        Generator convention: S = V * conj(I).  So I = conj(S / V).
        Then E' = V + j*X'd*I.  delta = angle(E'), |E'| is its magnitude.
        Pm = Re(V * conj(I)) ensures the swing equation is in balance.
        """
        Xdp = self.params["Xdp"]
        I = np.conj(S / V)
        E = V + 1j * Xdp * I
        self.params["Eint"] = float(abs(E))
        self.params["Pm"] = float((V * np.conj(I)).real)
        self.state["delta"] = float(np.angle(E))
        self.state["omega"] = 0.0
        self.inputs["V_terminal"] = V
        # Cache initial values for diagnostics.
        self._last_I = complex(I)
        self._last_Pe = self.params["Pm"]

    # ----- model interface --------------------------------------------

    def derivatives(self) -> dict:
        V = self.inputs["V_terminal"]
        delta = self.state["delta"]
        omega = self.state["omega"]
        H = self.params["H"]
        D = self.params["D"]
        Xdp = self.params["Xdp"]
        Eint = self.params["Eint"]
        Pm = self.params["Pm"]
        w0 = self.params["omega0"]

        E_int = Eint * np.exp(1j * delta)
        I = (E_int - V) / (1j * Xdp)
        Pe = float((V * np.conj(I)).real)

        ddelta = w0 * omega
        domega = (Pm - Pe - D * omega) / (2.0 * H)

        # Cache for current_injection / algebraic_output / plotting.
        self._last_I = complex(I)
        self._last_Pe = Pe

        return {"delta": ddelta, "omega": domega}

    def current_injection(self, V: complex) -> complex:
        """Current injected by the machine into its bus, given terminal V.

        Computed fresh from the current state and the supplied V so the
        simulator can call this independently of derivatives().
        """
        E_int = self.params["Eint"] * np.exp(1j * self.state["delta"])
        return (E_int - V) / (1j * self.params["Xdp"])

    def algebraic_output(self) -> dict:
        """Quantities exposed for plotting / downstream models."""
        V = self.inputs["V_terminal"]
        I = self._last_I
        # Standard 5 (P, Q, |V|, Id, Iq) at the machine TERMINAL, in MACHINE dq.
        S = V * np.conj(I)
        # Rotate into machine dq frame using the standard PSSE/Kundur
        # convention (Park's transformation, d-axis 90° behind q-axis in
        # the direction of rotation).  Decompose X·exp(-j·delta) so that:
        #   X_q = Re(X·exp(-j·delta))     parallel to rotor q-axis (E')
        #   X_d = -Im(X·exp(-j·delta))    along rotor d-axis (field direction)
        # At steady state for an exporting lagging-PF generator:
        #   I_q > 0 (torque-producing)
        #   I_d > 0 (demagnetising)
        I_rot = I * np.exp(-1j * self.state["delta"])
        I_q = float(I_rot.real)
        I_d = -float(I_rot.imag)
        return {
            "P": float(S.real),
            "Q": float(S.imag),
            "|V|": float(abs(V)),
            "Id": I_d,
            "Iq": I_q,
            "delta": self.state["delta"],
            "omega": self.state["omega"],
            "Pe": self._last_Pe,
            "Pm": self.params["Pm"],
            "Eint": self.params["Eint"],
        }

    # ----- Norton equivalent (for solver coupling) --------------------

    def norton_admittance(self) -> complex:
        """Y = 1 / (j*X'd) — the shunt admittance to fold into Ybus."""
        return 1.0 / (1j * self.params["Xdp"])

    def norton_source(self) -> complex:
        """I_source = E' * Y, the current source in the Norton equivalent."""
        E_int = self.params["Eint"] * np.exp(1j * self.state["delta"])
        return E_int * self.norton_admittance()
