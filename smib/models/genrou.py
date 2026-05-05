"""GENROU — round-rotor synchronous machine, **4-state two-axis transient
model** (no sub-transient detail).

This is a deliberately simplified Phase 2.0 implementation that captures
the full TRANSIENT-time-scale flux dynamics (the regime ~50 ms to a few
seconds after a disturbance) but skips the SUB-TRANSIENT regime
(~tens of ms) handled by full PSSE GENROU.

For the *physics* — what each parameter means, why X"_d < X'_d < X_d,
the equivalent circuits, the open-circuit time constants, the
saturation function — see ``docs/genrou_physical_foundations.md`` in
the repo root.

States (4 total)
----------------
Mechanical (same two as GENCLS):
    delta   [rad]    rotor angle, global synchronous DQ frame
    omega   [pu]     slip (omega_rotor - omega_0) / omega_0

Electrical (NEW):
    Eqp     [pu]     d-axis transient EMF, ∝ field winding flux linkage
    Edp     [pu]     q-axis transient EMF, ∝ slow q-axis damper flux

Differential equations
----------------------
::

    T'_d0 · dEqp/dt = Efd - Eqp - (Xd - X'd) · Id - Sat(|Eqp|) · Eqp
    T'_q0 · dEdp/dt = -Edp + (Xq - X'q) · Iq

    ddelta/dt       = omega0 · omega
    2H · domega/dt  = Pm - Pe - D · omega

with::

    Pe = Vd · Id + Vq · Iq

Output (algebraic) equations — terminal voltage in machine dq from
the transient EMF states and the stator currents (R_a = 0)::

    Vq = -X'd · Id + Eqp
    Vd = +X'q · Iq + Edp

These are derived from Park's stator equations at synchronous speed
(ω = 1 pu, dψ/dt = 0 within an integration step) using the transient
flux linkages ψ_d = -X'd · Id + Eqp and ψ_q = -X'q · Iq - Edp.

Saturation
----------
PSSE quadratic ``S_E(x) = B (x - A)^2`` for ``x > A``, fit from
``S_E(1.0) = S1`` and ``S_E(1.2) = S2``.  Acts on |Eqp| (effectively
the air-gap flux on the d-axis, since at steady state |Eqp| is the
unsaturated open-circuit voltage that the field flux linkage produces
on the stator).

dq sign convention
------------------
PSSE / Kundur, **Park's transformation with d-axis 90° behind q-axis**
in the direction of rotation::

    X_q = Re(X_global · exp(-j · delta))
    X_d = -Im(X_global · exp(-j · delta))

After this rotation, the synchronous EMF E_q lies on the positive real
axis (the q-axis) by definition of delta from the "phasor behind X_q"
construction.  At steady state for an exporting lagging-PF generator:

    Iq > 0  (torque-producing component, parallel to E_q)
    Id > 0  (demagnetising component, perpendicular to E_q)
    Eqp > 0 (dominant — close to |E_a|)
    Edp varies with saliency; small for round rotor, larger for salient

This convention matches GENCLS in the same repo and the standard
PSSE/Kundur references.

Differences vs full PSSE GENROU
-------------------------------
**Missing:** sub-transient detail.  Full PSSE GENROU adds two more
states (psidpp = ψ"_d, psiqpp = ψ"_q) and uses sub-transient reactances
X"_d, X"_q in the output equations instead of X'_d, X'_q.  This means:

- The first ~50 ms of any fault response will look slightly different
  in our 4-state model than in PSSE.  In PSSE, the machine impedance
  during the first 50 ms is X"_d (≈0.23 in Kundur params), so the
  fault current spike is ~30% larger than in our 4-state model
  (which uses X'_d ≈ 0.30).
- Beyond ~50 ms (i.e. for the rest of the swing) both models agree
  to the level of the constant-flux-linkage approximation.  The
  rotor-angle peak, oscillation period, and CCT should all match
  PSSE within a few percent.
- Saturation in our model acts only on the d-axis transient EMF.
  Full GENROU saturates the air-gap flux which couples both d- and
  q-axis sub-transient dynamics.
- D and Q axis damping (the "amortisseur" in physical hardware) is
  represented only by the slow time constant T'_q0 here.  Full GENROU
  has the fast T"_q0 ≈ 0.07 s in addition.

**Net effect on the canonical-5 traces** for a typical 100 ms fault:

- |V|, P, Q, delta, omega, Iq:  match PSSE within ~5 % everywhere except
  the first 50 ms, where X' vs X" gives a ~30 % discrepancy on the
  initial fault-current spike.
- Id (the d-axis "fault current" component): same ~30 % offset for
  the first 50 ms, then matches.

**When to upgrade to full GENROU (Phase 2.0b):** if you want PSSE-
accurate first-cycle fault current magnitudes (relevant for protection
coordination studies, generator differential protection, etc.).  For
the swing-stability and AVR/PSS pedagogy in Phase 2.1+, the 4-state
model is sufficient.
"""
from __future__ import annotations

import math

import numpy as np

from .base import Model


# ---------- saturation helpers -----------------------------------------

def _sat_AB(S1: float, S2: float) -> tuple[float, float]:
    """Fit the PSSE quadratic saturation parameters A, B so that
    S_E(1.0) = S1 and S_E(1.2) = S2.

    S_E(x) = B * (x - A)^2  for x > A.

    For S1 = S2 = 0 (saturation off), returns A = 0, B = 0.
    """
    if S1 <= 0 and S2 <= 0:
        return 0.0, 0.0
    if S1 <= 0 or S2 <= 0:
        raise ValueError("S(1.0) and S(1.2) must both be positive (or both zero).")
    sqrt_ratio = math.sqrt(S2 / S1)
    A = (1.2 - 1.0 * sqrt_ratio) / (1.0 - sqrt_ratio)
    B = S1 / (1.0 - A) ** 2
    return A, B


def _sat(x: float, A: float, B: float) -> float:
    """Quadratic saturation function S_E(x).  Operates on |x|."""
    ax = abs(x)
    if ax <= A or B == 0.0:
        return 0.0
    return B * (ax - A) ** 2


# ---------- model ------------------------------------------------------

class GENROU(Model):
    """Round-rotor synchronous machine, 4-state two-axis transient model.

    See module docstring for the full equation set, conventions, and
    the list of differences vs full PSSE GENROU.
    """

    name = "GENROU"
    state_keys = ("delta", "omega", "Eqp", "Edp")

    def __init__(self,
                 name: str = "GENROU",
                 # mechanical
                 H: float = 4.0,
                 D: float = 0.0,
                 f0: float = 60.0,
                 # d-axis (full GENROU has Xdpp; we ignore it in 4-state)
                 Xd: float = 1.81, Xdp: float = 0.30,
                 Tdo_p: float = 8.0,
                 # q-axis (full GENROU has Xqpp; we ignore it in 4-state)
                 Xq: float = 1.76, Xqp: float = 0.65,
                 Tqo_p: float = 1.0,
                 # saturation (acts on |Eqp|)
                 S1: float = 0.13, S2: float = 0.50):
        """Default values are Kundur Table 4.2 (thermal unit, round rotor).

        Reactances all in pu on machine MVA base, time constants in
        seconds.

        Parameters Xdpp, Xqpp, Xl, Tdo_pp, Tqo_pp from full PSSE GENROU
        are absent because we don't model sub-transient dynamics here —
        see module docstring for what this approximation costs.
        """
        # Sanity-check the reactance hierarchy.  Allow equality so the
        # GENCLS-equivalence-limit case (X_d = X'_d, T'_d0 -> infty) is
        # constructible for regression tests.
        assert 0 < Xdp <= Xd, f"Bad d-axis order: X'={Xdp}, X={Xd}"
        assert 0 < Xqp <= Xq, f"Bad q-axis order: X'={Xqp}, X={Xq}"

        Asat, Bsat = _sat_AB(S1, S2)

        params = {
            "H": H, "D": D, "omega0": 2 * math.pi * f0,
            "Xd": Xd, "Xdp": Xdp,
            "Tdo_p": Tdo_p,
            "Xq": Xq, "Xqp": Xqp,
            "Tqo_p": Tqo_p,
            "S1": S1, "S2": S2, "Asat": Asat, "Bsat": Bsat,
            # Init-time constants (filled by initialise()):
            "Efd": 0.0, "Pm": 0.0,
        }
        super().__init__(name, params)
        self.inputs = {"V_terminal": complex(1.0, 0.0)}

        # cached for plotting/diagnostics
        self._last_I = 0j
        self._last_Pe = 0.0
        self._last_Vd = 0.0
        self._last_Vq = 0.0
        self._last_Id = 0.0
        self._last_Iq = 0.0

    # ----- helpers -----------------------------------------------------

    def _to_machine_dq(self, X_global: complex) -> tuple[float, float]:
        """Project a global-DQ phasor into machine dq.

        Standard PSSE/Kundur convention:
            X_q = Re(X · exp(-j·delta))
            X_d = -Im(X · exp(-j·delta))
        """
        rot = X_global * np.exp(-1j * self.state["delta"])
        Xq = float(rot.real)
        Xd = -float(rot.imag)
        return Xd, Xq

    def _stator_currents(self, V_terminal: complex) -> tuple[float, float, float, float]:
        """Solve algebraically for (Id, Iq) given V_terminal and the
        current internal Eqp, Edp.  R_a = 0:

            Vq = -X'd · Id + Eqp     ⇒  Id = (Eqp - Vq) / X'd
            Vd =  X'q · Iq + Edp     ⇒  Iq = (Vd - Edp) / X'q
        """
        Xdp, Xqp = self.params["Xdp"], self.params["Xqp"]
        Eqp, Edp = self.state["Eqp"], self.state["Edp"]

        Vd, Vq = self._to_machine_dq(V_terminal)

        Id = (Eqp - Vq) / Xdp
        Iq = (Vd - Edp) / Xqp
        return Id, Iq, Vd, Vq

    # ----- Model interface --------------------------------------------

    def initialise(self, V: complex, S: complex, **kwargs) -> None:
        """Back-calculate every state and the constants Efd, Pm so that
        derivatives() returns zero at t=0.

        Five-step protocol:

        1. Terminal current from PF: I = conj(S/V).
        2. Rotor angle from "phasor behind X_q":
             E_a_apparent = V + j · X_q · I
             delta = angle(E_a_apparent).
           This is the canonical Kundur classical result for round-
           rotor with full saliency — it places the synchronous EMF on
           the q-axis in machine dq.
        3. Project V and I into machine dq (standard PSSE convention).
        4. Compute Eqp and Edp from the steady-state output equations:
             Eqp = Vq + X'd · Id
             Edp = Vd - X'q · Iq
           Verify that the differential equations also evaluate to
           zero (Edp_diff_ss = (Xq - X'q) · Iq must equal Edp from
           output).  Mismatch indicates a sign/convention error.
        5. Compute Efd and Pm:
             Efd = Eqp + (Xd - X'd) · Id + S_E(|Eqp|) · Eqp
             Pm  = Vd · Id + Vq · Iq
        """
        Xd, Xdp = self.params["Xd"], self.params["Xdp"]
        Xq, Xqp = self.params["Xq"], self.params["Xqp"]
        Asat, Bsat = self.params["Asat"], self.params["Bsat"]

        # ----- step 1 — terminal current
        I_term = np.conj(S / V)

        # ----- step 2 — rotor angle from phasor behind X_q
        E_a_apparent = V + 1j * Xq * I_term
        delta = float(np.angle(E_a_apparent))

        # ----- step 3 — project V and I into machine dq (with delta now set)
        self.state["delta"] = delta
        Vd, Vq = self._to_machine_dq(V)
        Id_rotated = -float((I_term * np.exp(-1j * delta)).imag)
        Iq_rotated = +float((I_term * np.exp(-1j * delta)).real)
        Id, Iq = Id_rotated, Iq_rotated

        # ----- step 4 — Eqp and Edp from output eqs (steady-state)
        Eqp = Vq + Xdp * Id
        Edp = Vd - Xqp * Iq

        # Cross-check: the differential equations should give the same
        # SS values.  These two consistency checks are what break when
        # sign conventions are inconsistent.
        Edp_from_diff_ss = (Xq - Xqp) * Iq
        if abs(Edp - Edp_from_diff_ss) > 1e-3:
            raise RuntimeError(
                f"GENROU init inconsistency on q-axis: "
                f"Edp from output = {Edp:.4f}, from diff = {Edp_from_diff_ss:.4f}. "
                "Check dq sign conventions."
            )

        # ----- step 5 — Efd and Pm
        Se = _sat(Eqp, Asat, Bsat)
        Efd = Eqp + (Xd - Xdp) * Id + Se * Eqp
        Pe = Vd * Id + Vq * Iq

        # Write everything into params/state.
        self.state["omega"] = 0.0
        self.state["Eqp"] = float(Eqp)
        self.state["Edp"] = float(Edp)
        self.params["Efd"] = float(Efd)
        self.params["Pm"] = float(Pe)
        self.inputs["V_terminal"] = V

        # Cache for diagnostics.
        self._last_I = complex(I_term)
        self._last_Pe = float(Pe)
        self._last_Vd, self._last_Vq = Vd, Vq
        self._last_Id, self._last_Iq = Id, Iq

    def derivatives(self) -> dict:
        V = self.inputs["V_terminal"]
        Xd, Xdp = self.params["Xd"], self.params["Xdp"]
        Xq, Xqp = self.params["Xq"], self.params["Xqp"]
        Tdo_p = self.params["Tdo_p"]
        Tqo_p = self.params["Tqo_p"]
        Asat, Bsat = self.params["Asat"], self.params["Bsat"]
        Efd, Pm = self.params["Efd"], self.params["Pm"]
        H, D = self.params["H"], self.params["D"]
        w0 = self.params["omega0"]

        omega = self.state["omega"]
        Eqp = self.state["Eqp"]
        Edp = self.state["Edp"]

        # Algebraic stator currents.
        Id, Iq, Vd, Vq = self._stator_currents(V)

        # Saturation on |Eqp|.
        Se = _sat(Eqp, Asat, Bsat)

        # Two electrical equations (transient time scale).
        dEqp_dt = (Efd - Eqp - (Xd - Xdp) * Id - Se * Eqp) / Tdo_p
        dEdp_dt = (-Edp + (Xq - Xqp) * Iq) / Tqo_p

        # Mechanical (same as GENCLS).
        Pe = Vd * Id + Vq * Iq
        ddelta_dt = w0 * omega
        domega_dt = (Pm - Pe - D * omega) / (2.0 * H)

        # Cache for plotting.
        self._last_Pe = float(Pe)
        self._last_Vd, self._last_Vq = Vd, Vq
        self._last_Id, self._last_Iq = Id, Iq

        return {
            "delta": ddelta_dt,
            "omega": domega_dt,
            "Eqp": dEqp_dt,
            "Edp": dEdp_dt,
        }

    def current_injection(self, V: complex) -> complex:
        """Current injected by the machine into its bus given terminal V."""
        Id, Iq, _, _ = self._stator_currents(V)
        # Inverse of standard rotation: machine dq → global DQ.
        # We have I_q = Re(I·exp(-j·δ)), I_d = -Im(I·exp(-j·δ)).
        # ⇒ I·exp(-j·δ) = I_q - j·I_d
        # ⇒ I = (I_q - j·I_d) · exp(j·δ)
        I_global = complex(Iq, -Id) * np.exp(1j * self.state["delta"])
        return I_global

    def algebraic_output(self) -> dict:
        """Quantities exposed for plotting and downstream models."""
        V = self.inputs["V_terminal"]
        I = self.current_injection(V)
        S = V * np.conj(I)
        return {
            "P": float(S.real),
            "Q": float(S.imag),
            "|V|": float(abs(V)),
            "Id": float(self._last_Id),
            "Iq": float(self._last_Iq),
            "delta": float(self.state["delta"]),
            "omega": float(self.state["omega"]),
            "Eqp": float(self.state["Eqp"]),
            "Edp": float(self.state["Edp"]),
            "Pe": float(self._last_Pe),
            "Pm": float(self.params["Pm"]),
            "Efd": float(self.params["Efd"]),
        }

    # ----- Norton equivalent for solver coupling ----------------------

    def norton_admittance(self) -> complex:
        """Y_norton = 1 / (j · X'_avg) where X'_avg is the average
        transient reactance (we use the mean of X'_d and X'_q to
        accommodate saliency).  For round rotor X'_d ≈ X'_q so this is
        a tight approximation; for salient pole consider the full
        2x2 admittance.
        """
        Xp_avg = 0.5 * (self.params["Xdp"] + self.params["Xqp"])
        return 1.0 / (1j * Xp_avg)

    def norton_source(self) -> complex:
        """Norton equivalent: I_source = E_p_global · Y_norton, where
        E'_global = (Eqp + j·Edp) rotated back to global DQ.

        With our convention I·exp(-j·δ) = I_q - j·I_d, the inverse
        rotation is I_global = (I_q - j·I_d) · exp(j·δ).  Same applies
        to the EMF phasor.
        """
        E_p_dq_complex = complex(self.state["Eqp"], -self.state["Edp"])
        E_p_global = E_p_dq_complex * np.exp(1j * self.state["delta"])
        return E_p_global * self.norton_admittance()
