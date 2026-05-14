"""ST1A — IEEE 421.5 Type ST1A static exciter.

Reference: IEEE Std 421.5-2016 §5.1.

Topology (simplified for smib pedagogy — drops the OEL/UEL inputs and
the rate feedback transformer KF/TF, which are off in our default
case):

    |V_term| ---> [ 1/(1+sTr) ] -----> Vc (sensed terminal voltage)
                                          |
                                          v  -
    Vref + Vpss --(+)----------------> error ---> [Ka·(1+sTc)/(1+sTb)] ---> [hard limit Vrmax/Vrmin] ---> Efd
                  +
                  ^
                  |
                  Vpss (0 here; non-zero in Phase 2.2)

States (3 total)
----------------
    Vc      [pu]   sensed terminal voltage (low-pass output of |V|)
    x_LL    [pu]   internal state of the (Tc, Tb) lead-lag
    x_int   [pu]   integrator state of the regulator output (only present
                   if we treat the regulator as a true PI; in pure ST1A
                   the regulator is a lag-lead WITHOUT an integrator and
                   the static gain Ka acts directly on the error.  We
                   implement the simpler one-state form here: regulator
                   output = Ka · (lead-lag of error), no integrator.
                   This matches PSSE's default ST1A when Ki = 0.)

We end up with **2 states** in this simplified form:

    Vc      — sensed terminal voltage, time constant Tr ~ 0.01-0.05 s
    x_LL    — lead-lag internal state, time constants Tc (zero) and Tb (pole).
              With the default Tb=10, Tc=1 this acts as transient gain
              reduction: full Ka at DC, Ka·Tc/Tb at frequencies above
              the pole.  See `__init__` docstring for the rationale.

Differential equations
----------------------
::

    Tr · dVc/dt = |V_term| - Vc

For the lead-lag block with input u = error and output y, the
canonical state-space realisation with one internal state x is:

    Tb · dx/dt = u - x
    y = (Tc/Tb) · u + (1 - Tc/Tb) · x

We then scale by Ka and clamp to get Efd.

Inputs
------
    V_terminal_mag  [pu]   magnitude of the generator terminal voltage
    Vref            [pu]   AVR setpoint (set at init)
    Vpss            [pu]   PSS contribution (0 here, non-zero in Phase 2.2)

Outputs (algebraic_output)
--------------------------
    Vc, error, Efd

Initialisation
--------------
At t=0 the field voltage Efd_init is given (it's whatever GENROU's
init computed to satisfy steady-state Eqp).  Working backwards through
the regulator:

    Efd_init / Ka = lead-lag steady-state output = error_init  (since
        at SS the lead-lag passes the input through with gain 1)
    error_init = Vref - Vc + Vpss
    Vc_init = |V_term_init|  (the sensor has zero error at SS)

So Vref = error_init + Vc_init - Vpss  (with Vpss=0)
       = Efd_init/Ka + |V_term_init|.

x_LL_init: with the lead-lag at steady state with input u_ss = error_init,
the internal state x = u_ss = error_init (because Tb·0 = u - x ⇒ x = u).
"""
from __future__ import annotations

import math

from .base import Model


class ST1A(Model):
    """Type ST1A static exciter, simplified (no OEL/UEL/rate feedback).

    Two states: Vc (sensed terminal voltage), x_LL (lead-lag internal).
    """

    name = "ST1A"
    state_keys = ("Vc", "x_LL")

    def __init__(self,
                 name: str = "ST1A",
                 Tr: float = 0.02,        # voltage transducer time constant [s]
                 Ka: float = 200.0,       # regulator gain (steady-state)
                 Ta: float = 0.0,         # placeholder (not used in this simplified form)
                 Tb: float = 20.0,        # lead-lag pole [s]  — transient gain reduction
                 Tc: float = 1.0,         # lead-lag zero [s]
                 Vrmax: float = 7.0,      # field voltage upper limit (ceiling)
                 Vrmin: float = -6.4):    # field voltage lower limit (negative ceiling)
        """Default parameters — typical thermal-unit ST1A values from
        IEEE 421.5 examples with **transient gain reduction (TGR)**.

        Tb > Tc gives the lead-lag a low-pass character: DC gain is 1
        (full Ka for steady-state voltage regulation) but at the rotor
        swing frequency (~0.7 Hz on this SMIB) the effective gain is
        Ka·Tc/Tb / |1 + jωTb| ≈ Ka·Tc/Tb = 20 instead of 200.  This
        is the textbook way to keep high steady-state accuracy without
        destabilising the electromechanical mode — exactly the role a
        Power System Stabiliser fills more aggressively in Phase 2.2.

        Set Tc = Tb (pure-gain lead-lag) to revert to the un-detuned
        form; it's pedagogically clean but produces poorly-damped
        post-fault rotor oscillations.
        """
        params = {
            "Tr": Tr, "Ka": Ka, "Ta": Ta, "Tb": Tb, "Tc": Tc,
            "Vrmax": Vrmax, "Vrmin": Vrmin,
            "Vref": 1.0,    # filled at init
        }
        super().__init__(name, params)
        self.inputs = {
            "V_terminal_mag": 1.0,   # |V_t|
            "Vpss": 0.0,             # PSS contribution (Phase 2.2)
        }
        self._last_error: float = 0.0
        self._last_Efd: float = 0.0

    # ----- helpers -----------------------------------------------------

    def _regulator_output(self, error: float) -> float:
        """Lead-lag block followed by gain Ka and limiter.  Returns the
        unlimited Efd (we apply the clamp outside)."""
        Tb, Tc = self.params["Tb"], self.params["Tc"]
        x = self.state["x_LL"]
        # y = (Tc/Tb) · error + (1 - Tc/Tb) · x_LL
        y = (Tc / Tb) * error + (1.0 - Tc / Tb) * x
        return self.params["Ka"] * y

    def _clamp(self, Efd: float) -> float:
        return max(self.params["Vrmin"], min(self.params["Vrmax"], Efd))

    # ----- Model interface --------------------------------------------

    def initialise(self, V: complex, S: complex = 0j, Efd_init: float = 0.0,
                   **kwargs) -> None:
        """Back-calculate Vref, Vc, and x_LL so that:
        - the AVR's steady-state Efd output equals Efd_init,
        - dVc/dt = 0 (so Vc = |V_t|),
        - dx_LL/dt = 0 (so x_LL = error at SS).

        ``Efd_init`` is the field voltage demanded by the generator at
        steady state (e.g. GENROU's params['Efd'] after its initialise()).
        """
        Vmag = abs(V)
        Vc_init = Vmag
        # error such that Ka · (lead-lag-passthrough(error)) = Efd
        # At SS the lead-lag passes the input through unchanged, so:
        error_init = Efd_init / self.params["Ka"]
        x_LL_init = error_init
        Vref = error_init + Vc_init - self.inputs.get("Vpss", 0.0)

        self.state["Vc"] = float(Vc_init)
        self.state["x_LL"] = float(x_LL_init)
        self.params["Vref"] = float(Vref)
        self.inputs["V_terminal_mag"] = float(Vmag)
        self._last_error = float(error_init)
        self._last_Efd = float(Efd_init)

    def derivatives(self) -> dict:
        Vmag = float(self.inputs["V_terminal_mag"])
        Vpss = float(self.inputs.get("Vpss", 0.0))
        Vref = self.params["Vref"]
        Tr, Tb = self.params["Tr"], self.params["Tb"]

        Vc = self.state["Vc"]
        x_LL = self.state["x_LL"]

        # 1. Voltage transducer
        dVc_dt = (Vmag - Vc) / Tr

        # 2. Error and lead-lag
        error = Vref + Vpss - Vc
        dx_LL_dt = (error - x_LL) / Tb

        # Cache for output.
        self._last_error = float(error)
        Efd_unlim = self._regulator_output(error)
        self._last_Efd = float(self._clamp(Efd_unlim))

        return {"Vc": dVc_dt, "x_LL": dx_LL_dt}

    def current_injection(self, V: complex) -> complex:
        """The AVR is not connected to the network; no current injection."""
        return 0j

    def algebraic_output(self) -> dict:
        return {
            "Vc": float(self.state["Vc"]),
            "Vref": float(self.params["Vref"]),
            "error": float(self._last_error),
            "Efd": float(self._last_Efd),
        }
