"""Run the Phase 1 deep-inductive-fault scenario in PSSE.

Reproduces the smib Phase 1 scenario: 100 ms three-phase shunt fault
at the gen bus with Z_f = j*0.10 pu, on a 2-bus SMIB with a GENCLS
machine.  Used to benchmark smib's Phase 1 notebook against PSSE.

Usage from a PSSE-aware Python environment:
    python run_phase1_fault.py

Output:
    smib_phase1_fault.out   — channel file (delta, omega, V, P, Q, ...)
    smib_phase1_fault.log   — text log

Comparison points (smib reference, h = 2 ms trapezoidal integrator):
    V1 from PF                : |V1| = 1.0178 pu, angle = +23.142 deg
    delta_0 (rotor angle)     : 35.495 deg
    |E'| (held constant)      : 1.1023 pu
    Mid-fault (t = 1.05 s)    : |V| = 0.348, P = +0.291, Q = +0.848
    Peak delta swing          : ~42 deg at t ~ 4.3 s (post-fault, undamped)
    Bolted-fault CCT (Z_f->0) : 191 ms (smib) vs 187 ms (analytic EAC)

PSSE results within ~3 % on the bolted-fault CCT and within ~2-5 %
on the post-fault period are a reasonable agreement (PSSE uses the
same trapezoidal integrator and the same swing equation; the only
expected difference is from numerical step size and how the shunt
fault is applied).
"""
import os, sys

# ---- locate PSSE (path varies by install; PSSE 35 default below) ----
PSSE_PATH = r"C:\Program Files\PTI\PSSE35\35.6\PSSPY39"
if PSSE_PATH not in sys.path:
    sys.path.append(PSSE_PATH)

import psspy
import psse35

PSSE_OUT = "smib_phase1_fault.out"


def main():
    psspy.psseinit(50)            # max bus count

    # ---- 1.  Load the power-flow case ----
    # Either save smib_phase1.sav from the GUI (see psse/phase1/README.md)
    # or build it from scratch with the function below.
    if os.path.exists("smib_phase1.sav"):
        psspy.case("smib_phase1.sav")
    else:
        build_case()
    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])      # full Newton solution

    # ---- 2.  Load dynamics ----
    psspy.dyre_new([1, 1, 1, 1], "smib_phase1.dyr", "", "", "")

    # ---- 3.  Set up channels for output ----
    # Channel ids matter for cross-referencing with smib's canonical 5.
    psspy.machine_array_channel([1,  1, 1], 101, "1 ")   # angle (delta)
    psspy.machine_array_channel([2,  2, 1], 101, "1 ")   # pu speed (omega)
    psspy.machine_array_channel([4,  4, 1], 101, "1 ")   # PELEC
    psspy.machine_array_channel([5,  5, 1], 101, "1 ")   # QELEC
    psspy.voltage_channel(       [6, -1,  1, 101], "Vmag_101")
    psspy.machine_array_channel([7,  7, 1], 101, "1 ")   # Efd

    # ---- 4.  Initialise dynamics, run pre-fault to t = 1.0 s ----
    psspy.dynamics_solution_param_2(realar3=0.001)        # dt = 1 ms
    psspy.strt_2([0, 0], PSSE_OUT)
    psspy.run(0, 1.0, 0, 1, 0)

    # ---- 5.  Apply 3-phase shunt fault Z_f = j*0.10 pu ----
    # PSSE shunt convention: G + jB in pu on system base.
    # Y_f = 1/(j*0.10) = -j*10 pu  →  G = 0, B = -10 pu (inductive).
    # NOTE: actual function name varies between psspy releases; in
    # PSSE 35 use shunt_data() to set the bus shunt directly.
    psspy.shunt_data(101, "1 ", 1, [0.0, -10.0])
    psspy.run(0, 1.10, 0, 1, 0)

    # ---- 6.  Clear the fault ----
    psspy.shunt_data(101, "1 ", 0, [0.0, 0.0])
    psspy.run(0, 5.0, 0, 1, 0)

    print(f"Done.  Channels written to {PSSE_OUT}")


def build_case():
    """Build the Phase 1 power-flow case from scratch.  Equivalent to
    the smib initial conditions: P=0.8, Q=0.2 at gen bus 101, slack
    at bus 102, line X = 0.5 pu, no resistance."""
    # System base
    psspy.base_frequency(50.0)

    # Bus data: bus number, name, base kV, code (1=PQ, 2=PV, 3=slack)
    psspy.bus_data_4([101, 1, 0, 0, 0], [   0.0,  18.0,  1.0, 0.0,  1.1, 0.9, 1.1, 0.9], "GEN BUS")
    psspy.bus_data_4([102, 1, 0, 0, 0], [   0.0,  18.0,  1.0, 0.0,  1.1, 0.9, 1.1, 0.9], "INF BUS")

    # Bus 102 = slack (code 3)
    psspy.bus_chng_4(102, 0, [3, 1, 1, 1], [_f]*8, "INF BUS")

    # Generator at 101 — PQ machine for GENCLS (no AVR), 100 MVA base.
    # P = 80 MW, Q = 20 MVAR, X'd in dynamics file, X_source for PF = X'd.
    psspy.plant_data_4(101, 0, [-1, 0], [1.0, 100.0])
    psspy.machine_data_4(101, "1 ",
        [1, 0, 0, 0, 0, 0, 0],
        [80.0, 20.0,  9999.0, -9999.0,  9999.0, -9999.0,
         100.0, 0.0, 0.30, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0])

    # Branch 101-102: R = 0, X = 0.5, B = 0
    psspy.branch_data_3(101, 102, "1 ",
        [1, 1, 1, 0, 0, 0, 0],
        [0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    print("Phase 1 case built from scratch.")


# Sentinel for "use default" in psspy array calls
_f = 0.0


if __name__ == "__main__":
    main()
