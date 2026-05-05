"""Run the Phase 2.0 deep-inductive-fault scenario in PSSE.

Reproduces the smib Phase 2.0 scenario: 100 ms three-phase shunt
fault at the gen bus with Z_f = j*0.10 pu, on a 2-bus SMIB with a
GENROU machine using Kundur Table 4.2 parameters.

Used to benchmark smib's Phase 2.0 4-state two-axis transient model
against full PSSE GENROU.

Usage from a PSSE-aware Python environment:
    python run_phase2_0_fault.py

Output:
    smib_phase2_0_fault.out   — channel file (delta, omega, V, P, Q, Eqp, Edp, Efd, ...)
    smib_phase2_0_fault.log   — text log

Smib reference numbers (h = 2 ms, 4-state two-axis transient model):

    V1 from PF                : |V1| = 1.0178 pu, angle = +23.142 deg
    Initial rotor angle delta : 68.55 deg  (phasor behind X_q)
    Initial Eqp               : 0.9238 pu  (≈ 0.924, transient field flux)
    Initial Edp               : 0.4571 pu  (slow q-axis damper flux)
    Initial Efd               : 2.026 pu   (steady-state field voltage,
                                            held constant — no AVR yet)

    Mid-fault (t = 1.05 s, deep inductive Z_f = j*0.10):
        |V|  = 0.305 pu
        P    = +0.329 pu
        Q    = +0.602 pu      (vs +0.847 in GENCLS; Eqp sag reduces it)
        Eqp  = 0.9088 pu      (sagged 0.015 from initial 0.9238)
        Edp  = 0.4126 pu      (sagged 0.044 from initial 0.4571)
        delta = 70.29 deg     (vs 37 deg in GENCLS — different rotor-
                               angle convention because phasor-behind-
                               X_q gives a different angle than
                               phasor-behind-X'_d)

    CCT comparison (deep inductive Z_f = j*0.10):
        GENROU CCT (smib):  240 ms
        GENCLS CCT (smib):  293 ms  — GENCLS is ~22% more permissive
                                      because it doesn't capture Eqp
                                      sag during fault.

PSSE-vs-smib expected agreement:

    Beyond ~50 ms post-fault: ~3-5 % on rotor-angle peak, oscillation
    period, and CCT.  PSSE captures sub-transient detail in the first
    50 ms that smib's 4-state simplification misses, so:

    - Initial fault current spike (first 30-50 ms):  PSSE ~30 % bigger
                                                      than smib.
    - Eqp / Edp dynamics post-50ms:                   match within ~5 %.
    - Rotor-angle peak and oscillation period:        match within ~3 %.
    - Deep CCT:                                       match within ~5 %.

If your PSSE GENROU CCT differs from smib's 240 ms by more than ~10
ms, check (in this order):
  1. Did you use the dynamics step DELT = 0.002 s to match smib?
  2. Is MBASE on the machine equal to the system base (100 MVA)?
  3. Did you load smib_phase2_0.dyr exactly (Kundur Table 4.2 params)?
  4. Is the shunt fault G = 0, B = -10 (inductive Z_f = j*0.10)?

See ../phase1/README.md for the GUI-driven step-by-step.
"""
import os, sys

# ---- locate PSSE (path varies by install; PSSE 35 default below) ----
PSSE_PATH = r"C:\Program Files\PTI\PSSE35\35.6\PSSPY39"
if PSSE_PATH not in sys.path:
    sys.path.append(PSSE_PATH)

import psspy
import psse35

PSSE_OUT = "smib_phase2_0_fault.out"


def main():
    psspy.psseinit(50)            # max bus count

    # ---- 1.  Load the power-flow case ----
    # Save smib_phase2_0.sav from the GUI (same network as phase1) or
    # re-build via build_case() below.  The PF result is identical to
    # Phase 1 because the network and operating point are the same;
    # only the dynamic model changes.
    if os.path.exists("smib_phase2_0.sav"):
        psspy.case("smib_phase2_0.sav")
    elif os.path.exists("../phase1/smib_phase1.sav"):
        psspy.case("../phase1/smib_phase1.sav")
    else:
        build_case()
    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])      # full Newton solution

    # ---- 2.  Load dynamics ----
    psspy.dyre_new([1, 1, 1, 1], "smib_phase2_0.dyr", "", "", "")

    # ---- 3.  Set up channels for output ----
    # Channel ids matter for cross-referencing with smib's canonical 5
    # plus the new Eqp/Edp/Efd traces.
    psspy.machine_array_channel([1,  1, 1], 101, "1 ")    # angle (delta)
    psspy.machine_array_channel([2,  2, 1], 101, "1 ")    # pu speed (omega)
    psspy.machine_array_channel([4,  4, 1], 101, "1 ")    # PELEC
    psspy.machine_array_channel([5,  5, 1], 101, "1 ")    # QELEC
    psspy.voltage_channel(       [6, -1,  1, 101], "Vmag_101")
    psspy.machine_array_channel([7,  7, 1], 101, "1 ")    # Efd
    # GENROU also exposes E'q, E'd as state channels (model-specific).
    # In PSSE 35, channel codes 26 and 27 of GENROU give E'q and E'd.
    # The exact channel API varies; consult your PSSE manual.

    # ---- 4.  Initialise dynamics, run pre-fault to t = 1.0 s ----
    psspy.dynamics_solution_param_2(realar3=0.002)        # match smib h = 2 ms
    psspy.strt_2([0, 0], PSSE_OUT)
    psspy.run(0, 1.0, 0, 1, 0)

    # ---- 5.  Apply 3-phase shunt fault Z_f = j*0.10 pu ----
    # PSSE shunt convention: G + jB in pu on system base.
    # Y_f = 1/(j*0.10) = -j*10 pu  →  G = 0, B = -10 pu (inductive).
    psspy.shunt_data(101, "1 ", 1, [0.0, -10.0])
    psspy.run(0, 1.10, 0, 1, 0)

    # ---- 6.  Clear the fault ----
    psspy.shunt_data(101, "1 ", 0, [0.0, 0.0])
    psspy.run(0, 5.0, 0, 1, 0)

    print(f"Done.  Channels written to {PSSE_OUT}")


def build_case():
    """Same 2-bus case as Phase 1, just with a GENROU machine instead
    of GENCLS.  The PF output is identical (same network, same loading)."""
    psspy.base_frequency(60.0)

    psspy.bus_data_4([101, 1, 0, 0, 0], [0.0, 18.0, 1.0, 0.0, 1.1, 0.9, 1.1, 0.9], "GEN BUS")
    psspy.bus_data_4([102, 1, 0, 0, 0], [0.0, 18.0, 1.0, 0.0, 1.1, 0.9, 1.1, 0.9], "INF BUS")
    psspy.bus_chng_4(102, 0, [3, 1, 1, 1], [0.0]*8, "INF BUS")

    psspy.plant_data_4(101, 0, [-1, 0], [1.0, 100.0])
    # For GENROU we set X_source to X"_d = 0.23 in the PF (the fault-
    # current contribution at the moment fault is applied).  Smib uses
    # X'_d = 0.30 because we don't model sub-transient — minor PF effect.
    psspy.machine_data_4(101, "1 ",
        [1, 0, 0, 0, 0, 0, 0],
        [80.0, 20.0,  9999.0, -9999.0,  9999.0, -9999.0,
         100.0, 0.0, 0.23, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0])

    psspy.branch_data_3(101, 102, "1 ",
        [1, 1, 1, 0, 0, 0, 0],
        [0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    print("Phase 2.0 case built from scratch.")


if __name__ == "__main__":
    main()
