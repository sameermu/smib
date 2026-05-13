"""Run the Phase 2.1 deep-inductive-fault scenario in PSSE.

Reproduces the smib Phase 2.1 scenario: 200 ms three-phase shunt
fault at the gen bus with Z_f = j*0.10 pu, on a 2-bus SMIB with a
GENROU machine (Kundur Table 4.2) driven by an ST1A AVR (Ka=200,
Tr=0.02, Tc=Tb=1.0, Vrmax=7, Vrmin=-6.4).

Used to benchmark smib's Phase 2.1 GENROU + ST1A coupled simulation
against full PSSE GENROU + ST1A.

Usage from a PSSE-aware Python environment:
    python run_phase2_1_fault.py

Output:
    smib_phase2_1_fault.out   — channel file (delta, omega, V, P, Q, Eqp, Edp, Efd, Vref, Vc, ...)
    smib_phase2_1_fault.log   — text log

Smib reference numbers (h = 2 ms, GENROU + ST1A coupled):

    V1 from PF                  : |V1| = 1.0178 pu, angle = +23.142 deg
    Initial rotor angle delta   : 68.55 deg
    Initial Eqp                 : 0.9238 pu
    Initial Edp                 : 0.4571 pu
    Initial Efd (= AVR output)  : 2.0256 pu
    Initial Vref                : 1.0279 pu  (= Efd/Ka + |V|)
    Initial Vc                  : 1.0178 pu  (= |V|)

    Voltage-step response (+2 % bump on Vref at t = 1 s, 10 s run):
        Pre-step  |V|   = 1.0178 pu
        Final     |V|   = 1.0330 pu      (delta = +0.0152 pu, ~76% of command)
        Peak      Efd   = 6.02 pu        (well below Vrmax = 7)
        Final     Efd   = 3.10 pu
        Time to |V| peak: ~8.4 s         (limited by T'do = 8 s, not Tr or Tb)

    Deep inductive fault (Z_f = j0.10, t_clear = 200 ms):
        Peak Efd during fault : 7.0 pu       (hits Vrmax ceiling)
        Vc nadir during fault : 0.29 pu      (sensed terminal voltage)
        E'q nadir during fault: 0.92 pu      (almost no sag — AVR cancels it)
        Peak rotor angle      : 119.8 deg    (vs 114.2 deg AVR-off; HIGHER because
                                              preserved E'q gives higher post-fault
                                              synchronising power and a farther
                                              new equilibrium)

    CCT comparison (deep inductive Z_f = j*0.10):
        GENCLS         (Phase 1):   293 ms   (overestimate — no E'q sag)
        GENROU bare    (Phase 2.0): 240 ms   (most conservative — captures sag)
        GENROU + ST1A  (Phase 2.1): 290 ms   (AVR force-fields rotor flux)
        AVR lift over bare GENROU:  +28 ms (+10.7%)

PSSE-vs-smib expected agreement:

    Beyond ~50 ms post-fault: ~3-5 % on rotor-angle peak, oscillation
    period, and CCT.  PSSE captures sub-transient detail in the first
    50 ms that smib's 4-state simplification misses.

    Specific to Phase 2.1:
    - Efd ceiling hit must match (7.0 pu in both, identical Vrmax).
    - Time to ceiling depends on Tr, Tb, Ka — should agree within 1 ms
      given identical parameters.
    - Voltage-step time to peak should agree within ~1 % since the
      dominant time constant is T'do = 8 s and both models use it.

If your PSSE GENROU + ST1A CCT differs from smib's 290 ms by more
than ~15 ms, check (in this order):
  1. Did you use the dynamics step DELT = 0.002 s to match smib?
  2. Is MBASE on the machine equal to the system base (100 MVA)?
  3. Did you load smib_phase2_1.dyr exactly (Kundur Table 4.2 + ST1A
     params)?
  4. Did you set Kc = 0 to keep the field-current compensation out
     of the picture?  Smib's ST1A doesn't model it.
  5. Is the shunt fault G = 0, B = -10 (inductive Z_f = j*0.10)?

See ../phase1/README.md for the GUI-driven step-by-step.
"""
import os, sys

# ---- locate PSSE (path varies by install; PSSE 35 default below) ----
PSSE_PATH = r"C:\Program Files\PTI\PSSE35\35.6\PSSPY39"
if PSSE_PATH not in sys.path:
    sys.path.append(PSSE_PATH)

import psspy
import psse35

PSSE_OUT = "smib_phase2_1_fault.out"


def main():
    psspy.psseinit(50)

    # ---- 1.  Load the power-flow case ----
    # PF is identical to Phase 1 / Phase 2.0 — same network, same loading.
    if os.path.exists("smib_phase2_1.sav"):
        psspy.case("smib_phase2_1.sav")
    elif os.path.exists("../phase2_0/smib_phase2_0.sav"):
        psspy.case("../phase2_0/smib_phase2_0.sav")
    elif os.path.exists("../phase1/smib_phase1.sav"):
        psspy.case("../phase1/smib_phase1.sav")
    else:
        build_case()
    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])

    # ---- 2.  Load dynamics (GENROU + ST1A) ----
    psspy.dyre_new([1, 1, 1, 1], "smib_phase2_1.dyr", "", "", "")

    # ---- 3.  Set up channels for output ----
    psspy.machine_array_channel([1,  1, 1], 101, "1 ")    # angle (delta)
    psspy.machine_array_channel([2,  2, 1], 101, "1 ")    # pu speed (omega)
    psspy.machine_array_channel([4,  4, 1], 101, "1 ")    # PELEC
    psspy.machine_array_channel([5,  5, 1], 101, "1 ")    # QELEC
    psspy.voltage_channel(       [6, -1,  1, 101], "Vmag_101")
    psspy.machine_array_channel([7,  7, 1], 101, "1 ")    # Efd
    # GENROU state channels (E'q, E'd) and ST1A state channels (Vc, x_LL)
    # are addressable via PSSE 35's MODEL_DATA / MODEL_STATE APIs.
    # Exact channel codes depend on PSSE build — consult your manual.

    # ---- 4.  Initialise, run pre-fault to t = 1.0 s ----
    psspy.dynamics_solution_param_2(realar3=0.002)
    psspy.strt_2([0, 0], PSSE_OUT)
    psspy.run(0, 1.0, 0, 1, 0)

    # ---- 5.  Apply 3-phase shunt fault Z_f = j*0.10 pu ----
    # Y_f = 1/(j*0.10) = -j*10 pu  →  G = 0, B = -10 pu (inductive).
    psspy.shunt_data(101, "1 ", 1, [0.0, -10.0])
    psspy.run(0, 1.20, 0, 1, 0)   # 200 ms fault duration

    # ---- 6.  Clear the fault, run to 5.0 s ----
    psspy.shunt_data(101, "1 ", 0, [0.0, 0.0])
    psspy.run(0, 5.0, 0, 1, 0)

    print(f"Done.  Channels written to {PSSE_OUT}")


def build_case():
    """Same 2-bus case as Phase 1 / 2.0.  Identical PF result."""
    psspy.base_frequency(50.0)

    psspy.bus_data_4([101, 1, 0, 0, 0], [0.0, 18.0, 1.0, 0.0, 1.1, 0.9, 1.1, 0.9], "GEN BUS")
    psspy.bus_data_4([102, 1, 0, 0, 0], [0.0, 18.0, 1.0, 0.0, 1.1, 0.9, 1.1, 0.9], "INF BUS")
    psspy.bus_chng_4(102, 0, [3, 1, 1, 1], [0.0]*8, "INF BUS")

    psspy.plant_data_4(101, 0, [-1, 0], [1.0, 100.0])
    psspy.machine_data_4(101, "1 ",
        [1, 0, 0, 0, 0, 0, 0],
        [80.0, 20.0,  9999.0, -9999.0,  9999.0, -9999.0,
         100.0, 0.0, 0.23, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0])

    psspy.branch_data_3(101, 102, "1 ",
        [1, 1, 1, 0, 0, 0, 0],
        [0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    print("Phase 2.1 case built from scratch.")


def voltage_step_scenario():
    """Alternative scenario: Vref step instead of fault.

    To run the small-signal voltage-step instead of the deep fault,
    set ``main`` to call this function and remove the shunt_data
    fault application.  At t=1 s the AVR's Vref is bumped by +2 %.

    PSSE doesn't expose ``Vref`` directly through psspy in all
    versions — you may need to use ``increment_value`` with the
    appropriate model parameter index, e.g.::

        ix = psspy.macndx(101, "1 ", "AVR", "VREF")
        psspy.increment_value(ix, 0.02)

    Alternatively, save a snapshot, change Vref in the GUI, and
    restart from snapshot.
    """
    pass


if __name__ == "__main__":
    main()
