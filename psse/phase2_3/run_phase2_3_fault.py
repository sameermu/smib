"""Run the Phase 2.3 deep-inductive-fault scenario in PSSE.

Reproduces the smib Phase 2.3 scenario: 200 ms three-phase shunt
fault at the gen bus with Z_f = j*0.10 pu, on a 2-bus SMIB with a
GENROU machine (Kundur Table 4.2, D=3 load damping) driven by an
ST1A AVR (Ka=200, Tb=20 for TGR), a PSS1A stabiliser (Tw=5,
T1=T3=0.5, T2=T4=0.05, Ks=20) and a TGOV1 turbine governor
(R=0.05, T1=0.5, T2=2.5, T3=7.5, Vmax=1.0, Vmin=0.0).

Used to benchmark smib's Phase 2.3 full-classical-stack coupled
simulation against PSSE.  The headline check is a NEGATIVE result:
the governor must NOT change first-swing behaviour (peak angle,
CCT) relative to Phase 2.2 — its only visible action on this
scenario is the small Pm dip/wiggle opposing the rotor slip
oscillation (PMECH channel).

Usage from a PSSE-aware Python environment:
    python run_phase2_3_fault.py

Output:
    smib_phase2_3_fault.out   — channel file
    smib_phase2_3_fault.log   — text log

Smib reference numbers (h = 2 ms, GENROU + ST1A + PSS1A + TGOV1
coupled, D=3 load damping):

    Initialisation (identical PF to Phase 2.0/2.1/2.2):
        V1 from PF                : |V1| = 1.0178 pu, angle = +23.142 deg
        delta_0                   : 68.55 deg
        Eqp_0, Edp_0              : 0.9238, 0.4571 pu
        Efd_0 (= AVR output)      : 2.0256 pu
        Vref                      : 1.0279 pu
        Pref (= Pm_0 = Pe_0)      : 0.8000 pu
        x_v_0, x_r_0              : 0.8000, 0.8000 pu
        PSS states                : all 0 (dormant at SS)

    Deep inductive fault (Z_f = j0.10, t_clear = 200 ms):
        Peak rotor angle          : 102.5 deg (vs 103.3 gov-off —
                                    the valve barely moves in 200 ms)
        Peak Efd                  : 7.0 pu (hits Vrmax)
        Pm excursion              : a few % dip during the fault
                                    (droop opposes rotor acceleration),
                                    then small phase-lagged wiggles
                                    opposing the slip oscillation
        CCT                       : 333 ms (+6.4 ms vs Phase 2.2's
                                    326 ms — the valve's HP fraction
                                    eases Pm off during a near-CCT
                                    fault; helpful direction only)

    Primary-frequency-response cross-check (smib only — see README):
        Grid event Δf = -0.2 Hz (Δω̄ = -0.004 pu on 50 Hz)
        Final slip                : -0.00400 pu (locks to grid)
        ΔPm                       : +0.0796 pu at t = 40 s,
                                    → +0.0800 pu asymptotically
                                    (= -Δω̄/R exactly)
        Valve position            : 0.80 → 0.88 pu

PSSE-vs-smib expected agreement on the fault case (post 50 ms):
    Rotor-angle peak, period:  ~3-5 %
    PMECH dip shape:           qualitative match (wiggle opposing slip,
                               phase-lagged by the T1=0.5 s servo)
    CCT:                       ~5-10 ms

If your PSSE result differs by more, check (in order):
  1. DELT = 0.002 s
  2. MBASE = 100 MVA
  3. .dyr exactly matches smib_phase2_3.dyr
  4. TGOV1 VMAX/VMIN = 1.0/0.0 — if PMECH pins at an unexpected
     level mid-run, the valve limits are mis-entered
  5. TGOV1 R = 0.05 on machine base — PSSE droop is entered on
     MBASE; with MBASE = 100 MVA = system base no conversion is
     needed, but double-check if you changed MBASE.
"""
import os, sys

PSSE_PATH = r"C:\Program Files\PTI\PSSE35\35.6\PSSPY39"
if PSSE_PATH not in sys.path:
    sys.path.append(PSSE_PATH)

import psspy
import psse35

PSSE_OUT = "smib_phase2_3_fault.out"


def main():
    psspy.psseinit(50)

    if os.path.exists("smib_phase2_3.sav"):
        psspy.case("smib_phase2_3.sav")
    elif os.path.exists("../phase2_1/smib_phase2_1.sav"):
        psspy.case("../phase2_1/smib_phase2_1.sav")
    elif os.path.exists("../phase1/smib_phase1.sav"):
        psspy.case("../phase1/smib_phase1.sav")
    else:
        build_case()
    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])

    psspy.dyre_new([1, 1, 1, 1], "smib_phase2_3.dyr", "", "", "")

    psspy.machine_array_channel([1,  1, 1], 101, "1 ")   # ANGLE
    psspy.machine_array_channel([2,  2, 1], 101, "1 ")   # PELEC
    psspy.machine_array_channel([3,  6, 1], 101, "1 ")   # PMECH  ← the Phase 2.3 channel
    psspy.machine_array_channel([4,  4, 1], 101, "1 ")   # ETERM
    psspy.machine_array_channel([5,  5, 1], 101, "1 ")   # EFD
    psspy.voltage_channel(       [6, -1,  1, 101], "Vmag_101")
    psspy.machine_array_channel([7,  7, 1], 101, "1 ")   # SPEED

    psspy.dynamics_solution_param_2(realar3=0.002)
    psspy.strt_2([0, 0], PSSE_OUT)
    psspy.run(0, 1.0, 0, 1, 0)

    # Z_f = j0.10 pu shunt → B = -10 pu on 100 MVA.
    psspy.shunt_data(101, "1 ", 1, [0.0, -10.0])
    psspy.run(0, 1.20, 0, 1, 0)

    psspy.shunt_data(101, "1 ", 0, [0.0, 0.0])
    psspy.run(0, 10.0, 0, 1, 0)

    print(f"Done.  Channels written to {PSSE_OUT}")
    print("Compare PMECH against smib's res.traces['Pm'] — expect the")
    print("same few-percent dip during the fault and phase-lagged")
    print("wiggles opposing the slip oscillation afterwards.")


def build_case():
    """Same 2-bus case as Phase 1/2.x.  Identical PF result."""
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
    print("Phase 2.3 case built from scratch.")


if __name__ == "__main__":
    main()
