"""Run the Phase 2.2 deep-inductive-fault scenario in PSSE.

Reproduces the smib Phase 2.2 scenario: 200 ms three-phase shunt
fault at the gen bus with Z_f = j*0.10 pu, on a 2-bus SMIB with a
GENROU machine (Kundur Table 4.2, D=3 load damping) driven by an
ST1A AVR (Ka=200, Tb=20 for TGR) with a PSS1A stabiliser
(Tw=5, T1=T3=0.5, T2=T4=0.05, Ks=20, Vstmax=±0.10).

Used to benchmark smib's Phase 2.2 GENROU + ST1A + PSS1A coupled
simulation against full PSSE GENROU + ST1A + PSS1A.

Usage from a PSSE-aware Python environment:
    python run_phase2_2_fault.py

Output:
    smib_phase2_2_fault.out   — channel file
    smib_phase2_2_fault.log   — text log

Smib reference numbers (h = 2 ms, GENROU + ST1A + PSS1A coupled,
D=3 load damping):

    Initialisation (identical PF to Phase 2.0 / 2.1):
        V1 from PF                : |V1| = 1.0178 pu, angle = +23.142 deg
        delta_0                   : 68.55 deg
        Eqp_0, Edp_0              : 0.9238, 0.4571 pu
        Efd_0 (= AVR output)      : 2.0256 pu
        Vref                      : 1.0279 pu
        Vc                        : 1.0178 pu
        x_w, x_LL1, x_LL2         : all 0.0 (PSS dormant at SS)
        Vpss                      : 0.0 pu

    Deep inductive fault (Z_f = j0.10, t_clear = 200 ms):
        Peak Efd                  : 7.0 pu (hits Vrmax)
        Peak |Vpss|               : 0.10 pu (hits ±Vstmax)
        E'q nadir                 : 0.91 pu (slightly worse than AVR-only
                                              because PSS modulates Efd
                                              during fault as well)
        Peak rotor angle          : 103.3 deg
        Late rotor osc (4-8 s pp) : 12.8 deg
            vs AVR-only:            27.9 deg
            damping reduction:      ~54%

    CCT comparison (deep inductive fault Z_f = j0.10, D=3):
        GENCLS (Phase 1)           : 339 ms
        GENROU bare (Phase 2.0)    : 275 ms
        GENROU + AVR (Phase 2.1)   : 325 ms
        GENROU + AVR + PSS (this)  : 325 ms

    The PSS doesn't lift CCT (it isn't meant to — first-swing
    stability is set by the AVR's ability to force E'q during the
    fault).  The PSS's job is in the *late* window: damping the
    swing-mode oscillation that the AVR alone cannot suppress.

PSSE-vs-smib expected agreement (post 50 ms of fault):
    Rotor-angle peak, period:  ~3-5 %
    Efd peak:                  ~1-2 %
    Vpss peak:                 ~1-2 %
    Late-window damping ratio: within ~10 %
    CCT:                       ~5-10 ms

If your PSSE result differs from smib's by more than ~10 % on
damping or 15 ms on CCT, check (in order):
  1. DELT = 0.002 s
  2. MBASE = 100 MVA
  3. .dyr exactly matches smib_phase2_2.dyr
  4. Kc = 0, Kf = 0 on the ST1A
  5. PSS sign — some PSSE versions take Δω in different units;
     verify Vpss starts at 0 and goes positive when ω̄ goes positive
     during the first post-fault accelerating swing.
"""
import os, sys

PSSE_PATH = r"C:\Program Files\PTI\PSSE35\35.6\PSSPY39"
if PSSE_PATH not in sys.path:
    sys.path.append(PSSE_PATH)

import psspy
import psse35

PSSE_OUT = "smib_phase2_2_fault.out"


def main():
    psspy.psseinit(50)

    if os.path.exists("smib_phase2_2.sav"):
        psspy.case("smib_phase2_2.sav")
    elif os.path.exists("../phase2_1/smib_phase2_1.sav"):
        psspy.case("../phase2_1/smib_phase2_1.sav")
    elif os.path.exists("../phase1/smib_phase1.sav"):
        psspy.case("../phase1/smib_phase1.sav")
    else:
        build_case()
    psspy.fnsl([0, 0, 0, 1, 1, 0, 0, 0])

    psspy.dyre_new([1, 1, 1, 1], "smib_phase2_2.dyr", "", "", "")

    psspy.machine_array_channel([1,  1, 1], 101, "1 ")
    psspy.machine_array_channel([2,  2, 1], 101, "1 ")
    psspy.machine_array_channel([4,  4, 1], 101, "1 ")
    psspy.machine_array_channel([5,  5, 1], 101, "1 ")
    psspy.voltage_channel(       [6, -1,  1, 101], "Vmag_101")
    psspy.machine_array_channel([7,  7, 1], 101, "1 ")
    # PSS output Vpss is exposed via the stabiliser model channels; the
    # exact API varies by PSSE version.  Consult your model channel
    # documentation for ST1A and PSS1A.

    psspy.dynamics_solution_param_2(realar3=0.002)
    psspy.strt_2([0, 0], PSSE_OUT)
    psspy.run(0, 1.0, 0, 1, 0)

    psspy.shunt_data(101, "1 ", 1, [0.0, -10.0])
    psspy.run(0, 1.20, 0, 1, 0)

    psspy.shunt_data(101, "1 ", 0, [0.0, 0.0])
    psspy.run(0, 5.0, 0, 1, 0)

    print(f"Done.  Channels written to {PSSE_OUT}")


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
    print("Phase 2.2 case built from scratch.")


if __name__ == "__main__":
    main()
