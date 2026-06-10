# PSSE benchmark — Phase 2.3 (GENROU + ST1A + PSS1A + TGOV1)

Same network and operating point as Phase 2.1/2.2 — now with a
TGOV1 turbine governor closing the mechanical loop.  Mechanical
power, a constant in every earlier phase, is now a live dynamic
channel (PMECH).

## What we are reproducing

- 2-bus SMIB, 100 MVA system base, 50 Hz
- GENROU at bus 101 (Kundur Table 4.2, D=3 load damping)
- ST1A AVR (Ka=200, Tb=20 for transient gain reduction)
- PSS1A stabiliser (Tw=5 s, T1=T3=0.5, T2=T4=0.05, Ks=20)
- TGOV1 governor (R=0.05, T1=0.5, T2=2.5, T3=7.5, Vmax=1.0, Vmin=0.0)
- Deep inductive fault Z_f = j*0.10 at bus 101 from t=1.0 s to
  t=1.20 s.

TGOV1 parameters in `smib_phase2_3.dyr` (PSSE CON order):

| CON | Parameter | Value | Meaning |
|---|---|---|---|
| J   | R    | 0.05  | droop (pu speed / pu power, machine base) |
| J+1 | T1   | 0.5 s | valve/servo time constant |
| J+2 | VMAX | 1.0   | max valve opening |
| J+3 | VMIN | 0.0   | min valve opening |
| J+4 | T2   | 2.5 s | reheat lead (HP fraction T2/T3 = 1/3) |
| J+5 | T3   | 7.5 s | reheat time constant |
| J+6 | Dt   | 0.0   | turbine damping |

## Smib reference numbers to compare against

### Initialisation (PF identical to Phase 2.0/2.1/2.2)

| Quantity | Smib value |
|---|---|
| V1 from PF                     | 1.0178 pu / +23.142° |
| delta_0                        | +68.55° |
| Eqp_0 / Edp_0                  | 0.9238 / 0.4571 pu |
| Efd_0                          | 2.0256 pu |
| Vref                           | 1.0279 pu |
| **Pref (= Pm_0 = Pe_0)**       | **0.8000 pu** |
| **x_v_0 / x_r_0**              | **0.8000 / 0.8000 pu** |
| PSS states                     | all 0 (dormant at SS) |

The governor back-solve is a unity pass-through: Pref = Pm_0
directly.  If your PSSE initialisation report shows PMECH ≠ 0.8 or
a suspect GREF, check the droop wiring — Pref through the 1/R gain
is the classic factor-of-20 bug.

### Deep inductive fault (Z_f = j0.10, t_clear = 200 ms)

| Quantity | Smib value | Notes |
|---|---|---|
| Peak rotor angle        | 102.5°    | vs 103.3° gov-off — valve barely moves in 200 ms |
| Peak Efd                | 7.0 pu    | hits Vrmax (unchanged) |
| PMECH during fault      | few-% dip | droop opposes rotor acceleration |
| PMECH after clearing    | small phase-lagged wiggles opposing slip | NOT a damping mechanism — see notebook §8 |
| CCT                     | 333 ms    | +6 ms vs Phase 2.2 — see below |

### CCT — five-way comparison (D=3 load damping)

| Model | CCT (ms) |
|---|---|
| GENCLS (Phase 1) | 339 |
| GENROU bare (Phase 2.0) | 275 |
| GENROU + ST1A AVR (Phase 2.1) | 326 |
| GENROU + ST1A + PSS1A (Phase 2.2) | 326 |
| GENROU + ST1A + PSS1A + TGOV1 (Phase 2.3) | 333 |

The headline check is a **bounded result**: the governor may move
CCT only by a few ms and only *upward*.  Transient stability is
decided in well under a second, but a near-CCT fault lasts 330+ ms
— long enough for the T1 = 0.5 s valve servo to track a meaningful
fraction of the droop command (which eases the valve off while
slip is positive) and for the HP fraction (T2/T3 = 1/3) to pass it
to the shaft.  Smib measures +6.4 ms.  If your PSSE run shows the
Phase 2.3 CCT *below* Phase 2.2, or more than ~15 ms above it,
suspect the TGOV1 .dyr entry (usually VMAX/VMIN or R units).

### Primary frequency response (smib-side check)

The droop steady state is validated in smib via the infinite-bus
angle-ramp scenario (`grid_frequency_ramp_schedule`), which has no
direct single-machine PSSE equivalent — in a PSSE 2-bus case the
swing bus pins the frequency, so a governor never sees a sustained
speed deviation.

| Quantity | Smib value |
|---|---|
| Grid event                  | Δf = −0.2 Hz (Δω̄ = −0.004 pu) ramped over 2 s |
| Final rotor slip            | −0.00400 pu (locks to grid exactly) |
| ΔPm at t = 40 s             | +0.0796 pu |
| ΔPm asymptotic              | +0.0800 pu = −Δω̄/R exactly |
| Valve                       | 0.80 → 0.88 pu |
| Measured droop (4-pt sweep) | R = 0.0500 (slope fit of ΔPm vs ω̄) |

To reproduce a governor frequency response in PSSE, the standard
trick is to replace the infinite bus with a large finite machine
(GENCLS, H = 50 s on a 10 GVA base) and apply a load step at bus
101 — the system frequency then genuinely sags and TGOV1 picks up
per its droop.  That is left as an optional extension; the analytic
droop relation is exact by construction and smib's measured slope
matches R to 4 decimal places.

## Caveats vs full PSSE TGOV1

Smib's TGOV1 matches the PSSE block structure (droop → servo with
non-windup valve limits → reheat lead-lag → Dt term).  Differences:

- PSSE applies the speed input as SPEED deviation in pu on system
  frequency; smib's ω̄ is identical by construction (both are slip).
- PSSE TGOV1 droop R is entered on the machine MBASE.  With
  MBASE = 100 MVA = system base here, no conversion arises.  If you
  rescale MBASE, convert R accordingly.
- Smib's non-windup valve limit freezes the servo state at the
  boundary (same as PSSE's non-windup convention); the trapezoidal
  corrector can overshoot the limit by O(h·dx/dt) on the crossing
  step — invisible at h = 2 ms.

## Running

GUI: load the Phase 1 `.sav`, read `smib_phase2_3.dyr` into a new
snapshot, add the seven channels (ANGLE, PELEC, **PMECH**, ETERM,
EFD, bus voltage, SPEED), set DELT = 0.002, run to 1.0 s, apply the
bus shunt B = −10 pu at bus 101, run to 1.20 s, remove it, run to
10 s.

psspy: `python run_phase2_3_fault.py` from this folder.
