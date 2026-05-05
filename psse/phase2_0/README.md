# PSSE benchmark — Phase 2.0 (GENROU bare)

This folder reproduces the smib Phase 2.0 deep-inductive-fault
scenario in PSSE so you can benchmark our 4-state two-axis transient
GENROU model against the full PSSE GENROU.

## What we are reproducing

Same scenario as Phase 1 (same network, same operating point, same
fault), but with a **GENROU machine** instead of GENCLS:

- 2-bus SMIB on a 100 MVA system base, 60 Hz
- Bus 101 = generator bus (GENROU machine), Bus 102 = infinite bus
- Line: R = 0, X = 0.5 pu
- Operating point: P = 0.8 pu, Q = 0.2 pu
- Slack: V = 1.0 pu / 0 deg
- Disturbance: 3-phase shunt fault at bus 101 with Z_f = j*0.10 pu
  applied at t = 1.0 s, cleared at t = 1.10 s.
- 5-second post-fault simulation horizon.

Machine parameters (Kundur Table 4.2, round rotor thermal unit) —
see `smib_phase2_0.dyr` for the exact values:

| Parameter | Value | Meaning |
|---|---|---|
| H        | 4.0 s   | inertia |
| D        | 0       | mechanical damping |
| Xd       | 1.81 pu | d-axis synchronous reactance |
| Xq       | 1.76 pu | q-axis synchronous reactance |
| X'd      | 0.30 pu | d-axis transient reactance |
| X'q      | 0.65 pu | q-axis transient reactance |
| X"d      | 0.23 pu | d-axis sub-transient reactance |
| Xl       | 0.16 pu | stator leakage |
| T'd0     | 8.0 s   | d-axis transient open-circuit time constant |
| T"d0     | 0.03 s  | d-axis sub-transient OC time constant |
| T'q0     | 1.0 s   | q-axis transient OC time constant |
| T"q0     | 0.07 s  | q-axis sub-transient OC time constant |
| S(1.0)   | 0.13    | quadratic saturation at 1.0 pu air-gap flux |
| S(1.2)   | 0.50    | quadratic saturation at 1.2 pu air-gap flux |

## Smib reference numbers to compare against

| Quantity | Smib value | Notes |
|---|---|---|
| V1 from PF                     | 1.0178 pu / +23.142°       | identical to Phase 1 |
| Initial rotor angle delta_0    | +68.55°                    | phasor behind X_q |
| Initial Eqp                    | 0.9238 pu                  | transient field flux |
| Initial Edp                    | 0.4571 pu                  | slow q-axis damper flux |
| Initial Efd (=AVR setpoint)    | 2.026 pu                   | held constant — no AVR yet |
| Mid-fault \|V\| (t=1.05 s)       | 0.305 pu                   |  |
| Mid-fault P                    | +0.329 pu                  |  |
| Mid-fault Q                    | +0.602 pu                  | less than GENCLS (0.847) — Eqp sag |
| Mid-fault Eqp                  | 0.9088 pu                  | sagged 0.015 from initial |
| Peak rotor angle               | ~73° around t = 1.3 s post-fault |  |
| **Deep-fault CCT (j*0.10)**    | **240 ms**                 | vs GENCLS 293 ms — Eqp sag costs 53 ms |
| Bolted-fault CCT (j*0.001)     | 153 ms                     | vs GENCLS 191 ms — Eqp sag costs 38 ms |

## The 4-state vs full GENROU caveat

Smib's Phase 2.0 GENROU is a **simplified 4-state model** that omits
sub-transient detail.  Specifically missing vs full PSSE GENROU:

- Sub-transient flux states `psidpp` (ψ"_d) and `psiqpp` (ψ"_q)
- Sub-transient reactances X"_d and X"_q in the output equations
- Sub-transient time constants T"_d0 and T"_q0
- Stator leakage reactance Xl in the compensated-current terms
- Saturation acts on |Eqp| only, not the full air-gap flux

**What this costs in the comparison:**

- **First ~30-50 ms of any fault**: PSSE's apparent reactance is X"_d
  (≈ 0.23) vs smib's X'_d (≈ 0.30).  PSSE shows ~30 % bigger initial
  fault-current spike.  By t = 50 ms post-fault start, both models
  are in the transient regime and converge.
- **Beyond 50 ms**: the swing dynamics are dominated by T'_d0 and
  T'_q0 which we DO model.  Rotor angle peak, oscillation period,
  and CCT match PSSE within ~3-5 %.
- **Saturation at high field forcing**: smib only saturates Eqp
  directly; PSSE saturates the full air-gap flux.  For our scenario
  with Efd held at its initial value of 2.0 (no AVR), saturation is
  barely active in either model.  Becomes more relevant in Phase 2.1
  when the AVR drives Efd to the ceiling.

**Expected agreement:** within 3-5 % on rotor angle, period, and
CCT.  Within 10 % on the headline canonical-5 traces beyond the
first 50 ms.  PSSE's first-50-ms fault-current spike will be
~30 % bigger.

## Files

- `smib_phase2_0.dyr` — GENROU dynamics with Kundur Table 4.2 params
- `run_phase2_0_fault.py` — psspy automation script

## Setup

Same overall flow as Phase 1.  The PF case is identical — you can
reuse `phase1/smib_phase1.sav` if you've saved it.  Only the
dynamics file changes.

### Option A: GUI

1. Load `smib_phase1.sav` (or rebuild the case, see Phase 1 README)
2. Solve PF → V1 should converge to 1.0178 pu / +23.142°, identical
   to Phase 1 (same network, same loading)
3. Replace the dynamics: load `smib_phase2_0.dyr` instead of
   `smib_phase1.dyr`
4. Activity sequence: `CONL`, `ORDR`, `FACT`, `TYSL`, `STRT`, `RUN to 1.0`,
   `SHUNT 101 1 G=0 B=-10`, `RUN to 1.10`, `SHUNT 101 1 G=0 B=0`,
   `RUN to 5.0`
5. Plot delta, omega, P, Q, Vmag, Efd, plus Eqp and Edp from the
   GENROU model (channel codes 26 and 27 in PSSE 35; check your
   manual for the exact API)

### Option B: psspy

Edit `PSSE_PATH` at the top of `run_phase2_0_fault.py` and run:

```
python run_phase2_0_fault.py
```

The script builds the case (or reuses `smib_phase2_0.sav` if present),
loads the GENROU dynamics, sets up channels, applies the fault, and
writes `smib_phase2_0_fault.out`.

## CCT comparison protocol

To bisect on CCT in PSSE matching the smib §9 result of 240 ms for
the deep inductive fault:

1. Run the script repeatedly with different fault durations.
2. Stability proxy: rotor stays within ±360° of starting angle.
3. Bisect between 100 ms and 1000 ms.

## Common gotchas

- **Convert loads first** (`CONL`) before dynamics initialisation.
- **MBASE = 100 MVA** must match the system base, otherwise PSSE
  silently rescales H, X, etc.
- **Shunt sign**: B = -10 pu is inductive (Z_f = j*0.10), B = +10 is
  capacitive.  The fault model in this scenario is purely inductive.
- **DELT = 0.002 s** to match smib's h = 2 ms integration step.
- **PSSE GENROU assumes X"_d = X"_q** for round-rotor.  In smib we
  don't model sub-transient at all in the 4-state form, so this
  doesn't directly apply, but it's worth knowing if you ever upgrade
  smib to a full 6-state GENROU.
