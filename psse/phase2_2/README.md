# PSSE benchmark — Phase 2.2 (GENROU + ST1A + PSS1A)

Same network and operating point as Phase 2.1 — now with a PSS1A
stabiliser in parallel to the AVR's voltage error.

## What we are reproducing

- 2-bus SMIB, 100 MVA system base, 50 Hz
- GENROU at bus 101 (Kundur Table 4.2, D=3 load damping)
- ST1A AVR (Ka=200, Tb=20 for transient gain reduction)
- PSS1A stabiliser (Tw=5 s, T1=T3=0.5, T2=T4=0.05, Ks=20,
  Vstmax=±0.10)
- Deep inductive fault Z_f = j*0.10 at bus 101 from t=1.0 s to
  t=1.20 s.

PSS parameters in `smib_phase2_2.dyr`:

| Parameter | Value | Meaning |
|---|---|---|
| Vstmin / Vstmax | ±0.10  | output limits |
| Tw  | 5.0 s  | washout time constant |
| T1  | 0.50 s | first lead zero |
| T2  | 0.05 s | first lag pole |
| T3  | 0.50 s | second lead zero |
| T4  | 0.05 s | second lag pole |
| Ks  | 20     | gain |

## Smib reference numbers to compare against

### Initialisation (identical to Phase 2.1)

| Quantity | Smib value |
|---|---|
| V1 from PF                     | 1.0178 pu / +23.142° |
| delta_0                        | +68.55° |
| Eqp_0 / Edp_0                  | 0.9238 / 0.4571 pu |
| Efd_0                          | 2.0256 pu |
| Vref                           | 1.0279 pu |
| Vc                             | 1.0178 pu |
| x_w, x_LL1, x_LL2              | all 0 (PSS dormant at SS) |
| Vpss_0                         | 0 pu |

### Deep inductive fault (Z_f = j0.10, t_clear = 200 ms)

| Quantity | Smib value | Notes |
|---|---|---|
| Peak Efd                | 7.0 pu     | hits Vrmax ceiling |
| Peak \|Vpss\|           | 0.10 pu    | hits ±Vstmax limit |
| E'q nadir               | 0.91 pu    | slightly lower than AVR-only (Vpss also rides Efd through fault) |
| Peak rotor angle        | 103.3°     | essentially same as AVR-only |
| **Late osc 4-8 s pp**   | **12.8°**  | **vs 27.9° AVR-only — 54% reduction** |
| Mean iters/step         | 2.5        | up from ~2.3 for AVR-only |

### CCT — four-way comparison (D=3 load damping)

| Model | CCT (ms) | vs GENCLS |
|---|---|---|
| GENCLS (Phase 1) | 339 | reference |
| GENROU bare (Phase 2.0) | 275 | −64 ms |
| GENROU + ST1A AVR (Phase 2.1) | 325 | −14 ms |
| GENROU + ST1A + PSS1A (Phase 2.2) | 325 | −14 ms |

The PSS does not lift CCT — it isn't meant to.  First-swing
stability is set by how much accelerating power the AVR force-
fields away during the fault.  The PSS's job is in the **post-fault
window**: damping the swing-mode oscillation that the AVR + load
damping alone cannot suppress.

## Caveats vs full PSSE PSS1A

Smib's PSS1A is the simplified one-input form.  Specifically
missing vs the full IEEE 421.5 §6.1 PSS1A:

- No input low-pass filter on Δω (PSSE sometimes adds a small Tr
  on the speed measurement).
- Only one washout stage (some PSSE variants double up).
- No multi-band lead-lag (single-band is fine for this SMIB swing
  mode at ~0.7 Hz).
- The 7-state PSS2A/PSS2B (dual-input, accel/speed combined) is
  not modelled — those land in a later phase if relevant.

**Expected agreement** with full PSSE PSS1A:
- Late-window damping ratio: within ~10 %
- Rotor-angle peak: within 3-5 %
- Vpss sign and timing: agreement at the kink-points (Vpss
  changes sign whenever Δω crosses zero — same in any
  implementation).
- CCT: within 5 ms.

## Files

- `smib_phase2_2.dyr` — GENROU + ST1A + PSS1A dynamics
- `run_phase2_2_fault.py` — psspy automation

## Setup

Same overall flow as Phase 2.1.  Replace the dynamics file with
`smib_phase2_2.dyr` and run the same activity sequence.

## CCT comparison protocol

Bisect on `t_clear` until rotor stays within ±360°.  Expect 325 ms
±10 ms.

## Common gotchas

- **Vpss sign convention** varies across PSSE versions / .dyr
  formats.  Smib uses the convention Δω > 0 (rotor accelerating) →
  Vpss > 0 (boost Efd to slow the rotor down via the damping
  torque).  Verify by checking that Vpss is positive in the first
  ~200 ms post-fault clear (when Δω is positive from the residual
  acceleration).
- **Kc = 0** on the ST1A (field-current compensation off; smib
  doesn't model it).
- **Tf cannot be 0** in some PSSE versions — leave Tf = 1.0, Kf = 0
  to disable the rate feedback.
- **DELT = 0.002** to match smib's h = 2 ms.
- **D = 3** on the GENROU — load damping, must match what the
  notebook uses or post-fault settling will diverge.
