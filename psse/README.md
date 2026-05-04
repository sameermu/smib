# PSSE benchmarks

Each smib phase ships a matching PSSE setup so you can independently
verify our trapezoidal integrator and per-model implementations
against the industry reference.

## Layout

- `phase1/` — GENCLS, deep inductive fault.  See `phase1/README.md`.
- `phase2/` — GENROU + ST1A + PSS1A + TGOV1.  Built incrementally as
  Phase 2 lands; see `phase2/README.md` once it exists.

## What to compare

Every phase ships:

1. **A `.dyr` dynamics file** with the model parameters in PSSE syntax.
2. **A psspy automation script** (`run_phaseN_fault.py`) that drives
   the same disturbance the smib notebook simulates.
3. **A README** with both GUI and psspy instructions, plus a table of
   smib reference numbers (V from PF, peak swing, CCT) for quick
   cross-check.

The expected agreement between smib and PSSE is within ~3 % on all
canonical traces, given:

- Both use the implicit trapezoidal integrator.
- Both use the same swing equation / model-equation set.
- Smib runs at h = 2 ms; PSSE defaults to ~4.17 ms (1/4 cycle).  Set
  PSSE's `DELT = 0.002` to match.

## Why we benchmark

Two distinct goals:

1. **Correctness floor.**  If smib disagrees with PSSE on the same
   setup, one of the two has a bug.  The disagreement points us at
   exactly which model, which init step, or which integrator detail
   to debug.
2. **Pedagogical confidence.**  Colleagues from the power-systems
   side trust PSSE.  Showing that smib reproduces PSSE on a curated
   case lets them trust smib too — and then use smib as a
   transparent way to *understand* what PSSE is doing under the hood.
