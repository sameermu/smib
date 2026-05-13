# PSSE benchmarks

Each smib phase ships a matching PSSE setup so you can independently
verify our trapezoidal integrator and per-model implementations
against the industry reference.

## Layout

- `phase1/` — GENCLS, deep inductive fault.  See `phase1/README.md`.
- `phase2_0/` — GENROU bare (no AVR/PSS/Gov), same fault scenario
  as Phase 1.  See `phase2_0/README.md`.
- `phase2_1/` — GENROU + ST1A AVR, voltage-step + deep fault scenarios.
  See `phase2_1/README.md`.
- `phase2_2/` — GENROU + ST1A + PSS1A (coming once Phase 2.2 lands).
- `phase2_3/` — GENROU + ST1A + PSS1A + TGOV1 (coming once Phase 2.3 lands).

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
- Smib runs at h = 2 ms; PSSE defaults to ~5 ms (1/4 cycle).  Set
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

## Planned companion — Dynawo (Phase 3 onward)

Once Phase 3 (IBR — REGC_A + REEC_A + REPC_A) lands we will add a
sibling `dynawo/` folder at the same level as `psse/`, using the
same per-phase template (one input deck, one run script, one README
with the smib reference numbers).  Dynawo (RTE + AIA) is built on
Modelica so its IBR libraries are *readable* — making it the
natural open-source complement to PSSE for inverter studies where
PSSE alone is unsatisfying.  Concrete role: independent
tie-breaker, read-the-equations sanity check on IBR limit logic,
and a credible reference for colleagues without a PSSE license.

We are *not* setting Dynawo up retroactively for Phase 1/2.  For
classical machines, smib + PSSE is the correctness floor that pays
the rent.  Phase 3 is where the second open-source leg earns its
keep.  See the "External-tool benchmarks" section of the top-level
`README.md` for the broader strategy.
