# PSSE benchmark — Phase 1 (GENCLS)

This folder contains everything needed to reproduce the smib Phase 1
deep-inductive-fault scenario in PSSE so you can benchmark our
trapezoidal integrator and swing-equation implementation against the
industry reference.

## What we are reproducing

The Phase 1 notebook simulates this scenario:

- 2-bus SMIB on a 100 MVA system base, 60 Hz
- Bus 101 = generator bus (GENCLS machine), Bus 102 = infinite bus (slack)
- Line: R = 0, X = 0.5 pu (no resistance, lossless)
- Operating point: P = 0.8 pu, Q = 0.2 pu injected at bus 101
- Slack: V = 1.0 pu / 0 deg
- Disturbance: 3-phase shunt fault at bus 101 with Z_f = j*0.10 pu
  (purely inductive, representative of transformer leakage / earth
  paths) applied at t = 1.0 s, cleared at t = 1.10 s.
- Run for 5 s after fault clearance.

## Smib reference numbers to compare against

| Quantity                              | Smib value     |
|---|---|
| V1 from PF (magnitude)                | 1.0178 pu      |
| V1 from PF (angle)                    | +23.142 deg    |
| Initial rotor angle delta_0           | +35.495 deg    |
| Internal EMF \|E'\| (held constant)   | 1.1023 pu      |
| Mid-fault \|V\| (t = 1.05 s)          | 0.348 pu       |
| Mid-fault P                           | +0.291 pu      |
| Mid-fault Q                           | +0.848 pu      |
| Peak rotor swing (post-fault)         | ~42 deg @ 4.3 s|
| Bolted-fault CCT (smib bisection)     | 191 ms         |
| Bolted-fault CCT (analytic EAC)       | 187 ms (2.3% off) |
| Deep inductive CCT (Z_f = j*0.10)     | 293 ms         |

PSSE results within ~3 % of these numbers should be considered a clean
match — the only expected differences are integrator step size and
exactly how the shunt fault is applied.

## Files

- `smib_phase1.dyr` — GENCLS dynamics model, Bus 101, machine ID '1 '.
  Parameters: H = 4.0 s, D = 0, X'd = 0.30 pu.
- `run_phase1_fault.py` — psspy automation script that loads the case,
  applies the fault, runs the simulation, and writes a channel file
  with delta, omega, P, Q, V, Efd traces for comparison.

## Setup — option A: PSSE GUI

If you prefer to drive PSSE interactively rather than via Python:

1. **Build the power-flow case.**
   - File → New
   - Add Bus 101 named "GEN BUS", base kV 18.0, code 1 (PQ)
   - Add Bus 102 named "INF BUS", base kV 18.0, code 3 (slack), V = 1.0
   - Add Generator at 101: MBASE = 100 MVA, P = 80 MW, Q = 20 MVAR.
     Leave X_source at 0.30 pu (matches X'd).
   - Add Branch 101 → 102: R = 0, X = 0.5 pu, MVA base 100
   - Solve power flow (full Newton).  V1 should converge to 1.0178 pu / +23.142 deg.
   - Save as `smib_phase1.sav` in this folder.

2. **Add dynamics.**
   - File → Open → `smib_phase1.dyr` (this folder)
   - Convert loads (none in this case but PSSE expects you to do it)
   - Order admittance matrix and factorise (CONL, ORDR, FACT, TYSL)

3. **Run the disturbance.**
   - In the Activity bar: `STRT` to initialise dynamics
   - `RUN` to t = 1.0 s
   - Apply shunt fault: `SHUNT 101 1` with G = 0, B = -10
     (B = -10 pu corresponds to admittance Y_f = -j*10, i.e. impedance
     Z_f = j*0.10 pu — this is the inductive shunt fault we want)
   - `RUN` to t = 1.10 s
   - Clear fault: `SHUNT 101 1` with G = 0, B = 0
   - `RUN` to t = 5.0 s
   - Plot channels: delta (deg), omega (pu), Vmag (pu), P (pu), Q (pu)

4. **Compare.**
   - Open the channel file in PSSPLT or Plot Editor
   - Cross-check the smib reference numbers in the table above.

## Setup — option B: psspy automation

If your PSSE install has the Python bindings:

1. Edit `run_phase1_fault.py` and update `PSSE_PATH` at the top to
   match your install (default points at `C:\Program Files\PTI\PSSE35\...`).
2. From a terminal in this folder:

   ```
   python run_phase1_fault.py
   ```

3. The script will:
   - Build the 2-bus case in memory (or load `smib_phase1.sav` if it
     exists)
   - Solve power flow
   - Load the dynamics file
   - Set up channels for delta, omega, P, Q, Vmag, Efd
   - Run pre-fault to 1.0 s, apply the shunt fault, run during-fault
     to 1.10 s, clear the fault, run post-fault to 5.0 s
   - Write `smib_phase1_fault.out` for comparison

## CCT comparison

To bisect on critical clearing time in PSSE (matching smib's §9):

1. Run the script repeatedly with different fault durations
   (`FaultDuration` variable in the activity sequence).
2. Stability proxy: rotor stays within ±360° of starting angle.
3. Bisect between 100 ms and 1000 ms.
4. Compare to smib's bracket of 293 ms for Z_f = j*0.10, or 191 ms
   for the bolted limit (Z_f = j*0.001).

## Common gotchas

- **Convert loads first.** PSSE's `CONL` activity must run before
  dynamics initialisation, even with zero load — it sets up the
  algebraic network for time-domain integration.
- **MBASE matches.** GENCLS parameters (H, X'd) are on machine MVA
  base. If MBASE ≠ system base, PSSE will silently rescale, and
  results won't match smib (which uses per-unit on the same base
  throughout).
- **Shunt sign.** PSSE shunt B is positive for capacitive, negative
  for inductive.  For a Z_f = j*0.10 pu fault, B = -10 pu (admittance
  -j*10, impedance j*0.10).
- **Step size.** smib uses h = 2 ms.  PSSE's default `DELT` is 1/4
  cycle ≈ 4.17 ms at 60 Hz.  Set `DELT = 0.002` or use the
  `dynamics_solution_param_2(realar3=0.001)` call in the psspy
  script to match smib.
