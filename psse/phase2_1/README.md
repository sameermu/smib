# PSSE benchmark — Phase 2.1 (GENROU + ST1A)

This folder reproduces the smib Phase 2.1 scenarios in PSSE so you
can benchmark our GENROU + ST1A coupled simulation against full
PSSE GENROU + ST1A.

## What we are reproducing

Same network and operating point as Phase 1 / 2.0, now with both a
GENROU machine **and** an ST1A static exciter driving its field:

- 2-bus SMIB on a 100 MVA system base, 50 Hz
- Bus 101 = generator bus (GENROU + ST1A), Bus 102 = infinite bus
- Line: R = 0, X = 0.5 pu
- Operating point: P = 0.8 pu, Q = 0.2 pu
- Slack: V = 1.0 pu / 0 deg
- Two disturbance scenarios:
  1. **Voltage-step (small-signal)**: at t = 1 s, bump `Vref` by
     +0.02 pu and watch `|V|`, `Vc`, `Efd`, `E'q` over 10 s.
  2. **Deep inductive fault (large-signal)**: 3-phase shunt fault at
     bus 101 with `Z_f = j*0.10` pu applied at t = 1.0 s, cleared at
     t = 1.20 s.  5-second post-fault simulation horizon.

AVR parameters (smib defaults, see `smib_phase2_1.dyr`):

| Parameter | Value | Meaning |
|---|---|---|
| Tr     | 0.02 s | voltage-transducer time constant |
| Ka     | 200    | regulator gain |
| Tb     | 1.0 s  | lead-lag pole |
| Tc     | 1.0 s  | lead-lag zero (= Tb means pure-gain LL) |
| Vrmax  | +7.0   | field voltage upper limit (ceiling) |
| Vrmin  | -6.4   | field voltage lower limit |
| Kc     | 0      | rate-of-change limit (off in smib) |
| Kf, Tf | 0, 1   | rate feedback (off in smib) |

## Smib reference numbers to compare against

### Initialisation (identical PF to Phase 2.0)

| Quantity | Smib value | Notes |
|---|---|---|
| V1 from PF                    | 1.0178 pu / +23.142°  | identical to Phase 1 / 2.0 |
| Initial rotor angle delta_0   | +68.55°               | phasor behind X_q |
| Initial Eqp                   | 0.9238 pu             | transient field flux |
| Initial Edp                   | 0.4571 pu             | slow q-axis damper flux |
| Initial Efd                   | 2.0256 pu             | demanded by GENROU init |
| Initial Vref                  | 1.0279 pu             | back-solved from Efd/Ka + |V| |
| Initial Vc                    | 1.0178 pu             | = |V_t,0| |

### Voltage-step response (+2 % bump on Vref at t = 1 s, 10 s horizon)

| Quantity | Smib value | Notes |
|---|---|---|
| Pre-step |V|       | 1.0178 pu      |  |
| Final |V| (t=10 s) | 1.0330 pu      | +0.0152 pu, ~76 % of command |
| Peak Efd           | 6.02 pu        | well below Vrmax = 7 |
| Final Efd          | 3.10 pu        | +1.08 pu from pre-step |
| Time to |V| peak   | ~8.4 s         | = T'do (dominant CL time const) |

### Deep inductive fault (Z_f = j0.10, t_clear = 200 ms)

| Quantity | Smib value | Notes |
|---|---|---|
| Peak Efd during fault  | 7.0 pu         | hits Vrmax ceiling |
| Vc nadir during fault  | 0.29 pu        | sensed terminal voltage |
| E'q nadir during fault | 0.92 pu        | nearly no sag — AVR cancels it |
| Peak rotor angle       | 119.8°         | vs 114.2° AVR-off (see below) |
| **Deep-fault CCT**     | **290 ms**     | vs 240 ms AVR-off |
| AVR lift over bare GENROU | **+28 ms (+10.7%)** | the headline number |

**Note on the rotor-angle peak**: with AVR on, the rotor peak at the
*sub-critical* 200 ms clearing time is actually *larger* than AVR-off
(119.8° vs 114.2°).  This is not a bug.  The AVR preserves `E'q`
through the fault, which means more synchronising power is available
post-fault and the new equilibrium angle is slightly farther out, so
the rotor overshoots farther on its way there.  The payoff is in
the **CCT** number: the AVR-on case stays *stable* at clearing
times where the AVR-off case loses synchronism.

## CCT comparison protocol across three machine models

This is the headline summary table from Phase 2.1 §8:

| Model | CCT (ms) | Notes |
|---|---|---|
| GENCLS (Phase 1)              | 293 ms | overestimate — no E'q sag |
| GENROU bare (Phase 2.0)       | 240 ms | most conservative — full sag |
| GENROU + ST1A AVR (Phase 2.1) | 290 ms | AVR force-fields rotor flux |

GENCLS and GENROU+AVR happen to agree closely on CCT for this case,
but for very different reasons: GENCLS doesn't model the sag at all
(it has no E'q state), while GENROU+AVR models the sag and then
cancels it via the field-forcing AVR.  The agreement is a coincidence
of this particular operating point — sweep the loading and the
GENCLS curve diverges from the GENROU+AVR curve.

## The simplified ST1A caveat

Smib's Phase 2.1 ST1A is a **2-state simplified model** that omits
several details of the full IEEE 421.5 §5.1 ST1A:

- No input limiter (`Vimax` / `Vimin` on the sensed-voltage path)
- No regulator time constant (`Ta`) — smib's regulator is
  instantaneous-gain past the lead-lag
- No rate feedback (`Kf`, `Tf`) — the rate-feedback path is
  pedagogically optional and smib leaves it out
- No field-current compensation (`Kc`) — smib's `Efd` is unloaded

**What this costs in the comparison:**

- **Voltage-step transient shape**: small differences (~1 % on peak
  Efd) from the missing `Ta` lag and rate feedback.  The dominant
  time constant is still `T'do = 8 s`, which both models capture
  identically.
- **Fault-ceiling behaviour**: PSSE's `Kc` makes `Vrmax` field-current
  dependent.  With Kc=0 (set in `smib_phase2_1.dyr`), the ceiling is
  fixed at 7.0 in both models.  Agreement on peak Efd should be
  better than 1 %.
- **Rate-of-change limit**: not modelled in smib; PSSE's `Kf` rate
  feedback prevents very fast Efd swings.  Set Kf = 0 in PSSE to
  match.

**Expected agreement:** within 3-5 % on rotor-angle peak and
oscillation period; within 1-2 % on Efd peak and steady-state Vref;
within 5-10 ms on CCT.

## Files

- `smib_phase2_1.dyr` — GENROU + ST1A dynamics, smib defaults
- `run_phase2_1_fault.py` — psspy automation (deep fault scenario)

## Setup

### Option A: GUI

1. Load `smib_phase1.sav` (or rebuild — same network as Phase 1 / 2.0)
2. Solve PF → V1 should converge to 1.0178 pu / +23.142°
3. Load `smib_phase2_1.dyr` (replaces any prior dynamics)
4. Activity sequence (deep-fault scenario):
   `CONL`, `ORDR`, `FACT`, `TYSL`, `STRT`, `RUN to 1.0`,
   `SHUNT 101 1 G=0 B=-10`, `RUN to 1.20`,
   `SHUNT 101 1 G=0 B=0`, `RUN to 5.0`
5. Plot delta, omega, P, Q, Vmag, Efd, plus the AVR states (Vc, x_LL)
   and Vref — channel codes depend on PSSE build.

For the voltage-step scenario, replace step 4 with a `Vref`
increment at t=1.0 s using either the GUI's "Change Model Data"
dialog or `psspy.increment_value` with `macndx(101, "1 ", "AVR",
"VREF")`.

### Option B: psspy

Edit `PSSE_PATH` at the top of `run_phase2_1_fault.py` and run:

```
python run_phase2_1_fault.py
```

The script builds the case (or reuses an existing `.sav` if
present), loads dynamics, sets up channels, applies the deep
inductive fault, and writes `smib_phase2_1_fault.out`.

## CCT comparison protocol

To bisect on CCT in PSSE matching the smib §8 result of 290 ms:

1. Run the script repeatedly with different fault durations.
2. Stability proxy: rotor stays within ±360° of starting angle.
3. Bisect between 100 ms and 500 ms.

## Common gotchas

- **Convert loads first** (`CONL`) before dynamics initialisation.
- **MBASE = 100 MVA** must match the system base.
- **Shunt sign**: B = -10 pu is inductive (Z_f = j*0.10).
- **DELT = 0.002 s** to match smib's h = 2 ms.
- **Kc = 0** in the .dyr — smib doesn't model field-current
  compensation, so use Kc = 0 in PSSE for the cleanest comparison.
- **Kf = 0, Tf = 1.0** — rate feedback off in PSSE to match smib.
  Tf cannot be zero in some PSSE versions; setting Kf = 0 makes the
  rate path inactive regardless of Tf.
- **Ta = 0.01 s** in the .dyr — smib uses Ta = 0 (instantaneous gain)
  but some PSSE versions reject Ta = 0; the tiny lag at 0.01 s is
  far below all other time constants and has negligible effect.
