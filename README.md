# smib

A transparent RMS/phasor-domain Single-Machine Infinite-Bus simulator,
built in NumPy. Every equation is in the code; nothing hides behind a
black box. Intended for teaching, benchmarking against PSSE, and
sparring with the question "can an LLM-authored transparent stack match
an industry simulator on a curated case?"

## Project status

Update this section as phases land.  Last update: Phase 2.0 in flight,
GENROU model + foundations doc + PSSE Phase 1 benchmark complete.

| Phase | Scope | Status |
|---|---|---|
| **0** — Skeleton | repo, solver, network, PF, scenarios, plotting, tests | ✅ done |
| **1** — GENCLS | classical machine on SMIB, full notebook + PSSE benchmark | ✅ done |
| **2.0** — GENROU bare | round-rotor 4-state two-axis transient model | 🔄 model done, simulator + notebook pending |
| **2.1** — + ST1A AVR | static exciter, demonstrates LV reactive support boost | ⏳ pending |
| **2.2** — + PSS1A | power system stabiliser, headline "with vs without PSS" plot | ⏳ pending |
| **2.3** — + TGOV1 | turbine governor, primary frequency response | ⏳ pending |
| **3** — IBR generic | REGC_A + REEC_A + REPC_A, weak-grid SCR sweep | ⏳ pending |
| **4** — IBR grid-forming | REGFM_A1, side-by-side with Phase 3 | ⏳ pending |
| **5** — Reactive support | SVSMO3 SVC, CSTATT STATCOM, SynCon | ⏳ pending |
| **6** — Scenario browser | dedicated ipywidgets explorer over all phases | ⏳ pending |
| **7** — Small-signal | numerical eigenvalues + Prony analysis | ⏳ pending |

**What "done" means for a phase**: model + tests + colleague-facing
notebook + PSSE benchmark folder all landed and pushed, with all 8
pedagogy rules applied.

**Phase 2.0 inner status** (since it's where we are):

- ✅ `smib/models/genrou.py` — 4-state two-axis transient model, init
  passes self-consistency at machine epsilon, flat-line drift < 4e-14°
  over 5 s.
- ✅ `docs/genrou_physical_foundations.md` — full physical explanation
  of every parameter with two SVG equivalent-circuit diagrams.
- ✅ `smib/models/gencls.py` — Iq/Id sign convention fixed to standard
  PSSE/Kundur (q-axis aligned with E', d-axis 90° behind).  All 8
  Phase 1 tests still pass.
- ⏳ `run_smib_genrou()` simulator helper with the 2×2 saliency-aware
  algebraic network solve.
- ⏳ `tests/test_genrou.py` 5-test battery.
- ⏳ `notebooks/phase2_0_genrou.ipynb` colleague-facing artefact.
- ⏳ `psse/phase2_0/` benchmark folder (parallel to `psse/phase1/`).
- ⏳ Phase 1 notebook narrative refresh — the dq sign convention
  changed under it, so §8.1 still reads "Iq surges" when it should
  now read "Id surges" (physics unchanged, label only).

## What this is not

- Not fast. State is held as named dicts for readability.
- Not a general N-bus package. The network is hard-coded as two buses
  (generator + infinite bus). The Ybus machinery generalises, but the
  scenarios do not.
- Not a replacement for PSSE. It is a reference implementation for
  model authors and students.

## Install

```bash
pip install numpy plotly ipywidgets pytest nbformat
```

From this folder:

```bash
pytest tests/
```

The frontend is a Jupyter notebook per phase, with embedded `ipywidgets`
sliders driving live Plotly plots via `smib.plotting.scenario_slider`.
No web server, no separate UI process — the notebook *is* the
dashboard. See `notebooks/phase1_gencls.ipynb` for the pattern.

## Repository layout (after Phase 1)

- `smib/powerflow.py` — 2-bus Newton–Raphson.
- `smib/network.py` — Ybus, shunt fault, slack-voltage step.
- `smib/solver.py` — implicit trapezoidal with fixed-point corrector.
- `smib/scenarios.py` — fault / voltage step / setpoint step schedules.
- `smib/simulator.py` — time-loop runner, Norton-augmented network solve.
- `smib/plotting.py` — Plotly traces, ride-through validator, slider helper.
- `smib/models/base.py` — abstract model class (state as dict).
- `smib/models/gencls.py` — classical machine, 2 states (delta, omega).
- `smib/models/pll.py` — MHI PSCAD SRF-PLL, 2 states (theta_pll, x_I).
- `tests/test_flatline.py` — power-flow & network round-trip.
- `tests/test_gencls.py` — flat-line, analytic eigenvalue, CCT bracket.
- `notebooks/phase1_gencls.ipynb` — Phase 1 colleague-facing artefact.

## Build order (details in the plan doc)

1. **Phase 1 — DONE** — GENCLS on SMIB (swing equation, CCT).
2. **Phase 2** — GENROU + ST1A + PSS1A + TGOV1 (classical full stack).
3. **Phase 3** — IBR generic (REGC_A + REEC_A + REPC_A).
4. **Phase 4** — Grid-forming (REGFM_A1), same test battery as Phase 3.
5. **Phase 5** — Reactive support (SVSMO3 SVC, CSTATT STATCOM, SynCon).
6. **Phase 6** — Scenario browser (ipywidgets in a Jupyter notebook).
7. **Phase 7** — Small-signal analysis (numerical eigenvalues + Prony).

Each phase delivers a Jupyter notebook as the colleague-facing artefact,
imports from a clean `smib.*` package, and must pass the correctness
floor before merging.

## Standard plots (every phase, every scenario)

Every simulation always plots these five canonical traces, regardless of
the technology under test:

- **P** — active power at the terminal (pu)
- **Q** — reactive power at the terminal (pu)
- **|V|** — terminal voltage magnitude (pu)
- **Id** — current d-axis component (pu, device dq frame)
- **Iq** — current q-axis component (pu, device dq frame)

Additional traces are added depending on the technology being assessed:

| Technology | Supplementary traces |
|---|---|
| GENCLS / GENROU + AVR + PSS + Gov | δ, ω, Efd, field current If, Pm, rotor fluxes E'q, E'd, ψ"d, ψ"q |
| Grid-following IBR (REGC + REEC + REPC) | PLL angle error, Iqinj, LVPL output, ride-through flag, limit flags, Pqflag |
| Grid-forming IBR (REGFM) | Internal EMF magnitude & angle, internal frequency, Qref tracking error, current-limit flag |
| SVC (SVSMO3) | Susceptance B, slope-characteristic operating point |
| STATCOM (CSTATT) | Iq command vs delivered Iq, saturation flag, DC-link voltage |
| SynCon (GENROU + ESST4B + limiters) | Efd, field current vs OEL thermal integrator, OEL/UEL/SCL flags |

## Ride-through validation rule (non-negotiable)

For any **undervoltage** event (voltage step down, fault clearing into a
depressed voltage, any scenario with |V| below nominal), **Q and Iq
must be increasing** in the direction that injects more reactive power
into the grid during the LV window. If they move the other way, the
model or the flags are wrong and the scenario cannot be trusted.

Mirror rule: for **overvoltage** events, Q and Iq must be **decreasing**
(absorbing reactive; for IBRs this means Iqinj flipping sign into HVRT
absorb-mode).

This rule is a sharper canary than the flat-line test and catches:

- Sign errors on controller feedback paths
- Reference-frame mismatches (global DQ vs device dq)
- LVRT / HVRT flag inversion bugs
- Current-limit anti-windup failures (P-priority strategies starving Iq)
- Droop-sign errors on GFM Q-V loops

A helper `validate_ride_through(t, V, Q, Iq, event_window, kind="LV")`
lives alongside the plotting module; every ride-through notebook calls
it at the end and prints PASS / FAIL with offending timestamps.

## Notebook pedagogy rules (non-negotiable)

Every phase notebook in `notebooks/` must hand-hold the colleague
through the math, not just call library functions. The primary metric
for the project is pedagogical clarity, and a notebook that just
invokes black boxes is a script, not a teaching artefact. Four rules:

1. **Show the equations before the call.** For every library function
   the notebook invokes (`two_bus_pf`, `gencls.initialise`,
   `run_smib_gencls`, ...), the cell above it must derive or restate
   the governing equation in LaTeX, motivate the inputs, and explain
   what the output means physically.

2. **Hand-hold iterative algorithms step-by-step.** When the notebook
   uses Newton-Raphson, fixed-point, bisection, FFT, or any iterative
   procedure, expose each iteration explicitly — print the state, the
   mismatch/residual, and the Jacobian or update direction at every
   step. Don't hide the loop inside a black-box call. For NR power
   flow specifically this means calling `two_bus_pf(..., verbose=True)`
   so the iteration table prints.

3. **Link state equations explicitly to every plot.** After every
   dynamics plot (fault, V step, setpoint step) include a markdown
   cell that walks through what each trace shows in terms of the state
   equations. For GENCLS this means tying $\delta(t)$ and
   $\bar\omega(t)$ back to $2H\,d\bar\omega/dt = P_m - P_e - D\bar\omega$
   and $d\delta/dt = \omega_0\bar\omega$ — what term dominates during
   the fault, what term dominates after clearing, why the oscillation
   has the period it does, why damping looks the way it does. Annotate
   the plots with vertical markers ("fault on", "fault clear", "first
   peak") and shaded regions where useful.

4. **Cross-check encapsulated methods against hand calculation.** When
   `model.initialise()` (or any black box) is called, the cell above
   must do the same calculation step-by-step, then the cell below must
   print BOTH the manual values and the encapsulated values
   side-by-side so the colleague can verify they match.

5. **Render an SLD after every load flow.** Every phase notebook
   includes a PSSE-style single-line diagram (`smib.plotting.plot_sld`)
   immediately after the power flow converges, showing each bus with
   $|V|$ and $\angle V$, the line with its impedance, $P$ and $Q$ flow
   with arrow direction, and the appropriate generator/slack symbols.
   This makes the operating point readable at a glance so colleagues
   never have to mentally reconstruct topology from a wall of code.

6. **Derive the algebraic power-transfer equations early and refer
   back to them.** For SMIB with a synchronous machine the canonical
   pair is $P_e = (|E'||V|/X_{\text{tot}})\sin\delta$ and
   $Q_e = (|E'|^2 - |E'||V|\cos\delta)/X_{\text{tot}}$, derived from
   $S = V \cdot I^*$ with $I = (E' - V)/(jX_{\text{tot}})$. Plot the
   $P$-$\delta$ curve with the operating point and the unstable
   equilibrium ($\pi - \delta_0$) marked. Compute the synchronising
   coefficient $K_s$, the maximum power $P_e^{\max}$, the two
   equilibria, and use the **equal-area criterion** to predict CCT
   analytically. Then compare the analytic CCT to the simulator's
   bracket — agreement to a few percent confirms the swing equation,
   the network, the fault application, and the integrator are
   internally consistent. For IBR and GFM phases, derive the
   equivalent algebraic relations (capability curves, droop lines)
   in the same style.

7. **Include a solver-flow visualisation in every phase notebook.**
   Add an inline SVG flowchart showing where the integrator evaluates
   `model.derivatives()` and where the trapezoidal rule adds those
   derivatives back into the state vector. Use the colour scheme from
   Phase 1 §10: orange = compute dx/dt, green = update state, blue =
   solve algebraic network, neutral = housekeeping. The flowchart
   must NOT reference Python file/line numbers — it should read as
   physics + numerics, not source-code archaeology. For Phase 2 onward,
   refresh the diagram so it shows the additional model `derivatives()`
   calls (GENROU, ST1A, PSS1A, TGOV1, IBR controllers) and how each
   model's outputs feed the others within the same timestep. The
   point of the diagram is to make the link
   *derivatives → trapezoidal rule → state changes* explicit and
   visual, so colleagues never have to read the integrator source to
   understand what is happening per step.

8. **Ship a PSSE benchmark for every phase.** Each phase delivers a
   `psse/phaseN/` folder containing (a) a `.dyr` dynamics file with
   the model parameters in PSSE syntax, (b) a psspy automation script
   (`run_phaseN_fault.py`) that drives the same disturbance the smib
   notebook simulates, and (c) a README with both GUI and psspy
   instructions plus a table of smib reference numbers (V from PF,
   $|E'|$, peak swing, mid-fault P/Q, CCT) so the user can cross-check
   independently. Expected agreement is within ~3 % on all canonical
   traces given matched integrator step size (smib $h = 2$ ms, set
   PSSE `DELT = 0.002`). Two purposes: (i) **correctness floor** — if
   smib and PSSE disagree on the same setup, one of them has a bug,
   and the disagreement points at exactly which model / init step /
   integrator detail to debug; (ii) **pedagogical confidence** —
   colleagues trust PSSE, so showing smib reproduces PSSE lets them
   trust smib too, and then use smib as a transparent way to
   understand what PSSE is doing. See `psse/phase1/` as the template.

These rules apply to every phase: GENCLS, GENROU + ST1A + PSS1A +
TGOV1, IBR, GFM, SVC/STATCOM/SynCon, scenario browser, small-signal.
Bake them in from the start of each new notebook rather than retrofit
them after review.

## PLL design (Phase 3 onward)

The PLL is the first control block that has to be right before any
grid-following IBR behaves. We use the Manitoba Hydro International
(MHI) PSCAD v5.0 SRF-PLL topology as the reference — see
`phase_locked_loop_pll.pdf` in the workspace folder.

MHI block diagram:

    Va,Vb,Vc -> Park(theta) -> Vd, Vq -> ATAN2(Vd, Vq) -> phi_err
             -> (subtract offset) -> PI(Kp, Ki/s) -> delta_omega
             -> (add omega0) -> frequency limiter -> 1/s -> theta

Three properties that matter and that many SRF-PLL implementations get
wrong:

- **Detector is `atan2(Vd, Vq)`**, not a normalised Vq. Linear over the
  full +/- pi range, so large-angle step tests converge. A sin-detector
  saturates near +/- pi/2 and loop gain collapses there.
- **Axis convention is swapped from the usual IEEE/PSSE form**: in MHI,
  Vq is cosine-aligned (top row of eq 1) and Vd is sine-aligned
  (bottom row). `atan2(Vd, Vq)` then returns `angle(V_terminal) - theta`
  directly. Module docstrings call this out to prevent sign-bug
  debugging days.
- **Base frequency is added downstream of the PI**, and a frequency
  limiter sits between the PI output and the integrator. The PI
  computes delta_omega, not absolute frequency, and the limiter is
  where conditional-integration anti-windup lives.

In our phasor simulator we do not carry abc signals. The Park + atan2
collapse to one line:

    phi_err = wrap(angle(V_terminal) - theta_pll - offset_angle)

The magnitude `|V|` cancels out, so the detector is inherently robust
at low voltage with no hand-rolled guard.

Two parameter presets are exposed (`smib.models.pll_preset(...)`):

| preset | origin | Kp (rad-input) | Ki (rad-input) | wn | zeta |
|---|---|---|---|---|---|
| `"bulk_AC"` | MHI Table 1 Model 1, deg->rad | 5.24 | 17.45 | 4.18 rad/s | 0.63 |
| `"ibr_fast"` | Textbook IBR vendor (wn=2*pi*30, zeta=1/sqrt(2)) | 266.6 | 35530 | 188 rad/s | 0.707 |

The SCR-sweep demo in Phase 3 switches between the two and shows the
weak-grid instability that disappears when the PLL is slowed down.

**Caveat (non-negotiable to flag in every PLL notebook).** MHI's
Example 3 shows a 120 Hz ripple on tracked frequency during a
single-line-to-ground fault, caused by negative-sequence content in
abc. Positive-sequence phasor simulation does not reproduce that
ripple. Our PLL traces will look cleaner than PSCAD EMT traces of the
same case; this is a limitation of RMS, not a modelling bug.
