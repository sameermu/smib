"""Disturbance scenarios.

Scenarios are plain functions that take the Network and simulator state and
schedule a change at a time. Keeping them as simple callables avoids an
over-engineered event system.

Supported for Phase 0..1:
- apply_three_phase_fault(network, t_start, t_clear, Z_fault)
- apply_voltage_step(network, t_start, V_new)
- apply_setpoint_step(target_model, attr, t_start, new_value)
"""
import math

import numpy as np


def three_phase_fault_schedule(t_start: float, t_clear: float, Z_fault: complex):
    """Returns a function(t, network) that sets/clears a bolted shunt fault.

    Z_fault = 0 gives a bolted fault. A small positive resistance can improve
    numerical conditioning if you see ill-conditioning during the fault.
    """
    Y_fault = 1.0 / Z_fault if Z_fault != 0 else 1e6 + 0j

    def apply(t_now, dt, network):
        if abs(t_now - t_start) < dt / 2:
            network.set_fault(Y_fault)
        elif abs(t_now - t_clear) < dt / 2:
            network.clear_fault()
    return apply


def voltage_step_schedule(t_start: float, V_new: complex):
    """Step the infinite bus voltage at t_start."""
    def apply(t_now, dt, network):
        if abs(t_now - t_start) < dt / 2:
            network.set_slack_voltage(V_new)
    return apply


def setpoint_step_schedule(target_model, attr: str, t_start: float, new_value: float):
    """Step a setpoint (e.g. Vref, Pref) on a model."""
    def apply(t_now, dt, network):
        if abs(t_now - t_start) < dt / 2:
            target_model.inputs[attr] = new_value
    return apply


def grid_frequency_ramp_schedule(t_start: float, delta_f_pu: float,
                                 ramp_time: float = 1.0, f0: float = 50.0):
    """Emulate a system frequency event by ramping the infinite-bus angle.

    Physics: our network phasors live in a frame rotating at exactly
    omega_0.  If the wider grid's frequency deviates by Δf, the infinite
    bus phasor rotates *relative to that frame* at dθ/dt = 2π·f0·Δω̄
    rad/s — i.e. the slack angle ramps.  A machine synchronised to that
    bus must, in the new steady state, rotate at the grid's speed, so
    its slip settles at ω̄ = Δω̄_grid and the governor sees a sustained
    speed deviation:  ΔPm = -ω̄/R.  This is the SMIB stand-in for
    "the rest of the system lost generation / gained load".

    Parameters
    ----------
    t_start    : event start [s].
    delta_f_pu : final grid frequency deviation, pu of f0.  Negative =
                 under-frequency (the common, load-gained case).
                 E.g. -0.004 pu on 50 Hz = -0.2 Hz.
    ramp_time  : seconds over which the deviation ramps linearly from
                 0 to delta_f_pu (instant steps in frequency are not
                 physical — system frequency moves on inertia).
    f0         : system base frequency [Hz] (50 Hz throughout smib).

    Implementation: each step we advance the slack phasor's CURRENT
    angle by this step's extra rotation, dθ = Δω̄ · ω0 · dt.  During the
    ramp Δω̄ grows linearly; after the ramp it is constant, so the slack
    angle ramps linearly forever (the machine's delta tracks it — plot
    delta relative to the slack angle if you want a settling trace).
    Because each schedule increments the live network angle rather than
    keeping a private copy, multiple frequency events compose: schedule
    a -Δf ramp and a later +Δf ramp to model an event plus recovery.
    """
    w0 = 2.0 * math.pi * f0

    def apply(t_now, dt, network):
        if t_now < t_start:
            return
        # Current per-unit speed offset contributed by this event.
        if t_now <= t_start + ramp_time:
            dw = delta_f_pu * (t_now - t_start) / ramp_time
        else:
            dw = delta_f_pu
        # Advance the slack angle by this step's extra rotation.
        theta = float(np.angle(network.V_slack))
        Vmag = abs(network.V_slack)
        network.set_slack_voltage(Vmag * np.exp(1j * (theta + dw * w0 * dt)))
    return apply
