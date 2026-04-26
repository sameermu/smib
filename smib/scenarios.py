"""Disturbance scenarios.

Scenarios are plain functions that take the Network and simulator state and
schedule a change at a time. Keeping them as simple callables avoids an
over-engineered event system.

Supported for Phase 0..1:
- apply_three_phase_fault(network, t_start, t_clear, Z_fault)
- apply_voltage_step(network, t_start, V_new)
- apply_setpoint_step(target_model, attr, t_start, new_value)
"""
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
