"""The flat-line test is the non-negotiable correctness floor.

Procedure: run any SMIB case with no disturbance for 10 seconds. Every
state variable in every model must stay flat to within 1e-5 pu. If it
drifts, the initialisation is wrong — stop and fix before anything else.

This file will grow as we add models. For Phase 0 it only exercises the
power flow + network + integrator plumbing (no dynamic model yet).
"""
import numpy as np

from smib.network import Network
from smib.powerflow import two_bus_pf


def test_powerflow_converges_simple_case():
    V1, iters = two_bus_pf(
        P_spec=0.5, Q_spec=0.1,
        V_slack_mag=1.0, V_slack_ang=0.0,
        R=0.0, X=0.2,
        bus_type="PQ",
    )
    assert iters < 10
    assert 0.9 < abs(V1) < 1.1


def test_powerflow_pv_bus():
    V1, iters = two_bus_pf(
        P_spec=0.5, Q_spec=0.0,  # Q_spec ignored
        V_slack_mag=1.0, V_slack_ang=0.0,
        R=0.0, X=0.2,
        bus_type="PV",
        V1_guess=1.02,
    )
    assert iters < 10
    assert abs(abs(V1) - 1.02) < 1e-8


def test_network_round_trip():
    """If we inject the current consistent with the PF solution, we should
    recover the same terminal voltage."""
    V1, _ = two_bus_pf(0.5, 0.1, 1.0, 0.0, 0.0, 0.2, bus_type="PQ")
    S = complex(0.5, 0.1)
    I_inj = np.conj(S / V1)
    net = Network(R=0.0, X=0.2)
    V1_solved = net.solve(I_inj)
    assert abs(V1_solved - V1) < 1e-8
