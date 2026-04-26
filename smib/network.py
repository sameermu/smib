"""Ybus and algebraic network solve.

For SMIB the network is trivial: one branch between the generator bus and
the infinite bus. We still build Ybus generally so later upgrades to N-bus
don't require rewriting. Shunt faults are applied by adding a shunt admittance
at the faulted bus during the fault interval.
"""
import numpy as np

from .powerflow import ybus_2bus


class Network:
    """Two-bus network for SMIB.

    Bus 0: generator bus (index 0 in Ybus, but labelled bus 1 in docs).
    Bus 1: infinite bus.
    """

    def __init__(self, R, X, V_slack_mag=1.0, V_slack_ang=0.0):
        self.R = R
        self.X = X
        self.V_slack = V_slack_mag * np.exp(1j * V_slack_ang)
        self._base_Ybus = ybus_2bus(R, X)
        self._shunt = 0.0 + 0.0j  # fault shunt at bus 0

    def set_fault(self, Y_fault):
        """Apply a shunt admittance at bus 0 (generator bus)."""
        self._shunt = Y_fault

    def clear_fault(self):
        self._shunt = 0.0 + 0.0j

    def set_slack_voltage(self, V):
        """Change the infinite bus voltage (for V step tests)."""
        self.V_slack = V

    def ybus(self):
        Y = self._base_Ybus.copy()
        Y[0, 0] += self._shunt
        return Y

    def solve(self, I_inj_gen):
        """Given current injection at bus 0, solve for V at bus 0.

        From Ybus: [V0; V_slack] -> [I0; I_slack]
          I0 = Y00 V0 + Y01 V_slack
        =>  V0 = (I_inj_gen - Y01 V_slack) / Y00
        """
        Y = self.ybus()
        V0 = (I_inj_gen - Y[0, 1] * self.V_slack) / Y[0, 0]
        return V0
