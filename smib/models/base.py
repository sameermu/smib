"""Abstract base class for dynamic models.

Design choice: state and params are named dicts, not numpy arrays. A dict
costs us a bit of speed, but it means every plot label, every debugger
inspection, and every notebook printout reads `genrou.state["Eqpp"]`, not
`x[14]`. For a teaching tool that is the correct trade.

At the integrator boundary we `flatten()` dict -> np.ndarray and reverse
with `unflatten()`. That is the only place the dict/array duality leaks.
"""
from __future__ import annotations

import numpy as np


class Model:
    name: str = ""
    # Ordered list of state variable names. Subclasses MUST set this.
    state_keys: tuple = ()

    def __init__(self, name: str, params: dict):
        self.name = name
        self.params = dict(params)
        self.state: dict = {k: 0.0 for k in self.state_keys}
        # Inputs filled by the simulator each step (terminal V, I, setpoints, etc.)
        self.inputs: dict = {}

    # ----- required interface -------------------------------------------
    def initialise(self, V: complex, S: complex, **kwargs) -> None:
        """Back-calculate internal states so that dx/dt = 0 at t=0.

        V : complex terminal voltage phasor (pu).
        S = P + jQ : complex power injection at terminal (pu, generator convention).
        """
        raise NotImplementedError

    def derivatives(self) -> dict:
        """Return dx/dt as a dict keyed by state_keys. Uses self.state and self.inputs."""
        raise NotImplementedError

    def current_injection(self, V: complex) -> complex:
        """Return the model's current injection into its bus (pu)."""
        raise NotImplementedError

    def algebraic_output(self) -> dict:
        """Optional outputs used by other models (e.g. Efd from AVR to machine)."""
        return {}

    # ----- helpers ------------------------------------------------------
    def flatten(self) -> np.ndarray:
        return np.array([self.state[k] for k in self.state_keys], dtype=float)

    def unflatten(self, x: np.ndarray) -> None:
        for i, k in enumerate(self.state_keys):
            self.state[k] = float(x[i])

    def assert_initialised(self, tol: float = 1e-6) -> None:
        """Sanity check: dx/dt should be ~0 at t=0 if initialisation is consistent."""
        d = self.derivatives()
        worst = max(abs(v) for v in d.values()) if d else 0.0
        if worst > tol:
            offenders = {k: v for k, v in d.items() if abs(v) > tol}
            raise AssertionError(
                f"{self.name}: initialisation inconsistent, dx/dt not zero: {offenders}"
            )
