"""Dynamic models.

Each model inherits from Model (see base.py) and implements four methods:
    initialise(V, I, **interface_inputs)
    derivatives(inputs) -> dict
    current_injection(V) -> complex
    algebraic_output() -> dict

Add new models here one at a time. Each model must pass the correctness
floor (see tests/) before being committed.
"""
from .base import Model
from .gencls import GENCLS
from .genrou import GENROU
from .pll import PLL, preset as pll_preset
from .pss1a import PSS1A
from .st1a import ST1A

__all__ = ["Model", "GENCLS", "GENROU", "PLL", "pll_preset", "PSS1A", "ST1A"]
