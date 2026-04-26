"""smib — a transparent RMS/phasor-domain SMIB simulator for teaching.

Not optimised for speed. Optimised for "I can read the code and understand
exactly what is happening."

Conventions
-----------
- pu throughout, machine base unless stated.
- Angles in radians. Frequency in pu of ws (= 2*pi*50 or 2*pi*60).
- Global synchronous reference frame is DQ (capital). Each synchronous
  machine rotates DQ into its own dq via its rotor angle delta.
- State vectors are named dicts per model, flattened only at the integrator
  boundary. Reader sees `genrou.Eqpp`, not `x[14]`.
"""
__version__ = "0.0.1"
