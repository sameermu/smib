"""Two-bus Newton-Raphson power flow for SMIB.

Bus 1: generator terminal (PV or PQ).
Bus 2: infinite bus (slack). |V|, angle fixed.

For an N-bus generalisation later, this file is the canonical reference
for how we apply Newton-Raphson to the power flow equations.

Power flow equations at bus k (polar form):
    P_k = sum_j |V_k||V_j| ( G_kj cos(th_k - th_j) + B_kj sin(th_k - th_j) )
    Q_k = sum_j |V_k||V_j| ( G_kj sin(th_k - th_j) - B_kj cos(th_k - th_j) )

Unknowns: (theta_1, |V_1|) if bus 1 is PQ; (theta_1) if PV.
"""
import numpy as np


def ybus_2bus(R, X):
    """Ybus for a single RL branch between bus 1 (gen) and bus 2 (inf bus)."""
    z = complex(R, X)
    y = 1.0 / z
    return np.array([[y, -y], [-y, y]], dtype=complex)


def _mismatch(V1_mag, th1, V2_mag, th2, Ybus, P_spec, Q_spec, bus_type):
    """Compute active and (if PQ) reactive power mismatches at bus 1."""
    G = Ybus.real
    B = Ybus.imag
    # bus 1 injection
    P1 = V1_mag * (
        V1_mag * (G[0, 0] * np.cos(0.0) + B[0, 0] * np.sin(0.0))
        + V2_mag * (G[0, 1] * np.cos(th1 - th2) + B[0, 1] * np.sin(th1 - th2))
    )
    Q1 = V1_mag * (
        V1_mag * (G[0, 0] * np.sin(0.0) - B[0, 0] * np.cos(0.0))
        + V2_mag * (G[0, 1] * np.sin(th1 - th2) - B[0, 1] * np.cos(th1 - th2))
    )
    dP = P_spec - P1
    dQ = Q_spec - Q1 if bus_type == "PQ" else 0.0
    return dP, dQ, P1, Q1


def two_bus_pf(P_spec, Q_spec, V_slack_mag, V_slack_ang, R, X,
               bus_type="PQ", V1_guess=1.0, th1_guess=0.0,
               tol=1e-10, max_iter=30, verbose=False):
    """Solve the 2-bus power flow by Newton-Raphson.

    Parameters
    ----------
    P_spec : float
        Active power injection at bus 1 (pu, generator convention: positive out).
    Q_spec : float
        Reactive power at bus 1 (ignored if bus_type == "PV").
    V_slack_mag, V_slack_ang : float
        Infinite bus magnitude (pu) and angle (rad).
    R, X : float
        Branch impedance in pu.
    bus_type : {"PQ", "PV"}
        "PQ" solves for (theta_1, V1). "PV" fixes |V_1| at V1_guess and solves theta_1.
    verbose : bool
        If True, prints a table of each Newton-Raphson iteration showing
        the current state, mismatch, and update.  Useful for notebooks
        and pedagogy; production callers should leave this False.

    Returns
    -------
    V1 : complex
        Bus 1 voltage phasor in pu.
    iters : int
        Number of Newton iterations.
    """
    Y = ybus_2bus(R, X)
    V1 = V1_guess
    th1 = th1_guess
    V2 = V_slack_mag
    th2 = V_slack_ang

    if verbose:
        if bus_type == "PQ":
            print(f"{'iter':>4}  {'|V1|':>9}  {'angle deg':>10}  "
                  f"{'mis P':>11}  {'mis Q':>11}  {'|mismatch|':>12}")
        else:
            print(f"{'iter':>4}  {'|V1|':>9}  {'angle deg':>10}  "
                  f"{'mis P':>11}  {'|mismatch|':>12}")

    for k in range(max_iter):
        dP, dQ, _, _ = _mismatch(V1, th1, V2, th2, Y, P_spec, Q_spec, bus_type)
        if bus_type == "PQ":
            mismatch = np.array([dP, dQ])
        else:
            mismatch = np.array([dP])

        max_mis = float(np.max(np.abs(mismatch)))

        if verbose:
            ang_deg = np.degrees(th1)
            if bus_type == "PQ":
                print(f"{k:>4}  {V1:>9.6f}  {ang_deg:>10.4f}  "
                      f"{dP:>11.3e}  {dQ:>11.3e}  {max_mis:>12.3e}")
            else:
                print(f"{k:>4}  {V1:>9.6f}  {ang_deg:>10.4f}  "
                      f"{dP:>11.3e}  {max_mis:>12.3e}")

        if max_mis < tol:
            if verbose:
                print(f">>> Converged in {k} iterations to |mismatch| < {tol:.0e}")
            return V1 * np.exp(1j * th1), k

        # Numerical Jacobian keeps this code readable. For a 2x2 case the
        # analytical Jacobian is short, but numerical is easier to inspect.
        eps = 1e-8
        J = np.zeros((len(mismatch), len(mismatch)))

        dP1, dQ1, _, _ = _mismatch(V1, th1 + eps, V2, th2, Y, P_spec, Q_spec, bus_type)
        J[0, 0] = -(dP1 - dP) / eps
        if bus_type == "PQ":
            J[1, 0] = -(dQ1 - dQ) / eps
            dP2, dQ2, _, _ = _mismatch(V1 + eps, th1, V2, th2, Y, P_spec, Q_spec, bus_type)
            J[0, 1] = -(dP2 - dP) / eps
            J[1, 1] = -(dQ2 - dQ) / eps

        dx = np.linalg.solve(J, mismatch)
        th1 += dx[0]
        if bus_type == "PQ":
            V1 += dx[1]

    raise RuntimeError(f"Power flow did not converge in {max_iter} iterations")
