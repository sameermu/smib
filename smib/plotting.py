"""Plotly traces + ipywidgets sliders for embedded notebook scenarios.

Why Plotly over matplotlib: every phase notebook serves as both the
report and the scenario browser, so hover/zoom on fault traces is worth
the dependency.  Matplotlib remains a fallback option for static PDF
export if/when needed.

Why a `scenario_slider` helper: the alternative is for every notebook
to re-write the same `interact(...)` boilerplate with `clear_output`,
figure caching, and trace updating.  This helper hides that so the
notebook reads as physics, not plumbing.
"""
from __future__ import annotations

from typing import Callable, Iterable, Mapping, Sequence

import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "plotly is required: pip install plotly"
    ) from e


# ---------- canonical trace order --------------------------------------

CANONICAL = ("P", "Q", "|V|", "Id", "Iq")
"""Five traces that every scenario must plot, in this order, regardless
of the technology under test.  See README "Standard plots" section."""


# ---------- core plotting ----------------------------------------------

def plot_traces(t: np.ndarray, traces: Mapping[str, np.ndarray],
                title: str = "", height_per_row: int = 140) -> "go.Figure":
    """One stacked subplot per trace, shared x-axis.

    `traces` is an ordered dict {label: 1d_array}.  Use Python 3.7+
    insertion order to control plot order; for the canonical 5,
    construct the dict as ``{k: arrays[k] for k in CANONICAL}``.
    """
    labels = list(traces.keys())
    n = len(labels)
    fig = make_subplots(rows=n, cols=1, shared_xaxes=True, vertical_spacing=0.025)
    for i, lbl in enumerate(labels, start=1):
        fig.add_trace(
            go.Scatter(x=t, y=traces[lbl], mode="lines", name=lbl, showlegend=False),
            row=i, col=1,
        )
        fig.update_yaxes(title_text=lbl, row=i, col=1)
    fig.update_xaxes(title_text="time [s]", row=n, col=1)
    fig.update_layout(
        title=title or None,
        height=height_per_row * n + 60,
        margin=dict(l=60, r=20, t=40 if title else 10, b=40),
        template="plotly_white",
    )
    return fig


def shade_event_window(fig: "go.Figure", t_start: float, t_end: float,
                       label: str = "fault",
                       color: str = "rgba(220,80,80,0.12)") -> None:
    """Add a coloured vertical band across all subplots to mark a fault
    or LV/HV event.  Mutates `fig` in place.
    """
    fig.add_vrect(
        x0=t_start, x1=t_end,
        fillcolor=color, line_width=0,
        annotation_text=label, annotation_position="top left",
    )


# ---------- ride-through validator -------------------------------------

def validate_ride_through(t: np.ndarray, V: np.ndarray, Q: np.ndarray, Iq: np.ndarray,
                          event_window: tuple[float, float], kind: str = "LV",
                          tol: float = 1e-3) -> dict:
    """Hard-rule check: during LV (V < 1.0 in window) Q and Iq must rise;
    during HV they must fall.  Returns a dict with PASS/FAIL and slope
    diagnostics.  Notebooks call this after the run and print the dict.

    The slope test is intentionally cheap (linear fit) and meant as a
    canary, not a substitute for inspecting the plot.

    See README "Ride-through validation rule" section.
    """
    t = np.asarray(t)
    mask = (t >= event_window[0]) & (t <= event_window[1])
    if not mask.any():
        return {"status": "SKIP", "reason": "event window outside time array"}

    Q_w, Iq_w, t_w = np.asarray(Q)[mask], np.asarray(Iq)[mask], t[mask]
    slope_Q = float(np.polyfit(t_w, Q_w, 1)[0])
    slope_Iq = float(np.polyfit(t_w, Iq_w, 1)[0])

    if kind == "LV":
        ok = (slope_Q > tol) and (slope_Iq > tol)
        sign = "+"
    elif kind == "HV":
        ok = (slope_Q < -tol) and (slope_Iq < -tol)
        sign = "-"
    else:
        raise ValueError(f"kind must be 'LV' or 'HV', got {kind!r}")

    return {
        "status": "PASS" if ok else "FAIL",
        "kind": kind,
        "expected_sign": sign,
        "slope_Q": slope_Q,
        "slope_Iq": slope_Iq,
        "window": event_window,
    }


# ---------- single-line diagram (SLD) ---------------------------------

def plot_sld(V_gen: complex, S_gen: complex,
             V_inf: complex = 1.0 + 0j,
             R_line: float = 0.0, X_line: float = 0.5,
             gen_label: str = "GEN",
             title: str = "SMIB single-line diagram (post-LF)") -> "go.Figure":
    """Render a PSSE-style single-line diagram of the post-LF state.

    Shows two buses, the line connecting them, the generator on bus 1
    and the infinite-bus / slack symbol on bus 2.  Each bus is labelled
    with its converged voltage magnitude and angle.  The line is
    labelled with its impedance.  The active and reactive flow on the
    line is annotated alongside the impedance, with an arrow indicating
    direction (gen-bus → infinite-bus is positive when the machine
    exports power).

    This deliberately mimics the look of a PSSE PSSPLT diagram: blocky
    bus rectangles, simple line, P/Q labelled with arrows, slack symbol
    distinct from generator symbol.

    Parameters
    ----------
    V_gen : complex
        Converged voltage phasor at the generator bus, pu.
    S_gen : complex
        Generator power injection, pu (P + jQ, generator convention).
    V_inf : complex
        Infinite-bus voltage phasor, pu.  Default 1.0 / 0 deg.
    R_line, X_line : float
        Branch impedance, pu.
    gen_label : str
        Label for the generator (e.g. "GENCLS", "GENROU", "REGC_A").
    """
    # Geometry — keep coordinates simple.
    bus1_x, bus1_y = 0.20, 0.50
    bus2_x, bus2_y = 0.80, 0.50
    bus_w, bus_h = 0.10, 0.18

    # Line flow: I = (V_gen - V_inf) / Z_line; S_line at sending end.
    Z = complex(R_line, X_line)
    I_line = (V_gen - V_inf) / Z if abs(Z) > 0 else 0j
    S_line_send = V_gen * np.conj(I_line)
    P_line = float(S_line_send.real)
    Q_line = float(S_line_send.imag)

    fig = go.Figure()

    # ---- Branch (transmission line) ----
    fig.add_shape(
        type="line",
        x0=bus1_x + bus_w / 2, y0=bus1_y,
        x1=bus2_x - bus_w / 2, y1=bus2_y,
        line=dict(color="#222", width=3),
    )

    # Line impedance + power flow label, mid-span.
    arrow = "→" if P_line >= 0 else "←"
    line_label = (
        f"<b>Z = {R_line:.3f} + j{X_line:.3f} pu</b><br>"
        f"P = {abs(P_line):.3f}  {arrow}<br>"
        f"Q = {abs(Q_line):.3f}  {arrow}"
    )
    fig.add_annotation(
        x=(bus1_x + bus2_x) / 2, y=bus1_y + 0.10,
        text=line_label, showarrow=False,
        font=dict(size=12, family="Courier New"),
        align="center",
        bgcolor="rgba(255,255,255,0.85)",
        bordercolor="#888", borderwidth=1, borderpad=4,
    )

    # ---- Bus 1 (generator bus) ----
    _draw_bus(
        fig, bus1_x, bus1_y, bus_w, bus_h,
        title="Bus 1 (gen)",
        V=V_gen,
        fill="#FFEFE0",
    )

    # ---- Bus 2 (infinite bus / slack) ----
    _draw_bus(
        fig, bus2_x, bus2_y, bus_w, bus_h,
        title="Bus 2 (slack)",
        V=V_inf,
        fill="#E0F0FF",
    )

    # ---- Generator symbol below bus 1 ----
    gen_cx, gen_cy = bus1_x, bus1_y - 0.28
    gen_r = 0.05
    fig.add_shape(
        type="circle",
        x0=gen_cx - gen_r, y0=gen_cy - gen_r,
        x1=gen_cx + gen_r, y1=gen_cy + gen_r,
        line=dict(color="#222", width=2),
        fillcolor="#FFFFFF",
    )
    fig.add_annotation(
        x=gen_cx, y=gen_cy, text="<b>~</b>",
        showarrow=False, font=dict(size=20),
    )
    # Lead from bus down to gen circle.
    fig.add_shape(
        type="line",
        x0=gen_cx, y0=bus1_y - bus_h / 2,
        x1=gen_cx, y1=gen_cy + gen_r,
        line=dict(color="#222", width=2),
    )
    # Gen label and injection.
    fig.add_annotation(
        x=gen_cx + 0.08, y=gen_cy,
        text=(f"<b>{gen_label}</b><br>"
              f"P = {S_gen.real:+.3f} pu<br>"
              f"Q = {S_gen.imag:+.3f} pu"),
        showarrow=False, font=dict(size=11, family="Courier New"),
        align="left", xanchor="left",
    )

    # ---- Slack symbol below bus 2 ----
    slk_x, slk_y = bus2_x, bus2_y - 0.28
    fig.add_shape(
        type="line",
        x0=slk_x, y0=bus2_y - bus_h / 2,
        x1=slk_x, y1=slk_y + 0.05,
        line=dict(color="#222", width=2),
    )
    # PSSE slack symbol: triangle pointing down with a horizontal bar
    fig.add_shape(
        type="path",
        path=f"M {slk_x-0.04},{slk_y+0.05} L {slk_x+0.04},{slk_y+0.05} L {slk_x},{slk_y-0.02} Z",
        line=dict(color="#222", width=2),
        fillcolor="#FFFFFF",
    )
    fig.add_annotation(
        x=slk_x + 0.08, y=slk_y,
        text=("<b>SLACK</b><br>"
              f"|V| = {abs(V_inf):.3f} pu<br>"
              f"∠   = {np.degrees(np.angle(V_inf)):+.2f}°"),
        showarrow=False, font=dict(size=11, family="Courier New"),
        align="left", xanchor="left",
    )

    # ---- Layout ----
    fig.update_xaxes(visible=False, range=[0, 1])
    fig.update_yaxes(visible=False, range=[0, 1], scaleanchor="x", scaleratio=1)
    fig.update_layout(
        title=title, height=420,
        margin=dict(l=20, r=20, t=50, b=20),
        plot_bgcolor="white",
        showlegend=False,
    )
    return fig


def _draw_bus(fig: "go.Figure", cx: float, cy: float, w: float, h: float,
              title: str, V: complex, fill: str) -> None:
    """Draw one bus rectangle with V magnitude and angle inside."""
    fig.add_shape(
        type="rect",
        x0=cx - w / 2, y0=cy - h / 2,
        x1=cx + w / 2, y1=cy + h / 2,
        line=dict(color="#222", width=2),
        fillcolor=fill,
    )
    fig.add_annotation(
        x=cx, y=cy + h / 2 - 0.025,
        text=f"<b>{title}</b>", showarrow=False,
        font=dict(size=11), yanchor="top",
    )
    fig.add_annotation(
        x=cx, y=cy,
        text=(f"|V| = {abs(V):.4f} pu<br>"
              f"∠   = {np.degrees(np.angle(V)):+.3f}°"),
        showarrow=False,
        font=dict(size=11, family="Courier New"),
        align="center",
    )


# ---------- ipywidgets slider helper -----------------------------------

def scenario_slider(run_fn: Callable[..., tuple],
                    sliders: Sequence[Mapping],
                    title: str = "",
                    extra_traces: Iterable[str] = ()) -> "object":
    """Wire an ipywidgets slider panel to a Plotly figure.

    `run_fn(**kwargs)` must return `(t, traces_dict)` where traces_dict
    contains at minimum the 5 canonical traces.  `sliders` is a list of
    dicts like:

        {"name": "H", "min": 1, "max": 10, "step": 0.5, "value": 4.0,
         "description": "Inertia H [s]"}

    Returns the assembled VBox widget.  Designed so a notebook cell
    boils down to:

        scenario_slider(my_run, sliders=[...], title="GENCLS swing")

    ipywidgets is imported lazily so plotting.py is importable in a
    plain Python REPL without the widgets dependency.
    """
    try:
        import ipywidgets as widgets
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "ipywidgets is required for scenario_slider: pip install ipywidgets"
        ) from e

    slider_widgets = []
    for s in sliders:
        sl = widgets.FloatSlider(
            value=s["value"], min=s["min"], max=s["max"], step=s.get("step", 0.1),
            description=s.get("description", s["name"]),
            continuous_update=False,
            style={"description_width": "initial"},
            layout=widgets.Layout(width="500px"),
        )
        sl._slider_name = s["name"]  # noqa: SLF001
        slider_widgets.append(sl)

    out = widgets.Output()

    def _replot(*_args):
        kwargs = {sl._slider_name: sl.value for sl in slider_widgets}  # noqa: SLF001
        t, traces = run_fn(**kwargs)
        ordered = {k: traces[k] for k in CANONICAL if k in traces}
        for k in extra_traces:
            if k in traces:
                ordered[k] = traces[k]
        fig = plot_traces(t, ordered, title=title)
        with out:
            out.clear_output(wait=True)
            fig.show()

    for sl in slider_widgets:
        sl.observe(_replot, names="value")

    panel = widgets.VBox([*slider_widgets, out])
    _replot()
    return panel
