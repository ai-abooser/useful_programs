"""Generate plots for the oblique projectile motion report.

Run from any directory with:
    python3 pages/solutions/04_sol_vscode/generuj_wykresy_rzut_ukosny.py

Requires matplotlib. The numerical integration uses only the Python standard
library and a fixed-step fourth-order Runge-Kutta method.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
OUT = HERE / "wykresy"
OUT.mkdir(parents=True, exist_ok=True)

V0 = 36.0
ALPHA_DEG = 64.0
ALPHA = math.radians(ALPHA_DEG)
G = 9.81
MASS = 0.145
DIAMETER = 0.073
RHO_AIR = 1.225
CD = 0.47
AREA = math.pi * DIAMETER**2 / 4.0
DRAG_FACTOR = 0.5 * RHO_AIR * CD * AREA
DT = 0.002
MAX_TIME = 20.0

VX0 = V0 * math.cos(ALPHA)
VY0 = V0 * math.sin(ALPHA)
IDEAL_FLIGHT_TIME = 2.0 * VY0 / G
IDEAL_RANGE = VX0 * IDEAL_FLIGHT_TIME
IDEAL_MAX_HEIGHT = VY0**2 / (2.0 * G)

BLUE = "#2463A6"
ORANGE = "#D47722"
GREEN = "#168578"
PURPLE = "#8455A4"
INK = "#233044"

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "lines.linewidth": 2.2,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    }
)


def acceleration(vx: float, vy: float) -> tuple[float, float]:
    """Return acceleration components for quadratic drag and gravity."""
    speed = math.hypot(vx, vy)
    drag_per_mass = DRAG_FACTOR * speed / MASS
    return -drag_per_mass * vx, -G - drag_per_mass * vy


def derivative(state: tuple[float, float, float, float]) -> tuple[float, ...]:
    x, y, vx, vy = state
    ax, ay = acceleration(vx, vy)
    return vx, vy, ax, ay


def rk4_step(
    state: tuple[float, float, float, float], step: float
) -> tuple[float, float, float, float]:
    """Advance the four-component state by one RK4 time step."""
    k1 = derivative(state)
    state2 = tuple(s + 0.5 * step * k for s, k in zip(state, k1))
    k2 = derivative(state2)
    state3 = tuple(s + 0.5 * step * k for s, k in zip(state, k2))
    k3 = derivative(state3)
    state4 = tuple(s + step * k for s, k in zip(state, k3))
    k4 = derivative(state4)
    return tuple(
        s + step * (a + 2.0 * b + 2.0 * c + d) / 6.0
        for s, a, b, c, d in zip(state, k1, k2, k3, k4)
    )


def simulate_with_drag() -> tuple[list[float], list[tuple[float, ...]]]:
    """Integrate until the projectile returns to ground level."""
    times = [0.0]
    states = [(0.0, 0.0, VX0, VY0)]
    time = 0.0
    state = states[0]

    while time < MAX_TIME:
        step = min(DT, MAX_TIME - time)
        next_state = rk4_step(state, step)
        next_time = time + step

        if next_state[1] < 0.0:
            fraction = state[1] / (state[1] - next_state[1])
            landing_time = time + fraction * step
            landing_state = tuple(
                start + fraction * (end - start)
                for start, end in zip(state, next_state)
            )
            landing_state = (*landing_state[:1], 0.0, *landing_state[2:])
            times.append(landing_time)
            states.append(landing_state)
            return times, states

        times.append(next_time)
        states.append(next_state)
        time, state = next_time, next_state

    raise RuntimeError(f"Projectile did not reach the ground within {MAX_TIME} s")


def save_figure(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


drag_times, drag_states = simulate_with_drag()
drag_x = [state[0] for state in drag_states]
drag_y = [state[1] for state in drag_states]
drag_vx = [state[2] for state in drag_states]
drag_vy = [state[3] for state in drag_states]
drag_speeds = [math.hypot(vx, vy) for vx, vy in zip(drag_vx, drag_vy)]
drag_ax, drag_ay = zip(*(acceleration(vx, vy) for vx, vy in zip(drag_vx, drag_vy)))
drag_accelerations = [math.hypot(ax, ay) for ax, ay in zip(drag_ax, drag_ay)]

ideal_count = 1501
ideal_times = [
    IDEAL_FLIGHT_TIME * i / (ideal_count - 1) for i in range(ideal_count)
]
ideal_x = [VX0 * time for time in ideal_times]
ideal_y = [VY0 * time - 0.5 * G * time**2 for time in ideal_times]

# Figure 1: the two trajectories.
fig, ax = plt.subplots(figsize=(9.2, 5.2))
ax.plot(ideal_x, ideal_y, color=BLUE, label="Bez oporu")
ax.plot(drag_x, drag_y, color=ORANGE, label="Z oporem kwadratowym")
ax.scatter([0.0], [0.0], color=INK, zorder=4, label="Punkt startu")
ax.scatter(
    [IDEAL_RANGE, drag_x[-1]],
    [0.0, 0.0],
    color=[BLUE, ORANGE],
    marker="x",
    s=70,
    zorder=4,
)
ax.set(
    xlabel="Położenie poziome x [m]",
    ylabel="Wysokość y [m]",
    title="Tor ruchu: model idealny i opór powietrza",
)
ax.legend()
save_figure(fig, "01_tor_ruchu")

# Figure 2: position components over time.
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True)
axes[0].plot(ideal_times, ideal_x, color=BLUE, label="Bez oporu")
axes[0].plot(drag_times, drag_x, color=ORANGE, label="Z oporem")
axes[0].set_ylabel("x [m]")
axes[0].set_title("Położenie poziome")
axes[0].legend()
axes[1].plot(ideal_times, ideal_y, color=BLUE, label="Bez oporu")
axes[1].plot(drag_times, drag_y, color=ORANGE, label="Z oporem")
axes[1].set(xlabel="Czas t [s]", ylabel="y [m]", title="Wysokość")
save_figure(fig, "02_polozenie_czas")

# Figure 3: velocity components and speed.
fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.4), sharex=True)
axes[0].plot(ideal_times, [VX0] * len(ideal_times), color=BLUE, label="vx, bez oporu")
axes[0].plot(ideal_times, [VY0 - G * t for t in ideal_times],
             color=GREEN, label="vy, bez oporu")
axes[0].plot(drag_times, drag_vx, color=ORANGE, ls="--", label="vx, z oporem")
axes[0].plot(drag_times, drag_vy, color=PURPLE, ls="--", label="vy, z oporem")
axes[0].set_ylabel("Składowa prędkości [m/s]")
axes[0].set_title("Składowe prędkości")
axes[0].legend(ncol=2, fontsize=8)
ideal_speeds = [math.hypot(VX0, VY0 - G * t) for t in ideal_times]
axes[1].plot(ideal_times, ideal_speeds, color=BLUE, label="Bez oporu")
axes[1].plot(drag_times, drag_speeds, color=ORANGE, label="Z oporem")
axes[1].set(xlabel="Czas t [s]", ylabel="Szybkość v [m/s]", title="Szybkość")
axes[1].legend()
save_figure(fig, "03_predkosc_czas")

# Figure 4: acceleration under gravity and drag.
fig, ax = plt.subplots(figsize=(9.2, 4.8))
ax.plot(drag_times, drag_ax, color=BLUE, label="ax, z oporem")
ax.plot(drag_times, drag_ay, color=ORANGE, label="ay, z oporem")
ax.plot(ideal_times, [0.0] * len(ideal_times), color=BLUE, ls=":",
        label="ax, bez oporu")
ax.plot(ideal_times, [-G] * len(ideal_times), color=ORANGE, ls=":",
        label="ay, bez oporu")
ax.set(xlabel="Czas t [s]", ylabel="Przyspieszenie [m/s²]",
       title="Składowe przyspieszenia")
ax.legend(ncol=2, fontsize=8)
save_figure(fig, "04_przyspieszenie_czas")

# Figure 5: mechanical energy and work done by drag.
initial_energy = 0.5 * MASS * V0**2
mechanical_energy = [
    0.5 * MASS * (vx**2 + vy**2) + MASS * G * y
    for _, y, vx, vy in drag_states
]
drag_power = [
    -DRAG_FACTOR * speed**3 for speed in drag_speeds
]
cumulative_drag_work = [0.0]
for i in range(1, len(drag_times)):
    elapsed = drag_times[i] - drag_times[i - 1]
    average_power = 0.5 * (drag_power[i] + drag_power[i - 1])
    cumulative_drag_work.append(cumulative_drag_work[-1] + elapsed * average_power)
energy_balance_error = max(
    abs(energy - (initial_energy + work))
    for energy, work in zip(mechanical_energy, cumulative_drag_work)
)
if energy_balance_error > 1e-3:
    raise RuntimeError(
        f"Energy balance check failed: maximum error {energy_balance_error:.3e} J"
    )

fig, ax = plt.subplots(figsize=(9.2, 4.8))
ax.plot(drag_times, mechanical_energy, color=BLUE, label="Energia mechaniczna")
ax.plot(
    drag_times,
    [initial_energy + work for work in cumulative_drag_work],
    color=ORANGE,
    ls="--",
    label="E₀ + praca oporu",
)
ax.axhline(initial_energy, color=GREEN, ls=":", label="E₀ (wartość początkowa)")
ax.set(
    xlabel="Czas t [s]",
    ylabel="Energia [J]",
    title="Spadek energii mechanicznej wskutek oporu",
)
ax.legend()
save_figure(fig, "05_energia")

print(f"v0 = {V0:.2f} m/s; alpha = {ALPHA_DEG:.1f} deg ({ALPHA:.6f} rad)")
print(f"Bez oporu: t_lotu = {IDEAL_FLIGHT_TIME:.6f} s")
print(f"Bez oporu: zasieg = {IDEAL_RANGE:.6f} m; h_max = {IDEAL_MAX_HEIGHT:.6f} m")
print(
    f"Z oporem: t_lotu = {drag_times[-1]:.6f} s; "
    f"zasieg = {drag_x[-1]:.6f} m; h_max = {max(drag_y):.6f} m"
)
print(f"Blad bilansu energii: {energy_balance_error:.3e} J")
print(f"Zapisano 5 wykresow w: {OUT}")
