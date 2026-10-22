# %% [markdown]
# # When does checking pay before an irreversible step?
#
# Code for the section "When gathering information uses up the chance to correct". It computes the value of checking
# in the paper's deadline model, draws Figure 1, and checks numerically the closed form for that value, each part of
# the proposition, and the statements about postponing the step.
#
# Run it as a script (`python recoverability_deadline_model.py [--out DIR]`) or cell by cell in an editor that reads
# `# %%` cells. It needs NumPy and Matplotlib.

# %% [markdown]
# ## Model
#
# - A yes-or-no fact $Y$ decides whether committing is right, with $P(Y=1)=\pi$.
# - The decision maker's own signal $X_0$ equals $Y$ with probability $q_0$, so on its own it errs with probability $e_0$.
# - A check returns a second signal. With probability $\iota$ (its independence) it is fresh and equals $Y$ with
#   probability $q_c$. Otherwise it repeats $X_0$, error included. The decision maker does not know which happened.
# - The check arrives before the deadline with probability $p$. If it is late, the decision is made on $X_0$ alone.
# - A wrong decision costs 1 if it cannot be corrected and $\kappa<1$ if it can. A correction needs one unit of material
#   and the check uses a share $c$ of it, so after checking a wrong decision costs
#   $L(c) = 1-(1-\kappa)\min(1,\max(0,1-c))$. The check itself costs $k$.
# - A fresh check that arrives in time moves the deadline with probability $\mu_y$ when $Y=y$, and the move is seen.
#   A move returns a share $g$ of the material, so a wrong decision then costs $L(c-g)$. A check that repeats $X_0$
#   never moves the deadline. With $\mu_0=\mu_1=0$ the deadline is fixed.
#
# With $e^*$ the smallest error after a timely check (using $X_0$, the check and whether the deadline moved) and
# $e^*_D$ the probability of a wrong decision in the cases where it moved, checking lowers the expected loss by
# $$V = e_0 L(0) - \bigl[p\,\bigl(e^* L(c) - e^*_D (L(c)-L(c-g))\bigr) + (1-p)\,e_0 L(c)\bigr] - k .$$
#
# All parameter values below are assumed and illustrative. The proposition does not depend on them.

# %%
import argparse
import functools
import itertools
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

PARAMS = dict(pi=0.5, q0=0.65, qc=0.90, kappa=0.30, k=0.02)
NO_MOVE = (0.0, 0.0)
TOL = 1e-12


def default_out():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


# %% [markdown]
# ## Exact errors and the value of checking

# %%
OBS = [(x, ch, d) for x in (0, 1) for ch in (0, 1) for d in (0, 1)]


def L(c, kappa):
    return 1 - (1 - kappa) * min(1.0, max(0.0, 1 - c))


def joint_timely(iota, mu, P):
    """Joint probability of (Y, X0, check, move) after a timely check."""
    J = {}
    for y in (0, 1):
        py = P["pi"] if y else 1 - P["pi"]
        for x in (0, 1):
            px = P["q0"] if x == y else 1 - P["q0"]
            for ch in (0, 1):
                for d in (0, 1):
                    fresh = iota * (P["qc"] if ch == y else 1 - P["qc"]) * (mu[y] if d else 1 - mu[y])
                    copy = (1 - iota) * (ch == x) * (d == 0)
                    J[(y, x, ch, d)] = py * px * (fresh + copy)
    return J


def error_without_check(P):
    return sum(min((P["pi"] if y else 1 - P["pi"]) * (P["q0"] if x == y else 1 - P["q0"]) for y in (0, 1))
               for x in (0, 1))


def errors(iota, mu, P):
    """Return e0, e* and e*_D."""
    return _errors(float(iota), float(mu[0]), float(mu[1]), P["pi"], P["q0"], P["qc"])


@functools.lru_cache(maxsize=None)
def _errors(iota, mu0, mu1, pi, q0, qc):
    P = dict(pi=pi, q0=q0, qc=qc)
    J = joint_timely(iota, (mu0, mu1), P)
    e_star = sum(min(J[(0,) + o], J[(1,) + o]) for o in OBS)
    e_moved = sum(min(J[(0,) + o], J[(1,) + o]) for o in OBS if o[2] == 1)
    return error_without_check(P), e_star, e_moved


def value(iota, p, c, P=PARAMS, g=0.0, mu=NO_MOVE):
    """Closed-form reduction in expected loss from checking (positive: checking pays)."""
    e0, e_star, e_moved = errors(iota, mu, P)
    Lc, L0, Lcg = L(c, P["kappa"]), L(0.0, P["kappa"]), L(c - g, P["kappa"])
    return e0 * L0 - (p * (e_star * Lc - e_moved * (Lc - Lcg)) + (1 - p) * e0 * Lc) - P["k"]


def condition_margin(iota, p, c, P=PARAMS, g=0.0, mu=NO_MOVE):
    """Left side minus right side of the condition for checking to pay."""
    e0, e_star, e_moved = errors(iota, mu, P)
    Lc, L0, Lcg = L(c, P["kappa"]), L(0.0, P["kappa"]), L(c - g, P["kappa"])
    return p * (e0 - e_star) * Lc + p * e_moved * (Lc - Lcg) - (P["k"] + e0 * (Lc - L0))


def value_by_search(iota, p, c, P=PARAMS, g=0.0, mu=NO_MOVE):
    """The same value, found by searching every deterministic decision rule."""
    J = joint_timely(iota, mu, P)
    kappa = P["kappa"]
    # loss of each rule = sum over observations of the probability of the wrong Y times the cost of an error there
    cost = np.array([L(c - g * o[2], kappa) for o in OBS])
    wrong_if = {a: np.array([J[(1 - a,) + o] for o in OBS]) * cost for a in (0, 1)}
    timely = min_over_rules(wrong_if, len(OBS))
    Jx = {a: np.array([(P["pi"] if 1 - a else 1 - P["pi"]) * (P["q0"] if x == 1 - a else 1 - P["q0"])
                       for x in (0, 1)]) for a in (0, 1)}
    late = min_over_rules({a: Jx[a] * L(c, kappa) for a in (0, 1)}, 2)
    none = min_over_rules({a: Jx[a] * L(0.0, kappa) for a in (0, 1)}, 2)
    return none - (p * timely + (1 - p) * late) - P["k"]


RULES = {n: np.array(list(itertools.product((0, 1), repeat=n))) for n in (2, len(OBS))}


def min_over_rules(wrong_if, n):
    """Smallest expected loss over all 2**n rules; wrong_if[a][j] is the loss of choosing a at observation j."""
    R = RULES[n]
    return float((R @ wrong_if[1] + (1 - R) @ wrong_if[0]).min())


# %% [markdown]
# ## Simulation against the exact values

# %%
def simulate(iota, p, c, P=PARAMS, g=0.0, mu=NO_MOVE, n=200_000, seed=0):
    rng = np.random.default_rng(seed)
    y = (rng.random(n) < P["pi"]).astype(int)
    x0 = np.where(rng.random(n) < P["q0"], y, 1 - y)
    fresh = rng.random(n) < iota
    check = np.where(fresh, np.where(rng.random(n) < P["qc"], y, 1 - y), x0)
    on_time = rng.random(n) < p
    moved = (fresh & on_time & (rng.random(n) < np.asarray(mu)[y])).astype(int)
    J = joint_timely(iota, mu, P)
    rule_timely = {o: int(J[(1,) + o] > J[(0,) + o]) for o in OBS}
    rule_late = {x: int(P["pi"] * (P["q0"] if x else 1 - P["q0"]) > (1 - P["pi"]) * (1 - P["q0"] if x else P["q0"]))
                 for x in (0, 1)}
    a_timely = np.array([rule_timely[o] for o in zip(x0, check, moved)])
    a_late = np.array([rule_late[x] for x in x0])
    a = np.where(on_time, a_timely, a_late)
    cost = np.where(moved == 1, L(c - g, P["kappa"]), L(c, P["kappa"]))
    loss_check = ((a != y) * cost).mean() + P["k"]
    a_none = np.array([rule_late[x] for x in x0])
    loss_none = (a_none != y).mean() * L(0.0, P["kappa"])
    return loss_none - loss_check


# %% [markdown]
# ## Postponing the step
#
# The check needs a minimum time `MIN_DELAY` and then a further wait, exponential with mean `mean_delay`, so it arrives
# by time $T$ with probability $p(T) = 1-e^{-(T-d)/m}$ for $T>d$ and 0 otherwise. Postponing the step by $\Delta$ costs
# $h$ per unit of time and leaves a share $e^{-\delta\Delta}$ of the material. Times are in arbitrary units.

# %%
MIN_DELAY = 20.0
MEAN_DELAY = 6.0
POSTPONE = dict(h=0.0005, delta=0.01)
WAITS = np.linspace(0, 72, 289)


def p_on_time(T, mean_delay=MEAN_DELAY):
    return np.where(T > MIN_DELAY, 1 - np.exp(-(np.asarray(T, dtype=float) - MIN_DELAY) / mean_delay), 0.0)


def loss_after_wait(iota, T, c, wait, P, Q, mean_delay=MEAN_DELAY):
    e0, e, _ = errors(iota, NO_MOVE, P)
    pw = p_on_time(T + wait, mean_delay)
    f = np.clip(np.exp(-Q["delta"] * wait) - c, 0.0, 1.0)
    return (pw * e + (1 - pw) * e0) * (1 - (1 - P["kappa"]) * f) + P["k"] + Q["h"] * wait


def best_wait(iota, T, c, P, Q, mean_delay=MEAN_DELAY):
    losses = loss_after_wait(iota, T, c, WAITS, P, Q, mean_delay)
    i = int(np.argmin(losses))
    return float(WAITS[i]), float(losses[i])


# %% [markdown]
# ## Figure 1. Where an independent check pays, as it uses more of the correction material (fixed deadline)

# %%
def figure_1(out):
    iotas = np.linspace(0, 1, 101)
    ps = np.linspace(0, 1, 101)
    cs = [0.0, 0.5, 1.0]
    Vs = {c: np.array([[value(i, p, c) for i in iotas] for p in ps]) for c in cs}
    norm = TwoSlopeNorm(vcenter=0.0, vmin=min(V.min() for V in Vs.values()), vmax=max(V.max() for V in Vs.values()))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    panels = {}
    for ax, c in zip(axes, cs):
        V = Vs[c]
        im = ax.imshow(V, origin="lower", extent=[0, 1, 0, 1], aspect="auto", cmap="RdBu", norm=norm)
        if V.max() > 0:
            ax.contour(iotas, ps, V, levels=[0], colors="k", linewidths=1.2)
        ax.set_title(f"c = {c}")
        ax.set_xlabel("independence of the check  ι")
        panels[str(c)] = {"share_of_grid_where_check_pays": float((V > 0).mean()), "max_value": float(V.max())}
    axes[0].set_ylabel("chance it arrives before the deadline  p")
    fig.colorbar(im, ax=axes, label="value of checking $V$ (blue: pays)")
    fig.savefig(os.path.join(out, "fig1_where_checking_pays.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    return panels


# %% [markdown]
# ## Check 1. The proposition and the postponement statements over a grid of assumed values
#
# Every assumed value is given three levels and all 2,187 combinations are checked. In each combination the
# proposition's parts (i) to (iv) are checked for a fixed deadline and for three moving-deadline cases, and the closed
# form is compared with the search over every decision rule and with the condition. The postponement statements are checked for a fixed
# deadline: the option to postpone never makes the best loss worse, and the best delay is at most $e_0L(0)/h$.

# %%
SWEEP = dict(q0=[0.55, 0.65, 0.75], qc=[0.80, 0.90, 0.97], kappa=[0.1, 0.3, 0.5], k=[0.01, 0.02, 0.05],
             mean_delay=[3.0, 6.0, 12.0], h=[0.0001, 0.0005, 0.002], delta=[0.0, 0.01, 0.05])
MOVES = [(NO_MOVE, 0.0), ((0.0, 1.0), 0.5), ((0.3, 0.6), 0.25), ((0.5, 0.5), 0.5)]


def proposition_holds(P, g, mu):
    grid = np.linspace(0, 1, 11)
    ok = {}
    ok["(i)"] = all(value(0.0, p, c, P, g, mu) <= -P["k"] + TOL for p in grid for c in (0.0, 0.5, 1.0))
    V = np.array([[value(i, p, 0.5, P, g, mu) for i in grid] for p in grid])
    ok["(ii) independence"] = bool((np.diff(V, axis=1) >= -TOL).all())
    ok["(ii) arrival"] = bool((np.diff(V, axis=0) >= -TOL).all())
    ok["(ii) deadline move"] = all(value(1.0, 1.0, 0.5, P, g1, mu) <= value(1.0, 1.0, 0.5, P, g2, mu) + TOL
                                   for g1, g2 in [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0)])
    ok["(iii)"] = all(value(1.0, 1.0, c1, P, g, mu) >= value(1.0, 1.0, c2, P, g, mu) - TOL
                      for c1, c2 in [(0.0, 0.25), (0.25, 0.5), (0.5, 1.0)])
    ok["(iv)"] = all(value(i, 0.0, c, P, g, mu) < 0 for i in grid for c in (0.0, 0.5, 1.0))
    ok["closed form = search"] = all(abs(value(i, p, c, P, g, mu) - value_by_search(i, p, c, P, g, mu)) <= TOL
                                     for i, p, c in [(0.3, 0.7, 0.2), (1.0, 1.0, 0.5), (0.8, 0.4, 0.9)])
    ok["condition = value"] = all(abs(value(i, p, c, P, g, mu) - condition_margin(i, p, c, P, g, mu)) <= TOL
                                  for i in grid for p in grid for c in (0.0, 0.5))
    return ok


def sweep():
    passed, n = {}, 0
    for vals in itertools.product(*SWEEP.values()):
        s = dict(zip(SWEEP.keys(), vals))
        P = dict(PARAMS, q0=s["q0"], qc=s["qc"], kappa=s["kappa"], k=s["k"])
        Q = dict(h=s["h"], delta=s["delta"])
        n += 1
        for mu, g in MOVES:
            label = "fixed deadline" if mu == NO_MOVE else f"moving deadline, mu={mu}, g={g}"
            for name, ok in proposition_holds(P, g, mu).items():
                passed[f"{label}: {name}"] = passed.get(f"{label}: {name}", 0) + int(ok)
        never_worse, bounded = True, True
        bound = error_without_check(P) * L(0.0, P["kappa"]) / Q["h"]
        for T in (5.0, 15.0, 25.0, 40.0):
            fixed = float(loss_after_wait(1.0, T, 0.0, 0.0, P, Q, s["mean_delay"]))
            wait, best = best_wait(1.0, T, 0.0, P, Q, s["mean_delay"])
            never_worse &= best <= fixed + TOL
            bounded &= wait <= bound
        passed["postponement: never worse"] = passed.get("postponement: never worse", 0) + int(never_worse)
        passed["postponement: best delay at most e0 L(0)/h"] = (
            passed.get("postponement: best delay at most e0 L(0)/h", 0) + int(bounded))
    return n, passed


# %% [markdown]
# ## Check 2. Random settings against the search over every decision rule
#
# 3,000 settings are drawn at random, with the deadline moving exactly when $Y=1$ in about 30% of them. In each, the
# closed form is compared with the search over every decision rule, and each part of the proposition is checked by a
# small step in the relevant quantity. A second set of 3,000 settings counts how often making a move more likely
# lowers the value, which is why the proposition makes no claim about $\mu$. When the move depends on $Y$ it is also
# a signal about $Y$, and raising $\mu_0$ and $\mu_1$ unevenly can make that signal noisier. A third set checks that
# scaling both chances up in the same proportion never lowers the value.

# %%
N_RANDOM = 3000


def random_setting(rng):
    q0 = rng.uniform(0.5, 0.95)
    P = dict(pi=rng.uniform(0.2, 0.8), q0=q0, qc=rng.uniform(q0, 1.0), kappa=rng.uniform(0, 1),
             k=rng.uniform(1e-4, 0.1))
    iota, p = rng.uniform(0, 1, 2)
    c, g = rng.uniform(0, 1.2, 2)
    mu = np.array([0.0, 1.0]) if rng.random() < 0.3 else rng.uniform(0, 1, 2)
    return P, iota, p, c, g, mu


def random_check():
    rng = np.random.default_rng(2)
    fails = dict.fromkeys(["closed form = search", "condition", "e* <= e0", "(i)", "(ii) arrival",
                           "(ii) independence", "(ii) deadline move", "(iii)", "(iv)"], 0)
    step = 1e-3
    for _ in range(N_RANDOM):
        P, iota, p, c, g, mu = random_setting(rng)
        V = value(iota, p, c, P, g, mu)
        e0, e_star, _ = errors(iota, mu, P)
        fails["closed form = search"] += abs(V - value_by_search(iota, p, c, P, g, mu)) > TOL
        fails["condition"] += abs(condition_margin(iota, p, c, P, g, mu) - V) > TOL
        fails["e* <= e0"] += e_star > e0 + TOL
        fails["(i)"] += value(0.0, p, c, P, g, mu) > -P["k"] + TOL
        fails["(ii) arrival"] += value(iota, min(p + step, 1), c, P, g, mu) < V - TOL
        fails["(ii) independence"] += value(min(iota + step, 1), p, c, P, g, mu) < V - TOL
        fails["(ii) deadline move"] += value(iota, p, c, P, g + step, mu) < V - TOL
        fails["(iii)"] += value(iota, p, c + step, P, g, mu) > V + TOL
        fails["(iv)"] += value(iota, 0.0, c, P, g, mu) >= 0
    lower = 0
    for _ in range(N_RANDOM):
        P, iota, p, c, g, mu = random_setting(rng)
        mu = rng.uniform(0, 0.9, 2)
        lower += value(iota, p, c, P, g, mu + 0.1 * rng.uniform(0, 1, 2)) < value(iota, p, c, P, g, mu) - TOL
    proportional_lower = 0
    for _ in range(N_RANDOM):
        P, iota, p, c, g, mu = random_setting(rng)
        mu = rng.uniform(0, 1, 2)
        scale = rng.uniform(0.5, 1.0)
        proportional_lower += value(iota, p, c, P, g, mu) < value(iota, p, c, P, g, scale * mu) - TOL
    return fails, lower, proportional_lower


# %%
def main(out):
    os.makedirs(out, exist_ok=True)
    results = {"params": PARAMS, "postponement": dict(POSTPONE, min_delay=MIN_DELAY, mean_delay=MEAN_DELAY)}

    e0, e1, _ = errors(1.0, NO_MOVE, PARAMS)
    results["e0"], results["e_independent_check"] = e0, e1
    print(f"e0 = {e0:.4f}, error with an independent check = {e1:.4f}")

    sims = []
    for iota, p, c, g, mu in [(1.0, 0.8, 0.0, 0.0, NO_MOVE), (0.5, 0.5, 0.5, 0.0, NO_MOVE),
                              (0.8, 0.9, 0.5, 0.5, (0.3, 0.6)), (1.0, 1.0, 0.5, 0.5, (0.0, 1.0))]:
        exact, sim = value(iota, p, c, PARAMS, g, mu), simulate(iota, p, c, PARAMS, g, mu)
        sims.append(dict(iota=iota, p=p, c=c, g=g, mu=list(mu), exact=exact, simulated=float(sim)))
        print(f"iota={iota}, p={p}, c={c}, g={g}, mu={mu}: exact {exact:.4f}, simulated {sim:.4f}")
    results["simulation"] = sims

    results["figure_1"] = figure_1(out)
    print("Figure 1:", json.dumps(results["figure_1"]))

    n, passed = sweep()
    results["sweep"] = {"combinations": n, "holds_in": {k: f"{v}/{n}" for k, v in passed.items()}}
    print(f"\nCheck 1: {n} combinations")
    for name, v in passed.items():
        print(f"  {v}/{n}  {name}")

    fails, lower, proportional_lower = random_check()
    results["random_settings"] = {"settings": N_RANDOM, "failures": {k: int(v) for k, v in fails.items()},
                                  "settings_where_a_likelier_move_lowers_the_value": int(lower),
                                  "settings_where_scaling_both_chances_up_lowers_the_value": int(proportional_lower)}
    print(f"\nCheck 2: {N_RANDOM} random settings, failures per check")
    for name, v in fails.items():
        print(f"  {v}  {name}")
    print(f"Settings where a likelier move lowers the value: {lower} of {N_RANDOM}")
    print(f"Settings where scaling both chances up lowers the value: {proportional_lower} of {N_RANDOM}")

    with open(os.path.join(out, "deadline_model_results.json"), "w") as f:
        json.dump(results, f, indent=1)
    print("\nSaved to", out)

    # A likelier move may lower the value (the proposition makes no claim about it); every other check must pass.
    failed = [name for name, v in passed.items() if v != n] + [name for name, v in fails.items() if v]
    if proportional_lower:
        failed.append("scaling both chances up lowers the value")
    if failed:
        sys.exit("FAILED checks: " + "; ".join(failed))
    print("All checks passed.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Deadline model: value of checking, Figure 1 and numerical checks.")
    ap.add_argument("--out", default=default_out(),
                    help="folder for the figure and the results file (default: output/ beside this script); "
                         "existing files there are overwritten")
    main(ap.parse_args().out)
