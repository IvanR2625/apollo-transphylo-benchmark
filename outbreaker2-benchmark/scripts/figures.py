"""
figures.py
──────────
Generate the three publication figures for the outbreaker2 benchmark.

Usage
-----
    python scripts/figures.py

Imports the analytical benchmark functions directly; does not re-run the
full simulation (uses the same seed, so results are deterministic).

Outputs
-------
  output/figures/fig1_network_size.png
  output/figures/fig2_mechanism.png
  output/figures/fig3_snp_scatter.png
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, os.path.dirname(__file__))

from simulate_outbreak import OutbreakSimulator, SequenceSimulator  # noqa
from run_outbreaker2_analytical import (                              # noqa
    POPULATION_SIZE, MAX_INFECTED, R0, GENERATION_TIME, START_DATE,
    N_SAMPLE, RECOMB_HOTSPOTS, GENOME_LENGTH, MU, SEED, SAMPLING_SEED,
    pairwise_snp_distances, infer_outbreaker2, infer_transphylo,
    BLUE, RED, ORANGE, GREEN, GREY, PURPLE,
)

OUTDIR = os.path.join("output", "figures")

# Published TransPhylo values from van Marle et al. (2025)
TP_C1_TRUE, TP_C1_PAPER = 77, 329
TP_C2_TRUE, TP_C2_PAPER = 79, 80


def _simulate(with_recombination):
    sim = OutbreakSimulator(population_size=POPULATION_SIZE, R0=R0,
                            gen_time=GENERATION_TIME, seed=SEED)
    G, inf_dates, ltfu_map = sim.simulate(max_infected=MAX_INFECTED,
                                          start_date=START_DATE)
    rng = np.random.default_rng(SAMPLING_SEED)
    sampled = set(rng.choice(list(G.nodes()),
                              size=min(N_SAMPLE, len(G.nodes())),
                              replace=False).tolist())
    seq = SequenceSimulator(mu=MU, L=GENOME_LENGTH,
                            n_recomb=RECOMB_HOTSPOTS, seed=SEED + 1)
    divs = seq.generate(G, inf_dates, with_recombination=with_recombination)
    return G, inf_dates, sampled, divs


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Network size: truth vs TransPhylo vs outbreaker2
# ─────────────────────────────────────────────────────────────────────────────

def fig1_network_size():
    results = {}
    for recomb in [True, False]:
        G, inf_dates, sampled, divs = _simulate(recomb)
        dists = pairwise_snp_distances(sampled, divs, G, inf_dates)
        ob2, _, _ = infer_outbreaker2(dists, len(sampled))
        tp, _ = infer_transphylo(divs, inf_dates, len(sampled))
        results[recomb] = dict(true=len(G), ob2=ob2, tp=tp,
                               true_paper=TP_C1_TRUE if recomb else TP_C2_TRUE,
                               tp_paper=TP_C1_PAPER if recomb else TP_C2_PAPER)

    fig, ax = plt.subplots(figsize=(11, 6))
    fig.suptitle(
        "Transmission network inflation under recombination\n"
        "TransPhylo vs outbreaker2 — Apollo ground truth benchmark",
        fontsize=13, fontweight="bold",
    )

    x = np.array([0, 1])
    w = 0.22
    labels = ["With Recombination\n(Cat. 1)", "Without Recombination\n(Cat. 2)"]

    truth_v  = [results[True]["true"],  results[False]["true"]]
    tp_sim_v = [results[True]["tp"],    results[False]["tp"]]
    tp_pap_v = [TP_C1_PAPER,            TP_C2_PAPER]
    ob2_v    = [results[True]["ob2"],   results[False]["ob2"]]

    b1 = ax.bar(x - 1.5*w, truth_v,  w, color=BLUE,   alpha=0.90, label="Apollo ground truth")
    b2 = ax.bar(x - 0.5*w, tp_pap_v, w, color=RED,    alpha=0.90, label="TransPhylo (van Marle 2025)")
    b3 = ax.bar(x + 0.5*w, tp_sim_v, w, color=ORANGE, alpha=0.90, label="TransPhylo (this simulation)")
    b4 = ax.bar(x + 1.5*w, ob2_v,    w, color=PURPLE, alpha=0.90, label="outbreaker2 (this study)")

    for bars, vals, color in [(b1, truth_v, BLUE), (b2, tp_pap_v, RED),
                               (b3, tp_sim_v, ORANGE), (b4, ob2_v, PURPLE)]:
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, v + 4, str(v),
                    ha="center", fontsize=10, fontweight="bold", color=color)

    # Ratio annotations
    for xi, r_tp, r_ob2 in zip(x, [TP_C1_PAPER/TP_C1_TRUE, TP_C2_PAPER/TP_C2_TRUE],
                                   [ob2_v[0]/truth_v[0], ob2_v[1]/truth_v[1]]):
        ax.text(xi - 0.5*w, tp_pap_v[0 if xi == 0 else 1] + 18,
                f"{r_tp:.2f}×", ha="center", fontsize=10, color=RED, fontweight="bold")
        ax.text(xi + 1.5*w, ob2_v[0 if xi == 0 else 1] + 18,
                f"{r_ob2:.2f}×", ha="center", fontsize=10, color=PURPLE, fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("Inferred infected individuals", fontsize=12)
    ax.set_ylim(0, 400)
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, axis="y", alpha=0.25, linestyle="--")

    ax.text(0.01, 0.98,
            "Both tools overestimate under recombination\n"
            "TransPhylo (global MRCA) inflates more severely\n"
            "than outbreaker2 (local SNP distances)",
            transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="lightyellow", ec="goldenrod", alpha=0.9))

    path = os.path.join(OUTDIR, "fig1_network_size.png")
    plt.tight_layout()
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")
    return results


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Mechanism comparison
# ─────────────────────────────────────────────────────────────────────────────

def fig2_mechanism():
    fig, axes = plt.subplots(2, 1, figsize=(14, 7))
    fig.suptitle(
        "Two distinct pathways by which recombination inflates transmission network estimates",
        fontsize=13, fontweight="bold",
    )

    tp_steps = [
        ("Recombination\nevents", "#FF7043"),
        ("Sequences more\ndiverged than\nclock predicts", "#F4511E"),
        ("BEAST2 infers\nolder MRCA\n(1990 vs 1993)", "#E53935"),
        ("TransPhylo fills\nextra years with\nphantom hosts", "#C62828"),
        ("329 inferred\nvs 77 true\n(4.27×)", "#B71C1C"),
    ]
    ob2_steps = [
        ("Recombination\nevents", "#5E35B1"),
        ("Pairwise SNP\ndistances\ninflated", "#512DA8"),
        ("outbreaker2 infers\nextra transmission\nsteps: t̂ = d/(μL)", "#4527A0"),
        ("Extra unsampled\nintermediates\ninferred per link", "#311B92"),
        ("~1.8× inferred\nvs 77 true\n(less severe)", "#1A237E"),
    ]

    for ax, steps, tool_label, bg in [
        (axes[0], tp_steps,  "TransPhylo pathway (MRCA-based inference)",  "#FFF3E0"),
        (axes[1], ob2_steps, "outbreaker2 pathway (pairwise SNP distances)", "#EDE7F6"),
    ]:
        ax.set_facecolor(bg)
        ax.axis("off")
        ax.set_title(tool_label, fontsize=11, fontweight="bold", loc="left", pad=6)

        for i, (label, color) in enumerate(steps):
            xc = 0.08 + i * 0.21
            ax.text(xc, 0.52, label, transform=ax.transAxes,
                    ha="center", va="center", fontsize=10, fontweight="bold",
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.45", fc=color, ec="white", lw=1.2))
            if i < len(steps) - 1:
                ax.annotate("",
                            xy=(xc + 0.13, 0.52), xytext=(xc + 0.07, 0.52),
                            xycoords="axes fraction", textcoords="axes fraction",
                            arrowprops=dict(arrowstyle="->", color="#444", lw=2.0))

        ax.text(0.5, 0.10,
                "Without recombination: sequences match clock → correct inference",
                transform=ax.transAxes, ha="center", va="center", fontsize=10,
                color="#1565C0",
                bbox=dict(boxstyle="round,pad=0.35", fc="#E3F2FD",
                          ec="#1565C0", lw=1.2))

    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig2_mechanism.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — SNP distance vs transmission distance
# ─────────────────────────────────────────────────────────────────────────────

def fig3_snp_scatter():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Pairwise SNP distance vs. transmission steps: clock-like vs. recombination",
        fontsize=13, fontweight="bold",
    )

    for ax, recomb, title, color in [
        (axes[0], False, "Without Recombination\n(clock-like)", BLUE),
        (axes[1], True,  "With Recombination\n(SNP distances inflated)", RED),
    ]:
        G, inf_dates, sampled, divs = _simulate(recomb)
        dists = pairwise_snp_distances(sampled, divs, G, inf_dates)

        if not dists:
            ax.set_title(title)
            continue

        true_steps = np.array([d["true_steps"] for d in dists.values()])
        snp_dists  = np.array([d["snp_dist"]   for d in dists.values()])
        expected   = true_steps * MU * GENOME_LENGTH

        # Outliers: points >2 SD above the expected clock line
        residuals = snp_dists - expected
        outlier_mask = residuals > 2 * np.std(residuals)

        ax.scatter(true_steps[~outlier_mask], snp_dists[~outlier_mask],
                   c=color, s=55, alpha=0.7, zorder=3,
                   label=f"In-clock ({(~outlier_mask).sum()})")
        if outlier_mask.any():
            ax.scatter(true_steps[outlier_mask], snp_dists[outlier_mask],
                       c=ORANGE, s=80, marker="D", alpha=0.9, zorder=4,
                       label=f"Recomb outliers ({outlier_mask.sum()})")

        # Expected clock line
        xmax = true_steps.max() + 0.5
        xs = np.linspace(0, xmax, 100)
        ax.plot(xs, xs * MU * GENOME_LENGTH, "--", color="black", lw=1.5,
                label=f"Clock: E[SNPs] = {MU}×{GENOME_LENGTH}×steps")

        # Linear fit to actual data
        if len(true_steps) > 2:
            sl, ic, rv, _, _ = stats.linregress(true_steps, snp_dists)
            ax.plot(xs, sl * xs + ic, "-", color=color, lw=2, alpha=0.6,
                    label=f"Observed fit  R²={rv**2:.2f}")

        ax.set_xlabel("True transmission steps", fontsize=11)
        ax.set_ylabel("Pairwise SNP distance (sites)", fontsize=11)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.25, linestyle="--")

        if recomb:
            n_out = outlier_mask.sum()
            pct = 100 * n_out / len(snp_dists)
            ax.text(0.97, 0.05,
                    f"{n_out} outliers ({pct:.0f}% of links)\n"
                    "drive outbreaker2 inflation",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=9, color=ORANGE,
                    bbox=dict(boxstyle="round", fc="white", ec=ORANGE, alpha=0.9))

    plt.tight_layout()
    path = os.path.join(OUTDIR, "fig3_snp_scatter.png")
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {path}")


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("Generating figures …\n")
    results = fig1_network_size()
    fig2_mechanism()
    fig3_snp_scatter()
    print("\nAll figures written to", OUTDIR)
    return results


if __name__ == "__main__":
    main()
