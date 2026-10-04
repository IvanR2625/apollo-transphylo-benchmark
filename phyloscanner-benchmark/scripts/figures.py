"""
figures.py
----------
Generates four publication figures for the phyloscanner vs TransPhylo benchmark.

Fig 1 — The trade-off: population size vs who-infected-whom accuracy
Fig 2 — Within-host diversity: why minimum distance is robust to recombination
Fig 3 — Mechanism diagram: parallel failure/success pathways
Fig 4 — Pairwise distance distributions: min vs mean under recombination
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from simulate_outbreak import OutbreakSimulator, SequenceSimulator
from simulate_within_host import (
    WithinHostSimulator,
    compute_pairwise_min_distances,
    compute_pairwise_mean_distances,
)

# ── colours ───────────────────────────────────────────────────────────────────
BLUE   = '#4C72B0'
RED    = '#C44E52'
GREEN  = '#55A868'
ORANGE = '#DD8452'
GREY   = '#888888'
LGREY  = '#DDDDDD'

# ── simulation setup (shared across figures) ──────────────────────────────────
POPULATION_SIZE = 300
MAX_INFECTED    = 90
R0              = 2.5
GENOME_LENGTH   = 701
MUTATION_RATE   = 3e-3
SEED            = 42
N_SAMPLE        = 55
DAYS_IN_YEAR    = 365.25
from datetime import datetime
START_DATE = datetime(1993, 5, 11)

TP_CAT1_TRUE, TP_CAT1_INF   = 77, 329
TP_CAT2_TRUE, TP_CAT2_INF   = 79, 80
TP_CAT1_F1, TP_CAT2_F1      = 0.0, 1.0   # from paper: 0% / ~100% who-infected-whom


def _setup():
    outbreak_sim = OutbreakSimulator(population_size=POPULATION_SIZE, R0=R0, seed=SEED)
    G, inf_dates, ltfu_map = outbreak_sim.simulate(
        max_infected=MAX_INFECTED, start_date=START_DATE
    )
    rng     = np.random.default_rng(99)
    sampled = set(rng.choice(list(G.nodes()), size=min(N_SAMPLE, len(G)), replace=False))
    wh_sim  = WithinHostSimulator()
    wh_with = wh_sim.generate(G, inf_dates, sampled, with_recombination=True)
    wh_no   = wh_sim.generate(G, inf_dates, sampled, with_recombination=False)
    return G, inf_dates, sampled, wh_with, wh_no


# ── Fig 1 — The trade-off ─────────────────────────────────────────────────────

def fig1_tradeoff(ps_f1_cat1, ps_f1_cat2):
    """
    Panel A: population size — TransPhylo drastically overestimates; phyloscanner N/A.
    Panel B: who-infected-whom F1 — phyloscanner outperforms TransPhylo under recombination.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        'phyloscanner vs TransPhylo: the fundamental trade-off\n'
        'Apollo ground-truth benchmark  |  van Marle et al. (2025) parameters',
        fontsize=12, fontweight='bold',
    )

    # ── panel A: network size ─────────────────────────────────────────────────
    ax = axes[0]
    x  = np.array([0, 1])
    w  = 0.28

    cats = ['With\nRecombination', 'Without\nRecombination']
    true_vals = [TP_CAT1_TRUE,   TP_CAT2_TRUE]
    tp_vals   = [TP_CAT1_INF,    TP_CAT2_INF]

    ax.bar(x - w, true_vals, w * 2, color=BLUE,  alpha=0.88, label='Apollo ground truth')
    ax.bar(x,     tp_vals,   w * 2, color=RED,   alpha=0.88, label='TransPhylo inferred')

    # phyloscanner: hatched "not estimated" bars
    ax.bar(x + w, [0.5, 0.5], w * 2, color=LGREY, alpha=0.5,
           hatch='///', edgecolor=GREY, linewidth=1.0, label='phyloscanner (N/A)')
    ax.text(x[0] + w, 25, 'N/A\n(not\nestimated)', ha='center', va='bottom',
            fontsize=8, color=GREY, style='italic')
    ax.text(x[1] + w, 25, 'N/A\n(not\nestimated)', ha='center', va='bottom',
            fontsize=8, color=GREY, style='italic')

    for xi, v in zip(x - w, true_vals):
        ax.text(xi, v + 5, str(v), ha='center', fontsize=12,
                fontweight='bold', color=BLUE)
    for xi, v in zip(x, tp_vals):
        ax.text(xi, v + 5, str(v), ha='center', fontsize=12,
                fontweight='bold', color=RED)

    ax.annotate('4.27× overestimate', xy=(x[0], TP_CAT1_INF),
                xytext=(0.4, 260), fontsize=10, color=RED, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=RED, lw=1.5))

    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=11)
    ax.set_ylabel('Network size (infected individuals)', fontsize=11)
    ax.set_title('(A)  Network size estimation', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='upper right')
    ax.set_ylim(0, 380)
    ax.grid(True, axis='y', alpha=0.25, linestyle='--')

    # ── panel B: F1 score ─────────────────────────────────────────────────────
    ax2 = axes[1]
    tp_f1s = [TP_CAT1_F1, TP_CAT2_F1]
    ps_f1s = [ps_f1_cat1,  ps_f1_cat2]

    b1 = ax2.bar(x - w, tp_f1s, w * 2, color=RED,   alpha=0.88,
                 label='TransPhylo (from paper)')
    b2 = ax2.bar(x + w, ps_f1s, w * 2, color=GREEN, alpha=0.88,
                 label='phyloscanner (analytical model)')

    for bar, v in zip(b1, tp_f1s):
        if v > 0.02:
            ax2.text(bar.get_x() + bar.get_width()/2, v + 0.02,
                     f'{v:.2f}', ha='center', fontsize=11, fontweight='bold', color=RED)
    for bar, v in zip(b2, ps_f1s):
        ax2.text(bar.get_x() + bar.get_width()/2, v + 0.02,
                 f'{v:.2f}', ha='center', fontsize=11, fontweight='bold', color=GREEN)

    ax2.annotate('TransPhylo: 0 correct\nlinks under recombination',
                 xy=(x[0] - w, 0.0), xytext=(-0.2, 0.3), fontsize=9, color=RED,
                 arrowprops=dict(arrowstyle='->', color=RED))
    ax2.annotate('phyloscanner\nmore robust',
                 xy=(x[0] + w, ps_f1_cat1), xytext=(0.55, ps_f1_cat1 + 0.15),
                 fontsize=9, color=GREEN,
                 arrowprops=dict(arrowstyle='->', color=GREEN))

    ax2.set_xticks(x); ax2.set_xticklabels(cats, fontsize=11)
    ax2.set_ylabel('Who-infected-whom F1-score', fontsize=11)
    ax2.set_title('(B)  Transmission link accuracy', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.set_ylim(0, 1.25)
    ax2.grid(True, axis='y', alpha=0.25, linestyle='--')

    out = os.path.join('output', 'figures', 'fig1_tradeoff.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


# ── Fig 2 — Within-host diversity ─────────────────────────────────────────────

def fig2_within_host_diversity(G, inf_dates, sampled, wh_with, wh_no):
    """
    Scatter: each host = vertical cluster of N_WITHIN_HOST divergence values.
    Shows why minimum distance is robust: min stays near the true clock even when
    recombinant sequences push max/mean upward.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=False)
    fig.suptitle(
        'Within-host viral diversity: why minimum distance is robust to recombination',
        fontsize=12, fontweight='bold',
    )

    root = next(h for h in G.nodes() if G.in_degree(h) == 0)

    for ax, wh_seqs, title, recomb in [
        (axes[0], wh_no,   'Without Recombination\n(tight within-host clouds)', False),
        (axes[1], wh_with, 'With Recombination\n(outliers elevate mean and max)', True),
    ]:
        hosts_plot = [h for h in sampled if h in wh_seqs]
        days_arr, min_arr, mea_arr, all_x, all_y = [], [], [], [], []

        for h in hosts_plot:
            d  = max(0, (inf_dates[h] - inf_dates[root]).days)
            seqs = wh_seqs[h]
            days_arr.append(d)
            min_arr.append(seqs.min())
            mea_arr.append(seqs.mean())
            for s in seqs:
                all_x.append(d)
                all_y.append(s)

        ax.scatter(all_x, all_y, alpha=0.25, s=10, color=BLUE, label='Within-host sequences')
        ax.scatter(days_arr, min_arr, s=40, color=GREEN, zorder=4,
                   label='Min per host (phyloscanner)')
        ax.scatter(days_arr, mea_arr, s=40, color=ORANGE, marker='^', zorder=4,
                   label='Mean per host (consensus tools)')

        # Clock expectation line
        xf    = np.linspace(0, max(days_arr), 100) if days_arr else np.linspace(0, 100, 100)
        clock = MUTATION_RATE * xf / DAYS_IN_YEAR
        ax.plot(xf, clock, color=RED, lw=2, linestyle='--', label='Clock prediction', zorder=5)

        ax.set_xlabel('Days since index case', fontsize=11)
        ax.set_ylabel('Divergence from root (subs/site)', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=8, loc='upper left')
        ax.grid(True, alpha=0.2, linestyle='--')

    out = os.path.join('output', 'figures', 'fig2_within_host_diversity.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


# ── Fig 3 — Mechanism diagram ─────────────────────────────────────────────────

def fig3_mechanism():
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.axis('off')
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    def box(ax, x, y, text, color, w=0.16, h=0.28):
        rect = mpatches.FancyBboxPatch(
            (x - w/2, y - h/2), w, h,
            boxstyle='round,pad=0.02', fc=color, ec='white', lw=1.5,
        )
        ax.add_patch(rect)
        ax.text(x, y, text, ha='center', va='center',
                fontsize=9, fontweight='bold', color='white',
                wrap=True, multialignment='center')

    def arrow(ax, x1, x2, y):
        ax.annotate('', xy=(x2 - 0.085, y), xytext=(x1 + 0.085, y),
                    arrowprops=dict(arrowstyle='->', color='#444', lw=2.0))

    # ── TransPhylo pathway (top) ──────────────────────────────────────────────
    y_tp = 0.72
    steps_tp = [
        (0.10, 'Consensus\nsequence\n(1 per host)', '#5C6BC0'),
        (0.29, 'Clock model\n+ BEAST2\nphylogenetics', '#3949AB'),
        (0.48, 'MRCA inferred\n~1990\n(3-yr error)', '#E53935'),
        (0.67, 'Birth-death\nfills 3 extra\nyears', '#C62828'),
        (0.86, '329 inferred\nvs 77 true\n(4.27×)', '#B71C1C'),
    ]
    for i, (x, txt, col) in enumerate(steps_tp):
        box(ax, x, y_tp, txt, col)
        if i < len(steps_tp) - 1:
            arrow(ax, x, steps_tp[i+1][0], y_tp)

    ax.text(0.01, y_tp, 'TransPhylo', ha='left', va='center',
            fontsize=10, fontweight='bold', color='#5C6BC0',
            rotation=90, transform=ax.transAxes)

    # ── phyloscanner pathway (bottom) ─────────────────────────────────────────
    y_ps = 0.28
    steps_ps = [
        (0.10, 'All within-host\nsequences\n(N per host)', '#2E7D32'),
        (0.29, 'Combined\nphylogenetic\ntree', '#388E3C'),
        (0.48, 'Minimum\npairwise\ndistance', '#43A047'),
        (0.67, 'Direct links\nonly (sampled\nhosts)', '#4CAF50'),
        (0.86, 'No unsampled-\nhost inflation\n✓', '#66BB6A'),
    ]
    for i, (x, txt, col) in enumerate(steps_ps):
        box(ax, x, y_ps, txt, col)
        if i < len(steps_ps) - 1:
            arrow(ax, x, steps_ps[i+1][0], y_ps)

    # Recombination side-effect on phyloscanner
    ax.annotate('Recombination\ncreates mosaic seqs\n→ some false positives',
                xy=(0.48, y_ps - 0.14), xytext=(0.48, 0.02),
                ha='center', fontsize=8, color='#E65100',
                bbox=dict(boxstyle='round', fc='#FFF3E0', ec='#E65100'),
                arrowprops=dict(arrowstyle='->', color='#E65100'))

    ax.text(0.01, y_ps, 'phyloscanner', ha='left', va='center',
            fontsize=10, fontweight='bold', color='#2E7D32',
            rotation=90, transform=ax.transAxes)

    # Recombination source (top-left, causes both)
    ax.text(0.5, 0.97,
            'Recombination present (14 hotspots, Category 1)',
            ha='center', va='top', fontsize=11, color='#4A148C', fontweight='bold',
            bbox=dict(boxstyle='round', fc='#EDE7F6', ec='#4A148C'))

    ax.set_title('Parallel causal pathways: why TransPhylo fails and phyloscanner is more robust',
                 fontsize=12, fontweight='bold', pad=4)

    out = os.path.join('output', 'figures', 'fig3_mechanism.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


# ── Fig 4 — Pairwise distance distributions ───────────────────────────────────

def fig4_distance_distributions(G, inf_dates, sampled, wh_with, wh_no):
    """
    Min vs mean pairwise distances for true transmission pairs vs non-pairs.
    Recombination makes mean distributions overlap; min distributions stay separated.
    """
    import networkx as nx
    true_edges = {(u, v) for u, v in G.edges() if u in sampled and v in sampled}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        'Pairwise sequence distances: minimum vs mean across recombination conditions',
        fontsize=12, fontweight='bold',
    )

    for ax, wh_seqs, title in [
        (axes[0], wh_no,   'Without Recombination'),
        (axes[1], wh_with, 'With Recombination'),
    ]:
        hosts_list = sorted(h for h in sampled if h in wh_seqs)
        min_true, min_non, mea_true, mea_non = [], [], [], []

        for i, ha in enumerate(hosts_list):
            for j, hb in enumerate(hosts_list):
                if i >= j:
                    continue
                sa, sb = wh_seqs[ha], wh_seqs[hb]
                diff   = np.abs(sa[:, None] - sb[None, :])
                mn     = diff.min()
                me     = abs(sa.mean() - sb.mean())
                if (ha, hb) in true_edges or (hb, ha) in true_edges:
                    min_true.append(mn); mea_true.append(me)
                else:
                    min_non.append(mn);  mea_non.append(me)

        bins = np.linspace(0, max(
            (max(min_non) if min_non else 0.05),
            (max(mea_non) if mea_non else 0.05)
        ), 30)

        ax.hist(min_non,  bins=bins, color=BLUE,   alpha=0.45, density=True,
                label='Min dist — non-pairs')
        ax.hist(min_true, bins=bins, color=GREEN,  alpha=0.70, density=True,
                label='Min dist — true pairs')
        ax.hist(mea_non,  bins=bins, color=RED,    alpha=0.35, density=True,
                linestyle='--', histtype='step', linewidth=1.8,
                label='Mean dist — non-pairs')
        ax.hist(mea_true, bins=bins, color=ORANGE, alpha=0.70, density=True,
                linestyle='--', histtype='step', linewidth=1.8,
                label='Mean dist — true pairs')

        ax.set_xlabel('Pairwise sequence distance (subs/site)', fontsize=11)
        ax.set_ylabel('Density', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2, linestyle='--')

        if len(min_true) > 1 and len(min_non) > 1:
            ax.text(0.97, 0.97,
                    f'Min dist separation:\n'
                    f'  true pairs:  {np.mean(min_true):.4f}\n'
                    f'  non-pairs:   {np.mean(min_non):.4f}\n'
                    f'Mean dist separation:\n'
                    f'  true pairs:  {np.mean(mea_true):.4f}\n'
                    f'  non-pairs:   {np.mean(mea_non):.4f}',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=8, bbox=dict(boxstyle='round', fc='white', alpha=0.85))

    out = os.path.join('output', 'figures', 'fig4_distance_distributions.png')
    plt.tight_layout()
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {out}')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    os.makedirs(os.path.join('output', 'figures'), exist_ok=True)

    print('Setting up simulation …')
    G, inf_dates, sampled, wh_with, wh_no = _setup()

    # Run phyloscanner analytical model to get F1 scores
    sys.path.insert(0, os.path.dirname(__file__))
    from run_phyloscanner_analytical import infer_phyloscanner
    _, _, ps_prec_c1, ps_rec_c1, ps_f1_c1 = infer_phyloscanner(
        G, inf_dates, wh_with, sampled, with_recombination=True)
    _, _, ps_prec_c2, ps_rec_c2, ps_f1_c2 = infer_phyloscanner(
        G, inf_dates, wh_no, sampled, with_recombination=False)

    print(f'phyloscanner F1: Cat1={ps_f1_c1:.2f}  Cat2={ps_f1_c2:.2f}')

    print('Generating Fig 1 …'); fig1_tradeoff(ps_f1_c1, ps_f1_c2)
    print('Generating Fig 2 …'); fig2_within_host_diversity(G, inf_dates, sampled, wh_with, wh_no)
    print('Generating Fig 3 …'); fig3_mechanism()
    print('Generating Fig 4 …'); fig4_distance_distributions(G, inf_dates, sampled, wh_with, wh_no)

    print('\nAll figures saved to output/figures/')


if __name__ == '__main__':
    main()
