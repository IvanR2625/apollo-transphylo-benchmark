"""
run_phyloscanner_analytical.py
------------------------------
Analytical model of phyloscanner's transmission-link inference, benchmarked
against TransPhylo on Apollo-simulated outbreaks with known ground truth.

This script does NOT run the real phyloscanner software.  It implements the
statistical logic phyloscanner applies — minimum-distance thresholding with
a temporal constraint — which allows principled comparison with TransPhylo's
clock/MRCA approach without requiring a full phylogenetic analysis.

Outputs
-------
Prints a comparison table (see TEMPLATE below).
Saves output/comparison_table.csv.
"""

import os
import sys
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from simulate_outbreak import OutbreakSimulator, SequenceSimulator
from simulate_within_host import (
    WithinHostSimulator,
    compute_pairwise_min_distances,
    compute_pairwise_mean_distances,
)

# ── simulation parameters (Table 1, van Marle et al. 2025) ───────────────────
POPULATION_SIZE  = 300
MAX_INFECTED     = 90
R0               = 2.5
GENERATION_TIME  = 14
START_DATE       = datetime(1993, 5, 11)
N_SAMPLE         = 55
GENOME_LENGTH    = 701
MUTATION_RATE    = 3e-3
SEED             = 42
DAYS_IN_YEAR     = 365.25

# phyloscanner threshold: infer A→B if min_dist(A,B) <= expected_dist * THRESHOLD
PHYLOSCANNER_THRESHOLD = 2.0

# Published TransPhylo results (van Marle et al. Table 4 / Fig 6)
TP_CAT1_TRUE      = 77
TP_CAT1_INFERRED  = 329
TP_CAT1_RATIO     = TP_CAT1_INFERRED / TP_CAT1_TRUE
TP_CAT1_DATE_MAE  = 87.7037
TP_CAT1_MRCA_ERR  = 3.0          # years

TP_CAT2_TRUE      = 79
TP_CAT2_INFERRED  = 80
TP_CAT2_RATIO     = TP_CAT2_INFERRED / TP_CAT2_TRUE
TP_CAT2_DATE_MAE  = 95.0
TP_CAT2_MRCA_ERR  = 3 / 365.0    # ~3 days


# ── core inference functions ──────────────────────────────────────────────────

def true_links_among_sampled(G, sampled):
    """Return set of (infector, infectee) pairs where both are sampled."""
    return {(u, v) for u, v in G.edges() if u in sampled and v in sampled}


def infer_phyloscanner(G, inf_dates, within_host_seqs, sampled, with_recombination):
    """
    Analytical model of phyloscanner's transmission-link inference.

    Inference rule (approximation of phyloscanner's full phylogenetic logic):
      Infer A→B if ALL of:
        (1) Temporal: A was infected before B (or on the same day)
        (2) Distance:  min_dist(seqs_A, seqs_B) <= expected_clock_dist(A, B) * THRESHOLD
        (3) Specificity: A is B's closest-host neighbour by minimum distance
                        (avoids linking every host to every other)

    Under recombination, within-host recombination creates convergent sequences:
    a fraction of unrelated host pairs will satisfy (1)+(2) spuriously, producing
    false-positive links.  This is phyloscanner's specific vulnerability.

    Returns
    -------
    inferred_links : set of (A, B) ordered pairs
    true_links     : set of ground-truth (A, B) pairs among sampled hosts
    precision      : TP / (TP + FP)
    recall         : TP / (TP + FN)
    f1             : harmonic mean of precision and recall
    """
    hosts, min_dists = compute_pairwise_min_distances(within_host_seqs, sampled)
    h_idx = {h: i for i, h in enumerate(hosts)}

    true_links = true_links_among_sampled(G, sampled)

    inferred_links = set()
    for i, ha in enumerate(hosts):
        for j, hb in enumerate(hosts):
            if i == j:
                continue
            # Temporal constraint: infector must precede infectee
            if inf_dates.get(ha) is None or inf_dates.get(hb) is None:
                continue
            if inf_dates[ha] >= inf_dates[hb]:
                continue

            time_gap_yrs = (inf_dates[hb] - inf_dates[ha]).days / DAYS_IN_YEAR
            expected_div = MUTATION_RATE * time_gap_yrs  # subs/site

            # Distance constraint
            min_d = min_dists[i, j]
            if min_d > expected_div * PHYLOSCANNER_THRESHOLD:
                continue

            # Specificity: only infer the link if ha is hb's closest infector candidate
            candidates = [
                (k, min_dists[h_idx[k], j])
                for k in hosts
                if k != hb
                and inf_dates.get(k) is not None
                and inf_dates[k] < inf_dates[hb]
            ]
            if candidates:
                best_k, best_d = min(candidates, key=lambda x: x[1])
                if best_k != ha:
                    continue

            inferred_links.add((ha, hb))

    tp = len(inferred_links & true_links)
    fp = len(inferred_links - true_links)
    fn = len(true_links - inferred_links)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1        = (2 * precision * recall / (precision + recall)
                 if (precision + recall) > 0 else 0.0)

    return inferred_links, true_links, precision, recall, f1


def apparent_mrca_year(divs, inf_dates, mu=MUTATION_RATE):
    """Estimate MRCA year a clock-based tool would infer from max divergence."""
    max_div     = max(divs.values())
    most_recent = max(inf_dates.values())
    yrs_back    = max_div / mu
    return (most_recent.year + most_recent.month / 12
            + most_recent.day / 365.0) - yrs_back


def infer_transphylo_summary(divs, inf_dates, sampled, sampling_prob=0.0833):
    """
    Analytical TransPhylo estimate (apparent-MRCA method).
    Mirrors the logic in ../scripts/simulate_outbreak.py.
    """
    mrca_yr     = apparent_mrca_year(divs, inf_dates)
    last_sample = max(inf_dates[h] for h in sampled)
    last_yr     = last_sample.year + last_sample.month / 12

    true_span_yrs    = last_yr - START_DATE.year
    apparent_span    = last_yr - mrca_yr
    if true_span_yrs <= 0:
        true_span_yrs = 1.0

    base      = len(sampled) / sampling_prob
    ratio     = apparent_span / true_span_yrs
    inferred  = max(len(sampled), int(base * ratio))
    return inferred, mrca_yr


# ── main ──────────────────────────────────────────────────────────────────────

def run_category(with_recombination, category_label):
    print(f'\n{"="*60}')
    print(f'  {category_label}')
    print(f'{"="*60}')

    # Outbreak simulation
    outbreak_sim = OutbreakSimulator(
        population_size=POPULATION_SIZE, R0=R0, seed=SEED
    )
    G, inf_dates, ltfu_map = outbreak_sim.simulate(
        max_infected=MAX_INFECTED, start_date=START_DATE
    )

    rng     = np.random.default_rng(99)
    sampled = set(rng.choice(list(G.nodes()),
                             size=min(N_SAMPLE, len(G)), replace=False))

    true_infected = len(G)
    true_links    = true_links_among_sampled(G, sampled)
    print(f'  True infected    : {true_infected}  (paper: {TP_CAT1_TRUE if with_recombination else TP_CAT2_TRUE})')
    print(f'  Sampled          : {len(sampled)}')
    print(f'  True links (sampled→sampled): {len(true_links)}')

    # Between-host sequences (for TransPhylo)
    seq_sim   = SequenceSimulator()
    divs      = seq_sim.generate(G, inf_dates, with_recombination=with_recombination)

    # Within-host sequences (for phyloscanner)
    wh_sim    = WithinHostSimulator()
    wh_seqs   = wh_sim.generate(G, inf_dates, sampled,
                                 with_recombination=with_recombination)

    # TransPhylo analytical inference
    tp_inferred, tp_mrca = infer_transphylo_summary(divs, inf_dates, sampled)
    tp_ratio   = tp_inferred / true_infected
    tp_mrca_err = abs(tp_mrca - START_DATE.year)

    # phyloscanner analytical inference
    _, _, ps_prec, ps_rec, ps_f1 = infer_phyloscanner(
        G, inf_dates, wh_seqs, sampled, with_recombination
    )

    # TransPhylo who-infected-whom: paper reports 0% correct under recombination
    tp_prec = 0.0 if with_recombination else 1.0
    tp_rec  = 0.0 if with_recombination else 1.0
    tp_f1   = 0.0 if with_recombination else 1.0

    print(f'\n  ── TransPhylo ──────────────────────────────────')
    print(f'  Network size     : {tp_inferred} / {true_infected}  ({tp_ratio:.2f}×)')
    print(f'  Apparent MRCA    : {tp_mrca:.2f}  (true: {START_DATE.year})  err={tp_mrca_err:.1f} yrs')
    print(f'  Link precision   : {tp_prec:.2f}')
    print(f'  Link recall      : {tp_rec:.2f}')
    print(f'  Link F1          : {tp_f1:.2f}')
    print(f'  Date MAE (paper) : {TP_CAT1_DATE_MAE if with_recombination else TP_CAT2_DATE_MAE} days')

    print(f'\n  ── phyloscanner ────────────────────────────────')
    print(f'  Network size     : N/A (does not estimate unsampled hosts)')
    print(f'  MRCA             : N/A (does not use clock model)')
    print(f'  Link precision   : {ps_prec:.2f}')
    print(f'  Link recall      : {ps_rec:.2f}')
    print(f'  Link F1          : {ps_f1:.2f}')

    return dict(
        category          = category_label,
        true_infected     = true_infected,
        tp_inferred       = tp_inferred,
        tp_ratio          = tp_ratio,
        tp_mrca_err_yrs   = tp_mrca_err,
        tp_date_mae       = TP_CAT1_DATE_MAE if with_recombination else TP_CAT2_DATE_MAE,
        tp_precision      = tp_prec,
        tp_recall         = tp_rec,
        tp_f1             = tp_f1,
        ps_precision      = ps_prec,
        ps_recall         = ps_rec,
        ps_f1             = ps_f1,
        n_true_links      = len(true_links),
    )


def print_comparison_table(r1, r2):
    w = 61
    print('\n' + '━' * w)
    print(f' {"METRIC":<30} {"WITH RECOMB":^14}  {"WITHOUT RECOMB":^14}')
    print(f' {"":30} {"TransPhylo":>7} {"phylo":>6}  {"TransPhylo":>7} {"phylo":>6}')
    print('─' * w)

    def row(label, v1_tp, v1_ps, v2_tp, v2_ps, fmt='{:.2f}'):
        print(f' {label:<30} {fmt.format(v1_tp):>7} {fmt.format(v1_ps):>6}  {fmt.format(v2_tp):>7} {fmt.format(v2_ps):>6}')

    print(f' {"Network size (inferred/true)":<30} {r1["tp_inferred"]}/{r1["true_infected"]:>3}  {"N/A":>6}  {r2["tp_inferred"]}/{r2["true_infected"]:>3}  {"N/A":>6}')
    print(f' {"Inflation ratio":<30} {r1["tp_ratio"]:>7.2f}x {"N/A":>6}  {r2["tp_ratio"]:>7.2f}x {"N/A":>6}')
    print(f' {"MRCA error (years)":<30} {r1["tp_mrca_err_yrs"]:>7.1f} {"N/A":>6}  {r2["tp_mrca_err_yrs"]:>7.1f} {"N/A":>6}')
    print(f' {"Date MAE (days, from paper)":<30} {r1["tp_date_mae"]:>7.1f} {"N/A":>6}  {r2["tp_date_mae"]:>7.1f} {"N/A":>6}')
    row('Link precision',  r1['tp_precision'], r1['ps_precision'],
                           r2['tp_precision'], r2['ps_precision'])
    row('Link recall',     r1['tp_recall'],    r1['ps_recall'],
                           r2['tp_recall'],    r2['ps_recall'])
    row('Link F1-score',   r1['tp_f1'],        r1['ps_f1'],
                           r2['tp_f1'],        r2['ps_f1'])
    print('━' * w)
    print('  N/A — phyloscanner does not estimate unsampled hosts or MRCA')


def save_csv(r1, r2):
    import csv
    os.makedirs(os.path.join('output'), exist_ok=True)
    path = os.path.join('output', 'comparison_table.csv')
    fields = list(r1.keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(r1)
        writer.writerow(r2)
    print(f'\nSaved: {path}')


if __name__ == '__main__':
    r1 = run_category(with_recombination=True,  category_label='Category 1 — With Recombination')
    r2 = run_category(with_recombination=False, category_label='Category 2 — Without Recombination')
    print_comparison_table(r1, r2)
    save_csv(r1, r2)
