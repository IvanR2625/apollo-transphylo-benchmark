"""
simulate_within_host.py
-----------------------
Extends the Apollo outbreak simulation to generate within-host viral sequence
diversity per sampled host.

TransPhylo and outbreaker2 discard within-host diversity, using only a single
consensus sequence per host. phyloscanner uses all within-host sequences.
This module generates the within-host sequence populations that phyloscanner
needs — the same data Apollo produces at full scale on GPU.

Outputs
-------
WithinHostSimulator.generate() → {host_id: [div_1, ..., div_N]}
    N_WITHIN_HOST divergence values (from root) per host

compute_pairwise_min_distances()  → 2D array of min pairwise distances
compute_pairwise_mean_distances() → 2D array of mean pairwise distances
"""

import os
import sys
import numpy as np
from datetime import datetime

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
from simulate_outbreak import OutbreakSimulator, SequenceSimulator

# ── parameters ────────────────────────────────────────────────────────────────
GENOME_LENGTH    = 701
MUTATION_RATE    = 3e-3        # subs/site/year — between-host clock
INTRA_HOST_MU   = 5e-3        # subs/site/year — within-host (slightly elevated)
N_WITHIN_HOST   = 20          # sequences sampled per host
DAYS_IN_YEAR    = 365.25
START_DATE      = datetime(1993, 5, 11)
SEED            = 42


class WithinHostSimulator:
    """
    For each sampled host, generates a population of N_WITHIN_HOST viral
    sequences that reflect within-host evolution since infection.

    Biology modelled:
      1. Founder divergence   — sequences inherited from infector carry the
                                 between-host clock divergence
      2. Within-host mutation — additional mutations accumulate during the
                                 within-host replication period
      3. Within-host recombination (optional) — creates mosaic sequences with
                                 elevated divergence relative to the founder

    This is the data Apollo generates at scale (millions of genomes per host)
    and that phyloscanner exploits.  TransPhylo uses only the per-host median.
    """

    def __init__(self, mu_between=MUTATION_RATE, mu_within=INTRA_HOST_MU,
                 L=GENOME_LENGTH, n_seqs=N_WITHIN_HOST, seed=SEED):
        self.mu_b  = mu_between
        self.mu_w  = mu_within
        self.L     = L
        self.n     = n_seqs
        self.rng   = np.random.default_rng(seed)

    def generate(self, G, inf_dates, sampled_hosts,
                 sample_date_offset_days=180, with_recombination=True):
        """
        Parameters
        ----------
        G               : nx.DiGraph transmission network
        inf_dates       : {host_id: datetime} infection date per host
        sampled_hosts   : set of sampled host ids
        sample_date_offset_days : days after infection when sample is taken
        with_recombination : whether to apply within-host recombination

        Returns
        -------
        {host_id: np.ndarray of shape (N_WITHIN_HOST,)}
            Each value is total divergence from root sequence.
        """
        import networkx as nx
        root     = next(h for h in G.nodes() if G.in_degree(h) == 0)
        root_dt  = inf_dates[root]

        within_host_seqs = {}

        for host in sampled_hosts:
            if host not in inf_dates:
                continue

            # Between-host clock divergence (founder sequence)
            elapsed_days  = max(0, (inf_dates[host] - root_dt).days)
            elapsed_years = elapsed_days / DAYS_IN_YEAR
            founder_div   = self.rng.poisson(
                max(0, self.mu_b * self.L * elapsed_years)
            ) / self.L

            # Within-host diversity: each of the N sequences evolves from founder
            within_days  = sample_date_offset_days
            within_years = within_days / DAYS_IN_YEAR
            expected_wh_muts = self.mu_w * self.L * within_years

            seqs = np.zeros(self.n)
            for i in range(self.n):
                intra_muts = self.rng.poisson(expected_wh_muts)
                intra_div  = intra_muts / self.L

                recomb_div = 0.0
                if with_recombination:
                    # Within-host recombination creates mosaic sequences.
                    # Each sequence independently has a chance to carry a
                    # recombinant segment from an older viral lineage.
                    if self.rng.random() < 0.25:
                        extra_yrs  = self.rng.uniform(0.5, 3.0)
                        recomb_div = self.rng.exponential(self.mu_b * extra_yrs)

                seqs[i] = founder_div + intra_div + recomb_div

            within_host_seqs[host] = seqs

        return within_host_seqs


def compute_pairwise_min_distances(within_host_seqs, sampled_hosts):
    """
    Minimum pairwise divergence distance between every ordered pair of hosts.

    phyloscanner uses minimum (not mean) distance to avoid recombinant
    outliers biasing the inference.

    Returns
    -------
    hosts  : list of host IDs (row/column index)
    dists  : 2D numpy array, dists[i, j] = min distance from host i to host j
    """
    hosts = sorted(h for h in sampled_hosts if h in within_host_seqs)
    n     = len(hosts)
    dists = np.zeros((n, n))

    for i, ha in enumerate(hosts):
        for j, hb in enumerate(hosts):
            if i == j:
                continue
            seqs_a = within_host_seqs[ha]
            seqs_b = within_host_seqs[hb]
            # Min absolute divergence difference across all cross-host pairs
            diff   = np.abs(seqs_a[:, None] - seqs_b[None, :])
            dists[i, j] = diff.min()

    return hosts, dists


def compute_pairwise_mean_distances(within_host_seqs, sampled_hosts):
    """
    Mean pairwise divergence distance — equivalent to consensus-based tools
    such as TransPhylo (which use only the median/consensus sequence).

    Returns
    -------
    hosts  : list of host IDs
    dists  : 2D numpy array of mean distances
    """
    hosts = sorted(h for h in sampled_hosts if h in within_host_seqs)
    n     = len(hosts)
    dists = np.zeros((n, n))

    for i, ha in enumerate(hosts):
        for j, hb in enumerate(hosts):
            if i == j:
                continue
            ma = within_host_seqs[ha].mean()
            mb = within_host_seqs[hb].mean()
            dists[i, j] = abs(ma - mb)

    return hosts, dists


# ── standalone test ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    import networkx as nx

    sim      = OutbreakSimulator(population_size=300, R0=2.5, seed=SEED)
    G, inf_dates, ltfu_map = sim.simulate(max_infected=90, start_date=START_DATE)

    rng      = np.random.default_rng(99)
    sampled  = set(rng.choice(list(G.nodes()), size=min(55, len(G)), replace=False))

    wh_sim   = WithinHostSimulator()

    seqs_with = wh_sim.generate(G, inf_dates, sampled, with_recombination=True)
    seqs_no   = wh_sim.generate(G, inf_dates, sampled, with_recombination=False)

    print(f'Hosts simulated : {len(seqs_with)}')
    print(f'Sequences/host  : {N_WITHIN_HOST}')
    print()

    h0     = list(seqs_with.keys())[0]
    d_with = seqs_with[h0]
    d_no   = seqs_no[h0]
    print(f'Example host {h0}:')
    print(f'  With recombination    — mean: {d_with.mean():.4f}  std: {d_with.std():.4f}  max: {d_with.max():.4f}')
    print(f'  Without recombination — mean: {d_no.mean():.4f}  std: {d_no.std():.4f}  max: {d_no.max():.4f}')

    hosts, min_d = compute_pairwise_min_distances(seqs_with, sampled)
    hosts, mea_d = compute_pairwise_mean_distances(seqs_with, sampled)
    print(f'\nPairwise distances (with recombination, first 3×3):')
    print('  min:\n', min_d[:3, :3].round(4))
    print('  mean:\n', mea_d[:3, :3].round(4))
