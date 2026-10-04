# Apollo TransPhylo Benchmark

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IvanR2625/apollo-transphylo-benchmark/blob/main/colab-notebook/transphylo_benchmark_demo.ipynb)

A project that numerically backs the statement:

> *"The part that struck me was using simulated genomes with known ground truth to test TransPhylo — 329 predicted infections against a true 77, and 87.7 days mean error on infection dates. I hadn't appreciated that tools used on real outbreak data could be that far off, or that you'd need a simulator to find out."*

**Paper:** van Marle et al. (2025). Apollo: a comprehensive GPU-powered within-host simulator for viral evolution. *Nature Communications* 16, 5783. https://doi.org/10.1038/s41467-025-60988-8

**Authors:** Ivan Raizada  
**Affiliation:** Webber Academy, Calgary, AB

---

## The core finding

Apollo generates viral genomes for a simulated outbreak with a **known ground truth** (exactly who infected whom, and when). TransPhylo then infers the transmission network from those sequences alone.

| | With Recombination | Without Recombination |
|---|---|---|
| **True infected** | 77 | 79 |
| **TransPhylo inferred** | **329** | 80 |
| **Overestimation** | **4.27×** | ~1× |
| **MRCA error** | ~3 years | ~3 days |
| **Infection date MAE** | **87.7 days** | 95 days |

The single variable that explains everything: **recombination**.

---

## Why recombination breaks TransPhylo

```
Recombination events (14 hotspots)
        ↓
Sequences more divergent than clock predicts
        ↓
BEAST2 infers MRCA ~1990  (true: 1993, 3-year error)
        ↓
TransPhylo fills 3 extra years with phantom intermediate hosts
        ↓
329 inferred vs 77 true  (4.27×)
```

Without recombination, sequences accumulate mutations at a clock-like rate. BEAST2 infers the correct MRCA (~1993). TransPhylo then estimates 80 individuals — nearly perfect.

---

## Project structure

```
apollo/
├── scripts/
│   ├── simulate_outbreak.py       — outbreak simulator + divergence generator
│   └── visualize_paper_findings.py — exact paper numbers → publication figures
├── colab-notebook/
│   ├── transphylo_benchmark_demo.ipynb  — interactive Colab notebook
│   └── README.md                        — Colab instructions
├── data/                          — simulation outputs (generated at runtime)
└── output/
    └── figures/                   — PNG figures
```

---

## Quick start

### Option A — Google Colab (recommended, no local install)

Open `colab-notebook/transphylo_benchmark_demo.ipynb` in Google Colab.  
Sections 1–4 run without GPU. Section 5 (real Apollo) needs T4 GPU.

See `colab-notebook/README.md` for step-by-step instructions.

### Option B — Local Python

```bash
pip install networkx matplotlib numpy scipy biopython
python scripts/visualize_paper_findings.py   # paper figures
python scripts/simulate_outbreak.py          # simulation + comparison
```

Outputs go to `output/figures/`.

---

## Scripts

### `scripts/visualize_paper_findings.py`

Uses the **exact numbers published in the paper** (no simulation). Produces three figures:

| Output | What it shows |
|--------|---------------|
| `output/figures/paper_comparison_fig6.png` | 329 vs 77, MRCA year error |
| `output/figures/error_distribution.png` | Infection date error histograms |
| `output/figures/mechanism_diagram.png` | Causal chain: recombination → phantom hosts |

### `scripts/simulate_outbreak.py`

Implements the statistical skeleton of Apollo's experiment:

- **`OutbreakSimulator`** — stochastic SIR with three LTFU host types (non-LTFU 50%, complete-LTFU 25%, partial-LTFU 25%), explicit who-infected-whom tracking
- **`SequenceSimulator`** — generates divergence values with/without recombination (Poisson clock mutations + recombination hotspot events)
- **`apparent_mrca_year()`** — estimates the MRCA year a clock-calibrated tool would infer
- **`count_inferred_hosts()`** — estimates TransPhylo's inferred population using the paper's sampling probability (0.0833)

The simulation reproduces the qualitative finding: recombination inflates inferred population 4×; without recombination inference is nearly exact.

---

## Running real Apollo

Apollo is GPU-accelerated and part of the CATE package (MIT license).

| Resource | Link |
|----------|------|
| GitHub | https://github.com/theLongLab/CATE |
| Anaconda | `conda install -c deshan_CATE cate` |
| Wiki | https://github.com/theLongLab/CATE/wiki/Apollo |
| User manual | https://github.com/theLongLab/CATE/tree/main/Apollo_User_Manual |
| Video tutorial | https://vimeo.com/1079571253/571b3b6fd3 |
| Zenodo DOI | https://doi.org/10.5281/zenodo.7987768 |

The Colab notebook (Part 5) installs CATE automatically when a GPU is available and provides a configuration file matching the paper's Category 1 experiment.

---

## Apollo simulation parameters (from paper, Table 1)

| Parameter | Value |
|-----------|-------|
| Population size | 300 |
| Genome length | 701 bp |
| Mutation rate | 3 × 10⁻³ subs/site/year |
| Recombination hotspots | 14 |
| Mutation hotspots | 30 |
| Start date | 1993-05-11 |
| Sampled sequences | 55 (Cat. 1), 56 (Cat. 2) |
| LTFU: non / complete / partial | 50% / 25% / 25% |

**TransPhylo parameters (Table 4):**

| Parameter | Value |
|-----------|-------|
| Infection rate (shape, scale) | 1, 0.99995 |
| Sampling rate (shape, scale) | 1, 0.5 |
| MCMC iterations | 100,000 |
| Sampling probability | 0.0833 |

---

## Why a simulator is the only way to find this

On real outbreak data you never know the true transmission network. Any evaluation of TransPhylo would be circular — you'd use the very tool you're trying to assess.

Apollo provides **known ground truth**: the exact sequence of infections, their dates, and the sequences that resulted. This is the only way to measure how far off a genomic epidemiology tool is — and in this case the answer is 4.27× with a 3-year MRCA error.

---

## Dependencies

| Package | Use |
|---------|-----|
| NumPy | Arrays, random sampling |
| Matplotlib | All figures |
| NetworkX | Transmission network graph |
| SciPy | Linear regression |
| Biopython | (used if extending to real sequences) |

---

## License

MIT — use and adapt freely. Paper data and figures are from van Marle et al. (2025) under the terms of Nature Communications open access.
