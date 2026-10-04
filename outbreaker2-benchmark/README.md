# outbreaker2 Benchmark — Cross-Tool Validation of Apollo Findings

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IvanR2625/apollo-transphylo-benchmark/blob/main/outbreaker2-benchmark/colab-notebook/outbreaker2_benchmark.ipynb)

**Author:** Ivan Raizada — Webber Academy, Calgary, AB

---

## What this is

Van Marle et al. (2025) showed that **TransPhylo over-estimates outbreak size 4.27× when genomes harbour recombination**, using the Apollo simulator for ground truth. They tested only TransPhylo.

This project extends the benchmark to **outbreaker2** — a fundamentally different tool that uses pairwise SNP distances rather than phylogenetic MRCA timing — and asks: *does the failure generalise?*

**Answer:** Yes. outbreaker2 also over-estimates under recombination, but less severely (~1.8×).

| | With Recombination | Without Recombination |
|---|---|---|
| **True infected** | 77 | 79 |
| **TransPhylo** (van Marle 2025) | **329 (4.27×)** | 80 (1.01×) |
| **outbreaker2** (this study) | **~140 (1.82×)** | ~83 (1.05×) |

The mechanisms differ: TransPhylo fails globally (MRCA timing), outbreaker2 fails locally (SNP distances). The global mechanism amplifies errors multiplicatively; the local one accumulates them additively.

---

## Why this is novel

No published study has benchmarked outbreaker2 against Apollo ground truth. This is the first cross-tool comparison under identical recombination conditions with a known true transmission network.

See the full research paper: [`paper/manuscript.md`](paper/manuscript.md)

---

## Project structure

```
outbreaker2-benchmark/
├── scripts/
│   ├── run_outbreaker2_analytical.py  — analytical benchmark (no R required)
│   └── figures.py                     — generate all three paper figures
├── colab-notebook/
│   └── outbreaker2_benchmark.ipynb    — self-contained Colab notebook
├── paper/
│   └── manuscript.md                  — complete research paper manuscript
├── output/
│   ├── comparison_table.csv           — results table (generated at runtime)
│   └── figures/                       — PNG figures (generated at runtime)
└── README.md
```

---

## How to run

### Option A — Google Colab (recommended)

Click the badge above, or open:  
`outbreaker2-benchmark/colab-notebook/outbreaker2_benchmark.ipynb`

No GPU required. All sections run on CPU. Section 10 (real outbreaker2 MCMC) is optional and takes ~10 min.

### Option B — Local Python

```bash
# From the apollo/ directory
pip install networkx matplotlib numpy scipy

# Run the analytical benchmark
python outbreaker2-benchmark/scripts/run_outbreaker2_analytical.py

# Generate all figures
python outbreaker2-benchmark/scripts/figures.py
```

Outputs go to `outbreaker2-benchmark/output/`.

---

## Key mechanism comparison

```
TransPhylo (4.27×)                    outbreaker2 (~1.82×)
──────────────────────────────────    ──────────────────────────────────
Recombination events                  Recombination events
       ↓                                     ↓
Max divergence inflated               Pairwise SNP distances inflated
       ↓                                     ↓
BEAST2: MRCA ~3 years early           t̂ = d/(μL) overestimates steps
       ↓                                     ↓
TransPhylo fills extra time           Extra intermediates per link
with phantom hosts                    accumulate additively
       ↓                                     ↓
GLOBAL amplification (severe)         LOCAL accumulation (moderate)
```

---

## Related project

This benchmark extends the parent project:  
[`../` — TransPhylo benchmark (van Marle et al. reproduction)](../README.md)

Parent Colab:  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IvanR2625/apollo-transphylo-benchmark/blob/main/colab-notebook/transphylo_benchmark_demo.ipynb)

---

## References

- van Marle G et al. (2025). Apollo. *Nature Communications* 16, 5783. https://doi.org/10.1038/s41467-025-60988-8
- Campbell F et al. (2019). outbreaker2. *PLOS Computational Biology* 15(3):e1006930.
- Didelot X et al. (2017). TransPhylo. *Mol Biol Evol* 34(4):997–1007.
