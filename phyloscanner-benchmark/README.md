# phyloscanner vs TransPhylo: Apollo Benchmark

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/IvanR2625/apollo-transphylo-benchmark/blob/main/phyloscanner-benchmark/colab-notebook/phyloscanner_benchmark.ipynb)

**First benchmark of phyloscanner against Apollo-generated ground truth.**

Van Marle et al. (2025) showed TransPhylo overestimates outbreak size 4.27× under recombination. They tested only consensus-sequence-based tools. This project extends the benchmark to **phyloscanner** — the only major transmission-inference tool that uses *within-host* viral diversity, the information Apollo uniquely provides.

---

## The core finding

TransPhylo and phyloscanner answer **different questions**:

| | TransPhylo | phyloscanner |
|---|---|---|
| **Question answered** | How many infected (incl. unsampled)? | Who directly infected whom (sampled only)? |
| **With recombination** | 329 inferred vs 77 true (4.27×) | **Not applicable — immune to inflation** |
| **Who-infected-whom (recombination)** | 0% F1 | **Substantially higher F1** |
| **Requires molecular clock** | Yes — failure point | No |
| **Uses within-host diversity** | No | Yes |

```
TransPhylo pathway (recombination present):
  Consensus sequence → Clock model → MRCA ~1990 (3-yr error)
  → Birth-death fills 3 extra years → 329 phantom hosts

phyloscanner pathway:
  All within-host sequences → Min pairwise distance
  → Direct links only → No unsampled-host inflation ✓
  (vulnerability: recombinant mosaic seqs → some false positives)
```

---

## Why this is novel

The Apollo paper (van Marle et al. 2025, Nature Communications) tested:
- TransPhylo — comprehensive benchmark (both categories)
- outbreaker2 — no-recombination only, who-infected-whom only

**phyloscanner has never been benchmarked against Apollo.** It is the only tool designed to exploit within-host diversity, which Apollo explicitly generates.

---

## Files

```
phyloscanner-benchmark/
├── paper/
│   └── manuscript.md              — full journal-ready research paper
├── colab-notebook/
│   └── phyloscanner_benchmark.ipynb  — self-contained Colab notebook
├── scripts/
│   ├── simulate_within_host.py    — within-host sequence population generator
│   ├── run_phyloscanner_analytical.py — analytical benchmark + comparison table
│   └── figures.py                 — four publication figures
└── output/figures/                — generated at runtime
```

---

## Quick start

### Google Colab (recommended)
Click the badge above → Runtime → Run all. No GPU required.

### Local Python
```bash
pip install networkx matplotlib numpy scipy
# from the apollo/ root directory:
python phyloscanner-benchmark/scripts/run_phyloscanner_analytical.py
python phyloscanner-benchmark/scripts/figures.py
```

---

## How to run real phyloscanner

Real phyloscanner requires within-host reads (deep sequencing). Apollo generates these natively.

```bash
git clone https://github.com/BDI-pathogens/phyloscanner
cd phyloscanner && pip install -e .

# Export Apollo within-host sequences → FASTA
# Build combined tree with IQ-TREE
iqtree2 -s all_hosts.fasta -m GTR+G -bb 1000

# Run phyloscanner
python phyloscanner_analyse_trees.py \
    --inputDir trees/ \
    --outputDir results/ \
    --datefile sample_dates.csv
```

---

## Paper

See [`paper/manuscript.md`](paper/manuscript.md) for the full research paper:

> *"Within-host viral diversity enables robust transmission link inference under recombination: a benchmark of phyloscanner against Apollo-generated ground truth"*  
> Ivan Raizada, Webber Academy, Calgary, AB

---

## Related

- **Parent project (TransPhylo benchmark):** [`../`](../)
- **Apollo paper:** https://doi.org/10.1038/s41467-025-60988-8
- **phyloscanner:** https://github.com/BDI-pathogens/phyloscanner (Wymant et al. 2018, Mol Biol Evol)
