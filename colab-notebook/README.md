# Apollo TransPhylo Benchmark — Google Colab Notebook

## What this is

`transphylo_benchmark_demo.ipynb` is a self-contained Google Colab notebook that backs the claim:

> *"329 predicted infections against a true 77, and 87.7 days mean error on infection dates.  
> I hadn't appreciated that tools used on real outbreak data could be that far off,  
> or that you'd need a simulator to find out."*

It does this in five progressive sections, from pure paper-number visualisation (no simulation needed) through to running the real Apollo GPU simulator.

---

## How to open it in Google Colab

**Option A — upload directly**

1. Go to [colab.research.google.com](https://colab.research.google.com)
2. File → Upload notebook → select `transphylo_benchmark_demo.ipynb`
3. Runtime → Run all  *(sections 1–4 work without GPU)*

**Option B — open from GitHub**

1. Go to [colab.research.google.com](https://colab.research.google.com)
2. File → Open notebook → GitHub tab
3. Paste your repo URL and select `apollo/colab-notebook/transphylo_benchmark_demo.ipynb`

**For Part 5 (real Apollo):** enable GPU first:  
Runtime → Change runtime type → Hardware accelerator → T4 GPU

---

## Sections at a glance

| Section | GPU needed? | What you get |
|---------|-------------|--------------|
| 1. Paper findings | No | Exact paper numbers visualised (329 vs 77, MRCA, date errors) |
| 2. Mechanism | No | Arrow-diagram of the recombination → phantom hosts causal chain |
| 3. Mini simulation | No | CPU outbreak simulator matching the paper's structure |
| 4. Inference failure | No | Numerical demonstration that recombination breaks MRCA estimation |
| 5. Real Apollo | Yes (T4) | Install CATE, configure Apollo, run the full simulation |

---

## Dependencies installed automatically

The first cell runs:
```
pip install biopython networkx matplotlib numpy scipy
```

No manual setup required for sections 1–4.

---

## Expected outputs

Running sections 1–4 produces four PNG files in the Colab working directory:

| File | What it shows |
|------|---------------|
| `paper_findings.png` | 3-panel figure: infected count, MRCA year, date error histogram |
| `mechanism.png` | Causal chain diagram |
| `ground_truth_network.png` | Transmission network coloured by LTFU host type |
| `divergence_comparison.png` | Clock-like vs recombination-inflated divergence scatter |
| `simulation_vs_paper.png` | Our simulation results vs paper's TransPhylo numbers |

---

## Paper reference

van Marle et al. (2025). *Apollo: a comprehensive GPU-powered within-host simulator for viral evolution.*  
**Nature Communications** 16, 5783.  
https://doi.org/10.1038/s41467-025-60988-8
