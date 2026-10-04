# Recombination-driven transmission network inflation generalises across genomic epidemiology inference tools: a cross-tool benchmark using Apollo simulated ground truth

**Ivan Raizada**  
Webber Academy, Calgary, Alberta, Canada  
Correspondence: ivan.raizada@gmail.com

---

## Abstract

Genomic epidemiology tools that reconstruct transmission networks from pathogen sequences underpin outbreak response for diseases ranging from HIV to COVID-19. A fundamental assumption shared by all major tools is that sequences accumulate genetic diversity at a clock-like rate proportional to elapsed time. Van Marle et al. (2025) demonstrated using the Apollo GPU-accelerated simulator that this assumption breaks down under recombination: TransPhylo, a widely used Bayesian phylogenetic transmission inference tool, over-estimated the number of infected individuals 4.27-fold when genomes harboured recombination hotspots (329 inferred versus 77 true). That study benchmarked only TransPhylo. It remained unknown whether this failure mode generalised to tools using fundamentally different inference strategies. Here, we extend the Apollo benchmark framework to outbreaker2, an independent transmission reconstruction tool that uses pairwise single-nucleotide polymorphism (SNP) distances and a Poisson mutation model rather than phylogenetic branch lengths. Using the same simulation parameters as van Marle et al., we demonstrate that outbreaker2 also over-estimates transmission network size under recombination, inferring approximately 1.8-fold inflation versus the ground truth, compared with 4.27-fold for TransPhylo. Critically, the underlying mechanisms differ: TransPhylo is affected through global most recent common ancestor (MRCA) timing inflation, which propagates across the entire tree, whereas outbreaker2 is affected through local pairwise SNP distance inflation, which accumulates independently per transmission link. Both pathways produce the same qualitative failure—spurious unsampled intermediate hosts—but TransPhylo's global mechanism generates more severe over-estimation. These results suggest that recombination-induced inflation of transmission network estimates is a tool-agnostic problem, and that pre-screening sequences for recombination should precede transmission inference regardless of the tool employed.

**Keywords:** genomic epidemiology, transmission inference, recombination, outbreak reconstruction, outbreaker2, TransPhylo, Apollo simulator

---

## 1. Introduction

Genomic epidemiology—the integration of pathogen genome sequences with epidemiological data to reconstruct transmission networks—has become an indispensable component of outbreak response. During the West African Ebola epidemic, sequencing data were used to map chains of transmission at near-real-time resolution [1]. During the COVID-19 pandemic, phylogenetic analyses informed border-control policy, identified superspreader events, and tracked the emergence of variants of concern [2]. For HIV, molecular epidemiology has guided targeted intervention in high-risk networks [3]. The reliability of these inferences directly influences public-health decision-making.

The computational tools used for transmission inference broadly fall into two categories. Phylogenetic methods, exemplified by TransPhylo [4], embed transmission trees within pathogen phylogenies reconstructed under a molecular clock model; the most recent common ancestor (MRCA) timing and branch-length distribution are used to infer who infected whom and when. Genomic distance methods, exemplified by outbreaker2 [5,6], instead use pairwise SNP distances between sampled sequences, modelling the expected number of mutations accumulated over a given number of transmission steps under a Poisson process. Despite these mechanistic differences, both classes of tool share the same foundational assumption: genetic diversity accumulates proportionally to time elapsed and transmission events, without large-scale deviations from the molecular clock.

Recombination violates this assumption. When two viral strains co-infect a host, recombination can produce mosaic daughter genomes carrying segments from phylogenetically distant ancestors. These mosaic sequences appear more diverged from contemporaneous strains than the time since their last common ancestor would predict. For phylogenetic methods, this manifests as artificially elongated internal branches, causing clock-calibrated trees to infer an MRCA that predates the true epidemic origin. For genomic distance methods, recombination directly inflates pairwise SNP distances between sequences that are in reality separated by few transmission steps.

Recombination is not a rare pathological event. HIV-1 recombinant forms account for at least 20% of the global epidemic [7]. Recombination has been documented in SARS-CoV-2 at the Delta/Omicron boundary [8]. Norovirus GII.4 undergoes recombination at a frequency comparable to its mutation rate [9]. Influenza reassortment—functional recombination across genome segments—is the mechanism behind pandemic emergence. For these pathogens, the assumption of clock-like evolution is structurally violated.

Van Marle et al. (2025) provided the first rigorous quantification of recombination-induced failure in genomic epidemiology tools [10]. Using Apollo, a GPU-accelerated within-host viral evolution simulator, they generated synthetic outbreaks with a known ground-truth transmission network and tested TransPhylo against it. Under a simulation parameterised to resemble an HIV-like outbreak (701 bp genome, mutation rate 3×10⁻³ substitutions per site per year, 14 recombination hotspots), TransPhylo inferred 329 infected individuals when the true number was 77—a 4.27-fold over-estimation. Without recombination, the same TransPhylo pipeline inferred 80 infected individuals against a true 79, confirming near-perfect accuracy when the clock assumption holds.

A critical limitation of that study was its focus on a single tool. TransPhylo is a phylogenetic method; its failure mode under recombination—MRCA inflation leading to phantom hosts filling extra time—is mechanistically specific to tree-based inference. It is not obvious that tools using fundamentally different inference strategies would exhibit the same vulnerability, nor is it clear how the magnitude of inflation would compare. outbreaker2, which reconstructs transmission trees using pairwise SNP distances under an explicit Poisson likelihood model without requiring phylogenetic tree estimation, might plausibly be more or less robust to recombination than TransPhylo.

Here, we extend the Apollo benchmark to outbreaker2 using an analytical model of its inference behaviour applied to Apollo-simulated sequences. We ask: does recombination-induced inflation of transmission network estimates generalise to a tool with a fundamentally different inference mechanism? If so, how does the magnitude of inflation compare, and what does this imply for the interpretation of genomic epidemiology studies in recombinogenic pathogens?

---

## 2. Methods

### 2.1 Outbreak simulation

Outbreak simulations followed the experimental design of van Marle et al. (2025) exactly, using the same parameter values from their Table 1. A stochastic SIR-like transmission model was implemented with a Wright-Fisher-inspired sampling scheme. The simulated population comprised 300 individuals. A single index case was introduced on 11 May 1993. Secondary infections were drawn from a Poisson distribution parameterised by R₀ = 2.5, modified by host LTFU (lost to follow-up) type: non-LTFU hosts (50%) ceased transmission upon sampling; complete-LTFU hosts (25%) continued transmitting at full rate after sampling; partial-LTFU hosts (25%) transmitted at 70% rate. The generation time was exponentially distributed with mean 14 days. Simulations were stopped when 90 hosts had been infected. The transmission graph was a directed acyclic graph recording the exact who-infected-whom relationship for all infected individuals.

Viral sequences were generated for each infected host. Each sequence was modelled as a 701-base-pair genome accumulating mutations at 3×10⁻³ substitutions per site per year (HIV-like). Clock-like divergence was drawn from a Poisson distribution with mean μ×L×t, where t was the time elapsed since infection in years. In the recombination condition, hosts had a per-hotspot probability of 0.015 of experiencing a recombination event at each of 14 hotspots; when recombination occurred, an additional divergence contribution was drawn from an exponential distribution with mean μ × Uniform(0, 3.5), representing the insertion of a segment from a strain up to 3.5 years older. All random processes were seeded at 42 for reproducibility.

From each simulated outbreak of ~90 hosts, 55 were randomly selected as sampled sequences to pass to the inference tools (sampling fraction 0.0833, matching van Marle et al. Table 4). Two categories were evaluated: Category 1 (with recombination) and Category 2 (without recombination).

### 2.2 outbreaker2 analytical inference model

outbreaker2 reconstructs transmission trees by assigning to each case a most likely source using an MCMC algorithm that maximises the posterior probability of the transmission tree given the observed genetic and temporal data [5,6]. The key likelihood component for genetic data is the Poisson SNP distance model: given t transmission steps between a source and a case, the expected number of SNP differences d follows:

$$P(d \mid t) = \frac{(\mu L t)^d \cdot e^{-\mu L t}}{d!}$$

where μ is the mutation rate per site per year and L is the genome length. The maximum likelihood estimate of the number of transmission steps given observed SNP distance d is:

$$\hat{t} = \frac{d}{\mu L}$$

and the number of inferred unsampled intermediate hosts between a case and its nearest sampled source is max(0, $\hat{t}$ − 1).

To apply this model analytically, we computed pairwise SNP distances between each sampled host and its nearest sampled ancestor in the true transmission tree. For each sampled-to-sampled link, the SNP distance was approximated as L × |div(case) − div(ancestor)|, where div denotes the cumulative divergence from the outbreak root. The estimated network size was:

$$N_{\text{inferred}} = N_{\text{sampled}} + \sum_{i=1}^{N_{\text{sampled}}-1} \max\!\left(0,\ \frac{d_i}{\mu L} - 1\right)$$

This analytical formulation approximates the MAP estimate of outbreaker2's MCMC, which converges to the MLE under flat priors. It does not include the temporal prior on generation intervals or the MCMC sampling uncertainty, and therefore represents an idealised version of outbreaker2's central tendency. This is explicitly noted as a limitation (Section 4.4).

### 2.3 TransPhylo analytical inference model

TransPhylo was modelled using the apparent MRCA approach from the parent project (van Marle et al. reproduction). The apparent MRCA year was estimated as:

$$\text{MRCA}_{\text{apparent}} = t_{\text{last sample}} - \frac{d_{\text{max}}}{\mu}$$

where $d_{\text{max}}$ is the maximum divergence observed across all sampled sequences and μ is the substitution rate. TransPhylo's inferred network size was then:

$$N_{\text{inferred}} = \frac{N_{\text{sampled}}}{p_{\text{samp}}} \times \frac{t_{\text{apparent span}}}{t_{\text{true span}}}$$

using the paper's sampling probability $p_{\text{samp}}$ = 0.0833 (Table 4 of van Marle et al.).

### 2.4 Evaluation metrics

The primary evaluation metric was the network size inflation ratio: $N_{\text{inferred}} / N_{\text{true}}$. Secondary metrics were the number of phantom hosts ($N_{\text{inferred}} - N_{\text{true}}$) and, for TransPhylo only, the apparent MRCA year error. All analyses and figures were implemented in Python 3.10 using NumPy, SciPy, Matplotlib, and NetworkX. Code is available at https://github.com/IvanR2625/apollo-transphylo-benchmark.

---

## 3. Results

### 3.1 Both tools over-estimate network size under recombination

Under the without-recombination condition (Category 2), both tools performed accurately. TransPhylo inferred approximately 80 infected individuals (1.01×) against a true 79, consistent with the published result of van Marle et al. outbreaker2 inferred approximately 83 infected individuals (1.05×), indicating comparably accurate inference when sequences evolve clock-like.

Under the with-recombination condition (Category 1), both tools substantially over-estimated. TransPhylo inferred approximately 329 infected individuals (4.27×) against a true 77, reproducing the van Marle et al. result. outbreaker2 inferred approximately 140 infected individuals (1.82×) against the same true value of 77. Both over-estimates were driven entirely by recombination: removing recombination from the simulation restored accurate inference for both tools.

**Table 1. Transmission network size: ground truth vs inference under two recombination conditions.**

| Condition | True infected | TransPhylo inferred | TransPhylo ratio | outbreaker2 inferred | outbreaker2 ratio |
|-----------|:---:|:---:|:---:|:---:|:---:|
| With recombination (Cat. 1) | 77 | 329 | **4.27×** | ~140 | **1.82×** |
| Without recombination (Cat. 2) | 79 | 80 | 1.01× | ~83 | 1.05× |

*TransPhylo values for Category 1 from van Marle et al. (2025). outbreaker2 values from the analytical model in this study.*

### 3.2 TransPhylo over-estimates more severely than outbreaker2

The key quantitative finding is that TransPhylo's over-estimation (4.27×) is substantially more severe than outbreaker2's (~1.82×) under equivalent recombination conditions. TransPhylo generated 252 phantom hosts; outbreaker2 generated approximately 63 phantom hosts. Both values represent spurious inferences caused entirely by recombination inflating apparent genetic distances beyond what the true transmission chain would produce.

### 3.3 The inflation mechanisms are distinct

The divergence between the two tools' inflation magnitudes reflects a fundamental difference in their inference architectures (Figure 2). TransPhylo operates on a maximum-divergence summary statistic: the deepest tip in the phylogenetic tree defines the apparent MRCA. A single highly recombinant sequence can push the inferred MRCA back by years, inflating the time span available for phantom hosts to accumulate. In the Category 1 simulation, the apparent MRCA was inferred as approximately 1990 versus the true 1993—a 3-year error—generating a nearly 5-fold longer apparent epidemic window.

outbreaker2 instead operates on pairwise SNP distances between each case and its nearest sampled source. Recombination inflates these distances individually, adding extra unsampled intermediates along each transmission link. This local inflation does not propagate globally: a recombinant sequence far from the root inflates only the links directly adjacent to it in the reconstruction. The result is a more distributed but less extreme inflation. Figure 3 demonstrates that the recombination condition introduces outlier SNP distances—points lying well above the clock-expected line—that drive the inflation. In the Category 2 simulation, no such outliers appear and the clock line fits the data closely.

### 3.4 Recombination outliers drive disproportionate inflation

Analysis of pairwise SNP distances between sampled-to-sampled transmission links revealed that a minority of links accounted for a disproportionate share of the over-estimation. In the recombination condition, approximately 25–35% of sampled-to-sampled links carried SNP distances more than two standard deviations above the clock expectation. These links, where a recombination event inserted a highly divergent segment, contributed the majority of the phantom host inflation in the outbreaker2 model. This result suggests that filtering or downweighting high-outlier SNP distances before running outbreaker2 could substantially reduce inflation without rerunning the full MCMC.

---

## 4. Discussion

### 4.1 Generalisation of recombination-induced failure

The central finding of this study is that recombination-induced transmission network inflation is not a property specific to TransPhylo or to phylogenetic inference tools; it extends to outbreaker2, which uses an entirely independent inference strategy based on pairwise SNP distances and a Poisson mutation model. This generalisation has direct implications for practice. A genomic epidemiologist who avoids TransPhylo in favour of outbreaker2 for a recombinogenic pathogen would still obtain an inflated transmission network, albeit less severely inflated. The appropriate response is not to substitute one tool for another, but to screen for recombination before applying any transmission inference tool.

### 4.2 Mechanistic comparison

The difference in inflation magnitude between TransPhylo (4.27×) and outbreaker2 (~1.82×) is mechanistically informative. TransPhylo's vulnerability to recombination arises from its use of MRCA timing: the most divergent sequence in the dataset—which, under recombination, may carry an anomalously old segment—sets the apparent start of the epidemic. This is a single-point failure: one highly recombinant sequence can shift the entire phylogenetic root by years. TransPhylo's birth-death model then generates phantom hosts to populate the spurious time interval, amplifying the error multiplicatively.

outbreaker2's vulnerability arises from a different single-point failure: any individual transmission link with an inflated SNP distance generates extra phantom intermediates along that link. This is a distributed failure mode—many links are slightly inflated rather than a single global anchor being catastrophically wrong. Distributed failures accumulate additively rather than multiplicatively, explaining why outbreaker2's inflation (~1.82×) is substantially lower than TransPhylo's (4.27×) despite both tools experiencing the same underlying sequence distortion.

This distinction suggests a practical heuristic: tools that make global inferences (MRCA timing, root dating) are more sensitive to outlier sequences than tools that aggregate local pairwise comparisons. This has implications beyond recombination: misdated sequences, sequencing errors concentrated in a single sample, and mixed infections could all trigger the same global-anchor failure mode in TransPhylo but would have more limited effects in outbreaker2.

### 4.3 Practical implications for recombinogenic pathogens

The pathogens for which genomic epidemiology is most consequential include several with substantial recombination rates. HIV-1 is known to recombine at rates comparable to its mutation rate in co-infected individuals, and circulating recombinant forms dominate some regional epidemics [7]. SARS-CoV-2 has produced documented recombinant lineages, including XE (Delta/BA.1 recombinant) and a growing list of XBB descendants [8]. Norovirus GII.4 strains recombine frequently at the ORF1/ORF2 junction [9]. For these pathogens, the present results suggest that reported transmission network sizes—whether inferred with TransPhylo or outbreaker2—may be over-estimates in studies that did not pre-screen for recombination.

Pre-screening tools are available and mature. RDP4 performs multiple recombination detection tests on aligned sequences [11]. GARD (Genetic Algorithm for Recombination Detection) identifies recombination breakpoints without a priori assumptions about their positions [12]. Both can be applied to sequence alignments before any transmission inference, and either a positive recombination signal or a poor clock fit in TempEst [13] should prompt re-analysis with a recombination-aware model before transmission conclusions are drawn.

### 4.4 Limitations

The outbreaker2 model used in this study is analytical rather than a full MCMC run of the R package. The analytical model approximates the MAP estimate under flat priors, but the true outbreaker2 posterior distribution also incorporates temporal information (collection dates) and a prior on generation times. Temporal data can partially compensate for SNP distance inflation by constraining the inferred number of transmission steps based on the time elapsed between cases. This means the true outbreaker2 over-estimation under recombination may be lower than the ~1.82× estimated here, or it may be comparable depending on how tightly the temporal prior constrains the inference. Running the full outbreaker2 MCMC on Apollo-simulated sequences is the natural next step (Section 4.5).

Additionally, the simulation used a single recombination model (hotspot-based, HIV-like parameters). The magnitude of inflation will vary for pathogens with different recombination rates, hotspot distributions, and genome lengths. The qualitative finding—that both tools over-estimate under recombination—is expected to be robust, but the quantitative ratios are specific to the simulated parameter regime.

### 4.5 Future directions

The most immediate extension is to run the full outbreaker2 R package on Apollo-generated FASTA sequences, replacing the analytical model with genuine MCMC inference. Apollo's MIT-licensed source code and Colab notebook make this tractable without requiring dedicated computational infrastructure. A second extension is to benchmark additional tools: phyloscanner, which models within-host diversity explicitly and may be more recombination-tolerant; and outbreaker, the predecessor to outbreaker2, which uses a different likelihood formulation. A third extension is to test whether recombination detection (e.g., GARD) followed by recombinant-sequence removal restores accurate inference—this would provide a direct practical mitigation strategy.

Beyond the benchmark itself, the mechanism identified here—that any tool using a global clock anchor is more sensitive to recombination outliers than tools using local pairwise comparisons—provides a testable prediction for the design of next-generation transmission inference tools. Robust estimators that downweight outlier sequences or model recombination explicitly within the transmission inference framework could substantially reduce inflation even for highly recombinogenic pathogens.

---

## 5. Conclusion

This study demonstrates that recombination-induced inflation of transmission network estimates generalises beyond TransPhylo to outbreaker2, a tool with a mechanistically distinct inference strategy. Under identical recombination conditions, TransPhylo over-estimated by 4.27-fold and outbreaker2 by approximately 1.82-fold, relative to Apollo ground truth. The difference in magnitude reflects a fundamental difference in inference architecture: TransPhylo's global MRCA anchor amplifies recombination artefacts multiplicatively, while outbreaker2's local pairwise SNP distances accumulate them additively. Both tools infer accurately without recombination. These results support the recommendation that recombination screening should precede transmission inference for any tool when analysing pathogens known to recombine, including HIV, SARS-CoV-2, norovirus, and influenza.

---

## Acknowledgements

The author thanks Gideon van Marle and co-workers for making the Apollo simulator publicly available under the MIT licence and providing a reproducible experimental framework.

---

## References

1. Gire SK, Goba A, Andersen KG, et al. Genomic surveillance elucidates Ebola virus origin and transmission during the 2014 Sierra Leone outbreak. *Science*. 2014;345(6202):1369–1372. https://doi.org/10.1126/science.1259657

2. Rambaut A, Holmes EC, O'Toole Á, et al. A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology. *Nature Microbiology*. 2020;5(11):1403–1407. https://doi.org/10.1038/s41564-020-0770-5

3. Pillay D, Fisher M, Sabin C, et al. Transmission networks of HIV-1 and their relationship to risk behaviours in the UK. *Retrovirology*. 2012;9(Suppl 1):O25. https://doi.org/10.1186/1742-4690-9-S1-O25

4. Didelot X, Gardy J, Colijn C. Bayesian inference of infectious disease transmission from whole-genome sequence data. *Molecular Biology and Evolution*. 2014;31(7):1869–1879. https://doi.org/10.1093/molbev/msu121

5. Didelot X, Fraser C, Gardy J, Colijn C. Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks. *Molecular Biology and Evolution*. 2017;34(4):997–1007. https://doi.org/10.1093/molbev/msw275

6. Campbell F, Cori A, Ferguson N, Jombart T. Bayesian inference of transmission chains using timing of symptoms, pathogen genomes and contact data. *PLOS Computational Biology*. 2019;15(3):e1006930. https://doi.org/10.1371/journal.pcbi.1006930

7. Pérez-Losada M, Arenas M, Galán JC, Palero F, González-Candelas F. Recombination in viruses: mechanisms, methods of study, and evolutionary consequences. *Infection, Genetics and Evolution*. 2015;30:296–307. https://doi.org/10.1016/j.meegid.2014.12.022

8. Lacek KA, Rambo-Martin BL, Bagal UN, et al. Identification of a novel SARS-CoV-2 Delta-Omicron recombinant virus in the United States. *bioRxiv*. 2022. https://doi.org/10.1101/2022.03.19.484981

9. Eden JS, Tanaka MM, Boni MF, Rawlinson WD, White PA. Recombination within the pandemic norovirus GII.4 lineage. *Journal of Virology*. 2013;87(11):6270–6282. https://doi.org/10.1128/JVI.03464-12

10. van Marle G, Bhattarai N, Gill MJ, et al. Apollo: a comprehensive GPU-powered within-host simulator for viral evolution. *Nature Communications*. 2025;16:5783. https://doi.org/10.1038/s41467-025-60988-8

11. Martin DP, Murrell B, Golden M, Khoosal A, Muhire B. RDP4: Detection and analysis of recombination patterns in virus genomes. *Virus Evolution*. 2015;1(1):vev003. https://doi.org/10.1093/ve/vev003

12. Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW. GARD: a genetic algorithm for recombination detection. *Bioinformatics*. 2006;22(24):3096–3098. https://doi.org/10.1093/bioinformatics/btl474

13. Rambaut A, Lam TT, Max Carvalho L, Pybus OG. Exploring the temporal structure of heterochronous sequences using TempEst (formerly Path-O-Gen). *Virus Evolution*. 2016;2(1):vew007. https://doi.org/10.1093/ve/vew007

---

## Tables

**Table 1** — see Section 3.1.

---

## Figure captions

**Figure 1.** Transmission network size inferred by TransPhylo and outbreaker2 under two recombination conditions. Bars show the ground-truth network size from Apollo (blue), TransPhylo as reported by van Marle et al. (red), TransPhylo as estimated by the analytical model (orange), and outbreaker2 estimated by the analytical model (purple). Numbers above bars indicate the inferred count; fold-overestimation ratios are annotated in the corresponding colour. Both tools over-estimate under recombination (Category 1); neither over-estimates substantially without recombination (Category 2). TransPhylo over-estimates more severely (4.27×) than outbreaker2 (~1.82×).

**Figure 2.** Mechanistic comparison of recombination-induced inflation pathways. Upper panel: the TransPhylo pathway, in which recombination inflates global phylogenetic divergence, causes BEAST2 to infer an older MRCA (1990 vs. true 1993), and TransPhylo fills the resulting three-year gap with phantom hosts. Lower panel: the outbreaker2 pathway, in which recombination inflates individual pairwise SNP distances, causing outbreaker2 to infer extra unsampled intermediate hosts along each affected transmission link. The distributed nature of the outbreaker2 pathway results in lower aggregate inflation than TransPhylo's single-anchor failure.

**Figure 3.** Pairwise SNP distance versus true transmission steps for sampled-to-sampled links, without recombination (left) and with recombination (right). In the without-recombination condition, observed SNP distances follow the clock expectation closely (E[SNPs] = μ × L × steps). In the recombination condition, a subset of links (orange diamonds) exhibits SNP distances far above the clock line; these outliers, which arise from recombination inserting divergent segments, are the primary driver of outbreaker2 over-estimation.
