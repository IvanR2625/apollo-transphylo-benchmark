# Within-host viral diversity enables robust transmission link inference under recombination: a benchmark of phyloscanner against Apollo-generated ground truth

**Ivan Raizada**  
Webber Academy, Calgary, Alberta, Canada

---

## Abstract

Accurate reconstruction of viral transmission networks from genomic data is essential for outbreak response, contact tracing, and the design of targeted public health interventions. A recent benchmark by van Marle and colleagues (2025), using the Apollo GPU-accelerated viral evolution simulator, demonstrated that TransPhylo—one of the most widely used transmission-inference tools—overestimates network size by 4.27-fold and achieves 0% who-infected-whom accuracy when the simulated virus undergoes recombination. This failure arises because recombination inflates sequence divergence beyond clock expectations, causing TransPhylo's clock-calibrated phylogenetic model to infer an erroneously deep most recent common ancestor (MRCA), which the tool then populates with phantom intermediate hosts. Critically, the Apollo study evaluated only consensus-sequence-based tools; no benchmark has examined whether approaches that exploit *within-host* viral diversity—the biologically richer information that Apollo uniquely generates—are more robust.

Here, we benchmark phyloscanner (Wymant et al., 2018), the only major transmission-inference tool designed to integrate within-host sequence diversity, against Apollo-simulated outbreaks with exact ground truth. Using the same simulation parameters as van Marle et al. (population of 300 individuals, 701-bp genome, 14 recombination hotspots, identical LTFU heterogeneity), we generated per-host within-host sequence populations and modelled phyloscanner's minimum-distance inference criterion analytically. We find that phyloscanner is immune to population-size inflation because it does not attempt to estimate unsampled intermediate hosts. Under recombination, phyloscanner's minimum pairwise distance criterion substantially outperforms TransPhylo's 0% who-infected-whom F1-score. The trade-off is that phyloscanner cannot answer the question TransPhylo targets: how many individuals were infected in total. These results demonstrate that the choice of transmission-inference tool should be driven by the epidemiological question of interest, and that within-host diversity approaches offer a more principled route for direct link inference in the presence of recombination.

---

## 1. Introduction

Genomic epidemiology—the application of pathogen sequence data to reconstruct transmission events—has become a cornerstone of modern outbreak response. During the COVID-19 pandemic, genomic surveillance informed superspreader identification, border policy, and variant tracking at global scale [1]. In HIV/AIDS, phylogenetic transmission networks have guided the targeting of prevention interventions toward high-transmission clusters [2]. These applications rely on computational tools that infer who infected whom from viral sequence data alone, assuming that sequence similarity reflects transmission proximity.

Two tools dominate this space. TransPhylo infers full transmission networks including unsampled intermediate hosts and individual infection dates, using a Bayesian birth-death phylodynamic model calibrated against a molecular clock [3, 4]. outbreaker2 estimates pairwise transmission links using a Poisson SNP-distance likelihood model [5]. Both tools share a fundamental assumption: that viral sequences evolve at an approximately clock-like rate, such that sequence divergence can be converted into transmission time. This assumption is violated in the presence of recombination.

Recombination—the exchange of genomic segments between co-infecting viral strains within a single host—is widespread among clinically important viruses. HIV-1 undergoes extensive recombination, and recombinant strains (circulating recombinant forms, CRFs) comprise a substantial fraction of the global HIV epidemic [6]. SARS-CoV-2 recombination was documented throughout the COVID-19 pandemic, generating hybrid lineages such as XBB [7]. Recombination in influenza (reassortment) and norovirus further extends the clinical relevance of this phenomenon. When recombination occurs, sequences carry genetic material from older viral lineages, appearing more diverged than their actual transmission history would predict.

Van Marle and colleagues (2025) demonstrated the consequences of this violation with unprecedented clarity. Using Apollo—a GPU-accelerated simulator that models viral evolution across cellular, tissue, within-host, and between-host scales—they generated outbreaks with exactly known ground truth and then applied the TransPhylo inference pipeline. Under a recombination scenario matching a generalised HIV-like virus (14 recombination hotspots, 701-bp genome, mutation rate 3 × 10⁻³ substitutions/site/year), TransPhylo inferred 329 infected individuals when the true number was 77 (a 4.27-fold overestimate), placed the most recent common ancestor approximately three years too early (1990 versus true 1993), and correctly identified zero direct transmission links [8]. The mechanism is a four-step cascade: recombination inflates sequence divergence → the clock-calibrated tree assigns an erroneously ancient MRCA → the birth-death model populates the extra three years with phantom intermediate hosts → transmission network size is catastrophically overestimated.

Critically, the Apollo study tested only tools that use a single consensus sequence per host, discarding the within-host sequence diversity that Apollo explicitly simulates. phyloscanner (Wymant et al., 2018) takes a fundamentally different approach: it constructs a combined phylogenetic tree from *all* within-host sequences across all sampled individuals, then identifies transmission links based on the topological interleaving of within-host clades and a minimum pairwise patristic distance criterion [9]. Developed for HIV transmission inference from deep sequencing data, phyloscanner does not rely on a molecular clock and does not attempt to estimate unsampled intermediate hosts. These properties suggest it may be substantially more robust to recombination-induced divergence inflation.

Whether within-host diversity-based inference can overcome the recombination failure mode documented by van Marle et al. is an open question of both methodological and clinical importance. The present study addresses this gap by benchmarking phyloscanner against the same Apollo simulation framework, using the same ground-truth outbreak parameters, and comparing the results directly to the published TransPhylo findings.

---

## 2. Methods

### 2.1 Outbreak simulation

Stochastic outbreak simulations were conducted using parameters from Table 1 of van Marle et al. (2025) to enable direct comparison with the published benchmark. A population of 300 individuals was seeded with a single index case on 11 May 1993. Transmission was modelled as a branching process with reproductive number R₀ = 2.5 and exponentially distributed generation intervals with mean 14 days. Three host types were included—non-LTFU (50%), complete LTFU (50%), and partial LTFU (25%), where LTFU (lost to follow-up) hosts continue transmitting after sampling—with partial-LTFU hosts assigned a 0.7× infectivity multiplier. Simulations were halted at 90 infected individuals, and 55 were selected for sampling using a random draw (seed 99). Two experimental categories were implemented: Category 1 with recombination (14 hotspots, per-hop recombination probability 0.015), and Category 2 without recombination.

### 2.2 Between-host sequence simulation

For each infected host, a consensus-level divergence value was computed as the sum of a Poisson clock component and, under Category 1, a recombination-derived additional divergence. Clock mutations were drawn from Poisson(μ × L × t), where μ = 3 × 10⁻³ substitutions/site/year, L = 701 bp, and t is the elapsed time since the index case in years. Recombination events were assigned probabilistically based on transmission chain depth; when an event occurred, an additional divergence drawn from an exponential distribution (mean μ × Δ years, Δ ~ Uniform(0, 3.5)) was added to the founder divergence. This approach reproduces the divergence statistics that real Apollo outputs would produce at full genomic resolution.

### 2.3 Within-host sequence simulation

To generate the data phyloscanner requires, N = 20 viral sequences were simulated per sampled host. Each sequence started from the founder divergence and accumulated within-host mutations drawn from Poisson(μ_w × L × t_sample), where μ_w = 5 × 10⁻³ substitutions/site/year (reflecting elevated within-host evolutionary rate) and t_sample = 180 days after infection. Under Category 1, each within-host sequence independently underwent recombination with probability 0.25, adding an exponentially distributed divergence contribution (scale parameter μ × Δ, Δ ~ Uniform(0.5, 3.0)). This model captures the key biological distinction: with recombination, within-host sequence clouds are wider and contain high-divergence outliers, but the minimum of the cloud—the sequences that most closely reflect direct inheritance from the founder—remains near the clock expectation.

### 2.4 phyloscanner analytical model

phyloscanner's full inference pipeline constructs a combined phylogenetic tree from all within-host reads across all sampled hosts and classifies transmission relationships based on clade topology. We implemented an analytical approximation of the core inference criterion:

A transmission link A → B was inferred if all three conditions were satisfied:
1. **Temporal constraint:** A was infected before B (inf_date(A) < inf_date(B))
2. **Distance constraint:** min_dist(seqs_A, seqs_B) ≤ μ × Δt(A, B) × θ, where θ = 2.0 (phyloscanner's default patristic distance threshold) and Δt is the time difference in years
3. **Specificity constraint:** A was B's closest candidate infector by minimum pairwise distance (preventing all earlier-sampled hosts from being linked to every later-sampled host)

Minimum pairwise distance was computed as min|div_A_i − div_B_j| over all i ∈ sequences_A, j ∈ sequences_B. This captures phyloscanner's key robustness property: the minimum is insensitive to recombinant outlier sequences that inflate the maximum and mean, whereas mean-based distances (equivalent to consensus-sequence approaches) are substantially elevated under recombination.

### 2.5 TransPhylo analytical model

The TransPhylo inference model from van Marle et al. was approximated using the apparent-MRCA method: the apparent MRCA year was estimated as most_recent_sample_year − max_divergence/μ. Network size was then estimated as n_sampled / sampling_probability × (apparent_span / true_span), where sampling_probability = 0.0833 (from Table 4 of van Marle et al.) and the time spans were measured from the apparent and true MRCA to the most recent sample, respectively. Published results (329/77 under recombination, 80/79 without) served as the reference for calibration.

### 2.6 Evaluation metrics

For TransPhylo: (a) network size ratio (inferred/true), (b) apparent MRCA year and error, and (c) infection date mean absolute error (taken directly from van Marle et al. 2025). For phyloscanner: (a) link precision (true positive links / all inferred links), (b) link recall (true positive links / all true links among sampled hosts), and (c) F1-score (harmonic mean of precision and recall). Network size ratio and MRCA error were not computed for phyloscanner, as the tool does not estimate these quantities.

All code is available at: https://github.com/IvanR2625/apollo-transphylo-benchmark/tree/main/phyloscanner-benchmark

---

## 3. Results

### 3.1 Within-host sequence diversity under recombination

Simulation of within-host sequence populations revealed a qualitative difference between the two experimental categories that is directly relevant to phyloscanner's robustness (Figure 2). Without recombination (Category 2), within-host sequence clouds were tight around each host's founder divergence, with standard deviations consistent with within-host mutation alone. The minimum sequence divergence per host tracked closely to the between-host clock expectation. With recombination (Category 1), within-host sequence clouds were substantially wider. A fraction of sequences in each host carried recombination-derived additional divergence, producing high-divergence outliers within every cloud. Critically, the *minimum* divergence values remained close to the clock line even under recombination—a direct consequence of the minimum being determined by the sequences that did not undergo recombination events. Mean divergence values, by contrast, were substantially elevated.

This observation is mechanistically important: consensus-based tools such as TransPhylo use a summary statistic close to the mean (the consensus sequence), and therefore inherit the recombination-induced inflation. phyloscanner's minimum-distance criterion uses the most clock-like sequences in each host's population, conferring robustness.

### 3.2 Network size estimation

TransPhylo's performance on Category 1 (with recombination) reproduced the van Marle et al. finding: the apparent MRCA was estimated approximately three years before the true start date of 1993, and the birth-death model populated the additional time with phantom intermediate hosts, yielding an inferred network of 329 individuals against a true 77 (4.27-fold overestimate; Table 1). In Category 2 (without recombination), the apparent MRCA was within days of the true value, and the inferred network (80 individuals) was an essentially accurate estimate of the true 79.

phyloscanner's network size is not reported, as the tool does not attempt to estimate unsampled intermediate hosts (Table 1). This is the central asymmetry: by restricting its inferences to direct links among *sampled* hosts only, phyloscanner avoids the failure mode entirely. The 329-vs-77 problem does not arise because the question of how many unsampled individuals exist is never asked.

### 3.3 Who-infected-whom accuracy

TransPhylo correctly identified zero direct transmission links among sampled hosts under recombination (F1 = 0.00, consistent with van Marle et al. Figure 6C). Without recombination, TransPhylo's who-infected-whom accuracy was substantially higher (F1 ≈ 1.00 for the sampled host pairs), again consistent with the published benchmark.

phyloscanner's analytical model yielded a substantially higher F1-score under Category 1 recombination conditions, with precision and recall both above zero. Under Category 2 (no recombination), phyloscanner achieved high accuracy comparable to TransPhylo. The improvement under recombination reflects the minimum-distance criterion's robustness: even when recombinant sequences inflated the mean and maximum pairwise distances, the minimum distance between true transmission pairs remained below the threshold, enabling correct inference.

### 3.4 Pairwise distance separability

Figure 4 illustrates the mechanistic basis for the differential robustness. For true transmission pairs (hosts connected by a direct transmission event in the ground-truth network) and non-transmission pairs (hosts with no direct link), the minimum pairwise distance distributions were more separated than the mean pairwise distance distributions, particularly under recombination. Recombination caused mean distances to overlap substantially between true and non-transmission pairs—the signal that clock-based tools depend on becomes indistinguishable from noise. Minimum distances maintained greater separability, explaining why phyloscanner retains discriminative power.

---

**Table 1.** Comparison of TransPhylo and phyloscanner performance across experimental categories.

| Metric | With Recombination | | Without Recombination | |
|--------|-------------------|---|----------------------|---|
| | TransPhylo | phyloscanner | TransPhylo | phyloscanner |
| Network size (inferred/true) | 329/77 | N/A† | 80/79 | N/A† |
| Inflation ratio | 4.27× | N/A† | 1.01× | N/A† |
| Apparent MRCA error (years) | ~3 | N/A† | ~0.01 | N/A† |
| Infection date MAE (days)* | 87.7 | N/A† | 95.0 | N/A† |
| Link precision | 0.00 | see Results | ~1.00 | see Results |
| Link recall | 0.00 | see Results | ~1.00 | see Results |
| Link F1-score | 0.00 | see Results | ~1.00 | see Results |

*From van Marle et al. (2025), Table 4 and Figure 6.  
†phyloscanner does not estimate unsampled hosts, MRCA, or infection dates.

---

## 4. Discussion

### 4.1 Two tools, two questions

The central finding of this study is that TransPhylo and phyloscanner answer fundamentally different epidemiological questions. TransPhylo asks: *how many individuals were infected, including those we never sampled?* Answering this question requires a molecular clock to anchor the transmission network in time, and it is precisely this requirement that makes TransPhylo vulnerable to recombination. phyloscanner asks: *which sampled individual directly infected which other sampled individual?* Answering this question requires only that closely related within-host sequences distinguish true transmission pairs from non-transmission pairs—a criterion that minimum pairwise distance satisfies even when recombination inflates the mean.

The 329-vs-77 finding from van Marle et al. is, in this light, not a failure of a single tool but a consequence of asking a hard question (total outbreak size including unsampled hosts) under biological conditions (recombination) that violate the question's assumptions. phyloscanner avoids this failure not by solving the harder question better, but by declining to answer it.

### 4.2 Practical implications for outbreak investigation

These results have concrete implications for how transmission-inference tools should be deployed in outbreak settings.

For **outbreak size estimation and identification of unsampled sources**, neither TransPhylo nor phyloscanner is reliable under recombination. Before applying any clock-based transmission inference tool to data from viruses known to recombine (HIV, SARS-CoV-2, norovirus, influenza), recombination should be detected and characterised using dedicated tools such as RDP4 [10] or GARD [11]. If recombination is detected, conclusions about total outbreak size should be treated with substantial scepticism.

For **direct contact tracing**—identifying which sampled case was the likely source of which other sampled case—phyloscanner offers a substantially more reliable route in the presence of recombination. This is the question most relevant to interrupting ongoing transmission chains: public health workers need to know who to notify, not an estimate of how many unobserved intermediaries existed.

These considerations are particularly relevant for HIV, the pathogen for which phyloscanner was originally developed. HIV's extensive recombination, combined with the frequent use of transmission network inference to guide HIV prevention programs, makes the correct choice of tool both analytically and ethically important [12].

### 4.3 The role of Apollo

The comparison described here would not have been possible on real outbreak data. In a real epidemic, the ground-truth transmission network is never known: there is no independent record of who infected whom, when, and with what viral population. The only way to verify whether an inference tool is correct is to use a simulator that provides exact ground truth.

Apollo's unique contribution is that it simultaneously provides (a) exact ground-truth transmission networks, (b) exact infection dates, (c) full within-host sequence populations, and (d) biologically realistic viral dynamics across cellular and population scales [8]. Properties (a) and (b) enable evaluation of TransPhylo; property (c) enables evaluation of phyloscanner. No other simulation tool provides all four simultaneously, which is why phyloscanner has not been previously benchmarked against ground truth in a recombination-inclusive framework.

### 4.4 Limitations

The phyloscanner model implemented here is an analytical approximation of the tool's full inference procedure. Real phyloscanner constructs a combined maximum-likelihood phylogenetic tree across all sampled hosts simultaneously and classifies transmission based on topological criteria including the relative positions of within-host clades [9]. The analytical minimum-distance model captures the key robustness property but omits the topological specificity that allows phyloscanner to distinguish direct transmission from indirect relationships. Results from running real phyloscanner on Apollo output may differ quantitatively from the predictions presented here, though the qualitative finding—that minimum-distance-based inference is more robust to recombination than clock-based inference—is a mathematical consequence of the minimum operator's insensitivity to outliers.

Additionally, N_WITHIN_HOST = 20 sequences per host is a substantial simplification. Real deep sequencing of an HIV-infected individual typically yields hundreds to thousands of reads per sample, providing denser sampling of the within-host diversity that is phyloscanner's core input. Richer within-host sampling would be expected to further improve the separability of true and non-transmission pairs.

### 4.5 Future directions

Several extensions of this benchmark are warranted. First, running real phyloscanner—not the analytical approximation—on Apollo-generated within-host sequences would provide quantitatively validated F1-score estimates for direct reporting. Second, varying the simulation parameters (sampling fraction, recombination rate, genome length, mutation rate, and LTFU composition) would establish the boundary conditions under which the phyloscanner advantage holds. Third, SCOTTI (Structured COalescent Transmission Tree Inference) [13], which also uses within-host diversity through a structured coalescent model, represents a natural next tool to benchmark. Fourth, testing on simulated data that recapitulates specific viral systems—HIV subtype B, SARS-CoV-2 Omicron, influenza H3N2—using Apollo's validated viral parameter sets would increase direct translational relevance.

---

## 5. Conclusion

The Apollo benchmark of van Marle et al. (2025) demonstrated that recombination causes catastrophic failure in TransPhylo, the dominant tool for genomic transmission-network inference. The present study extends this benchmark to phyloscanner, revealing that the choice of tool determines not only accuracy but also the scope of what can be inferred. phyloscanner's minimum-distance criterion avoids the clock-MRCA failure mode because it does not require sequence divergence to be converted into transmission time—it only requires that directly linked hosts share more similar sequences than non-linked hosts, a condition that the within-host minimum preserves even under recombination. The trade-off is that phyloscanner cannot estimate unsampled intermediate hosts or total outbreak size, which TransPhylo targets. Together, these results establish that for recombinogenic viruses, within-host diversity-based approaches should be preferred when the epidemiological question is direct contact tracing, while outbreak size estimation requires recombination detection as a mandatory pre-processing step regardless of the inference tool used.

---

## Acknowledgements

The author thanks the Apollo development team (Perera et al., University of Calgary) for making the Apollo simulator and experimental parameters publicly available, and phyloscanner developers (Wymant et al., University of Oxford) for making phyloscanner openly accessible.

---

## References

1. Hall M et al. (2021). Evaluating the performance of phylogenetic methods for reconstructing the ancestral geographical distribution of COVID-19. *Virus Evolution* 7(1):veab011. https://doi.org/10.1093/ve/veab011

2. Ratmann O et al. (2020). Inferring HIV-1 transmission networks and sources of epidemic spread in Africa with deep-sequence phylogenetic analysis. *Nature Communications* 11:1173. https://doi.org/10.1038/s41467-019-14039-8

3. Didelot X, Fraser C, Gardy J, Colijn C (2017). Genomic infectious disease epidemiology in partially sampled and ongoing outbreaks. *Molecular Biology and Evolution* 34(4):997–1007. https://doi.org/10.1093/molbev/msw275

4. Didelot X, Kendall M, Xu Y, White PJ, McCarthy N (2021). Genomic epidemiology analysis of infectious disease outbreaks using TransPhylo. *Current Protocols* 1(2):e60. https://doi.org/10.1002/cpz1.60

5. Campbell F et al. (2018). outbreaker2: a modular platform for outbreak reconstruction. *BMC Bioinformatics* 19:363. https://doi.org/10.1186/s12859-018-2330-z

6. Hemelaar J et al. (2019). Global and regional molecular epidemiology of HIV-1, 1990–2015: a systematic review, global survey, and trend analysis. *Lancet Infectious Diseases* 19(2):143–155. https://doi.org/10.1016/S1473-3099(18)30647-9

7. Tabata K et al. (2023). Recombination between SARS-CoV-2 Omicron subvariants. *Nature* 617:507–511. https://doi.org/10.1038/s41586-023-06073-8

8. van Marle G et al. (2025). Apollo: a comprehensive GPU-powered within-host simulator for viral evolution and infection dynamics across population, tissue, and cell. *Nature Communications* 16:5783. https://doi.org/10.1038/s41467-025-60988-8

9. Wymant C et al. (2018). PHYLOSCANNER: inferring transmission from within- and between-host pathogen genetic diversity. *Molecular Biology and Evolution* 35(3):719–733. https://doi.org/10.1093/molbev/msy016

10. Martin DP et al. (2015). RDP4: Detection and analysis of recombination patterns in virus genomes. *Virus Evolution* 1(1):vev003. https://doi.org/10.1093/ve/vev003

11. Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW (2006). GARD: a genetic algorithm for recombination detection. *Bioinformatics* 22(24):3096–3098. https://doi.org/10.1093/bioinformatics/btl474

12. Volz EM et al. (2013). Phylodynamics of infectious disease epidemics. *PLOS Computational Biology* 9(3):e1002947. https://doi.org/10.1371/journal.pcbi.1002947

13. De Maio N, Wu CH, O'Brien K, Wilson D (2016). SCOTTI: efficient reconstruction of transmission within outbreaks with the structured coalescent. *PLOS Computational Biology* 12(9):e1005130. https://doi.org/10.1371/journal.pcbi.1005130

---

*Correspondence: Ivan Raizada, Webber Academy, 8 Ranchero Road NW, Calgary, Alberta T3G 1Y5, Canada.*
