# Research report draft

Status: provisional; see run eligibility and convergence evidence. Not an approval record.

## Detecting ancestry from unsampled lineages: introduction

Recombination makes a genome a mosaic of histories rather than a single family tree. When formerly separated populations interbreed, descendants can inherit segments whose ancestry differs from that of neighboring sequence. Introgression from an unsampled donor is particularly difficult to detect because no donor genome is available for direct matching. The word ghost describes that absence of a sampled reference. It does not establish that the donor was extinct, a separate species, or identifiable from a fossil.

Two computational approaches can use different aspects of the same descendant sequences. ArchIE combines population-genetic summaries in a simulation-trained binomial classifier. Its features include allele frequencies, distances among haplotypes, divergence from a comparison population, and patterns of private derived variation. Training links those features to known simulated ancestry. Consequently, the fit between training and test histories is part of the scientific question, not merely an installation detail.

TRACE instead examines genealogical information in ancestral recombination graphs. Deep branches that persist across linked sequence can support an introgression interpretation. The method combines information along the chromosome with a hidden Markov model. In this project, SINGER reconstructs graphs from the observed genomes before TRACE runs. These graphs are estimates with uncertainty, so a detector can lose information before ancestry classification begins. A true-graph diagnostic separates that reconstruction step from the detector itself.

Durvasula and Sankararaman (2020) used ArchIE alongside a separate allele-frequency analysis to investigate archaic ancestry in African populations. Zhang and colleagues (2026) introduced TRACE using inferred genealogies. The present comparison concerns ancestry-segment detection with adapted software. It does not reproduce the 2020 paper's separate demographic-model analysis, and it does not presume either method will be superior.

## Introduction: question and data

The precise analysis question is how sensitivity, achieved false discovery rate, and computational cost change when donor ancestry becomes older or rarer and when assumptions about a comparison population or demographic parameters fail. A secondary question asks how much TRACE changes when it receives inferred rather than true graphs. The primary comparison freezes thresholds using independent validation data, then evaluates paired methods on exactly the same held-out sequences.

The main dataset is generated with msprime. Each independent region contains ten diploid target and ten diploid reference individuals across one million bases. The donor has no sampled genomes. Migration records identify ancestry along the simulated chromosomes, including intervals ancestral to multiple descendants. Phased biallelic sequence data and sample order are exported separately from the true graphs and ancestry labels. Training, validation and test seeds occupy disjoint ranges.

The controlled baseline uses effective size 100,000, mutation and recombination rates of 1e-8 per base per generation, donor divergence 200,000 generations ago, target/reference divergence 50,000 generations ago, and five percent admixture 2,000 generations ago. These are controlled software-study parameters, not estimates of a mosquito history. The data contain neutral simulated sequence and known migration histories; they do not contain selection, functional annotations, real sequencing errors, or direct evidence of a particular extinct organism.

Public Ag1000G genomes motivate a conditional nonhuman demonstration, but empirical ancestry has no equivalent exact truth labels. Documented sample quality, relatedness, ancestral allele orientation and biological parameters are prerequisites. In particular, a reference allele is not automatically ancestral. The simulation study can quantify accuracy; empirical method agreement can only establish consistency. The assignment permits simulation projects, while approval of this specific proposal remains a separate course requirement.

## Results: evidence available and validity

Only diagnostic results are available. All numerical values below use an unfrozen 0.5 threshold and do not estimate held-out calibrated performance.

The convergence pilot did not complete: it was manually interrupted after about 23 minutes, and a PermissionError during cleanup caused the run to be recorded as failed. Therefore, convergence has not been assessed. The planned checks use two independent chains and examine split R-hat and effective sample size for diversity fit, incompatible mutation sites, and TRACE ancestry amount on a predeclared focal haplotype. Constant statistics alone cannot establish mixing. Passing this pilot would support a choice of settings; it would not prove convergence in every test history.

The primary table pools exact callable-base counts. It must be read together with the per-condition tables and paired bootstrap results: pooling different conditions can conceal failure at low ancestry prevalence or old admixture. Undefined FDR or precision means there are no called bases, not perfect control of false discoveries. True-ARG TRACE is a diagnostic and is never a deployable empirical alternative.

## Results: comparisons, tests and resource cost

The primary endpoint is recall at thresholds targeting estimated validation FDR no greater than five percent. The planned calibration will use six independent one-megabase validation regions for each method. If no nonempty call set qualifies, that method will abstain at this operating target. The test result will be achieved FDR alongside recall; a validation target cannot guarantee the same error rate under changed demographic histories.

Uncertainty resamples independent simulation replicates and preserves method pairing. Five or more usable paired replicates are required for confidence intervals; smaller sets remain descriptive. Failed, incomplete and corrupted outputs are counted in the run-status table rather than imputed as zero ancestry. An alternative-history scenario contains ancestral migration, so its calls are not interpreted as ordinary false positives.

The first principal figure summarizes recall and achieved FDR. The second shows an ancestry-localization example chosen by a deterministic largest-error rule for eligible test data; a diagnostic-only draft instead uses the first target haplotype explicitly. The third compares inference cost and available robustness conditions. ArchIE training is recorded once and excluded from its per-region inference time; SINGER reconstruction is included in TRACE's total. CPU and memory measurements include availability notes because the local sandbox may restrict process-tree sampling.

Correctness checks cover exact mixed-base overlap, masks and adjacent intervals, null data, unattainable operating points, stale hashes, validation/test separation, paired resampling, and process-group cancellation. These checks verify scientific bookkeeping; they do not establish that the demographic model captures a natural population.

## Discussion: strengths, limitations and unanswered questions

ArchIE turns interpretable population summaries into inexpensive repeated predictions, but its training distribution, comparison population and fixed window resolution can limit transfer. TRACE uses richer genealogical information and finer output intervals, but reconstructing those histories costs time and introduces uncertainty before detection. The true-graph diagnostic is useful precisely because it separates these two limitations.

The controlled simulation supplies exact ancestry and excludes several natural processes. A high score can show success under specified assumptions, not species-wide validity. Reference admixture and demographic mismatch stress tests probe only two departures. A short chain, few independent regions, or no qualifying FDR threshold can make the answer inconclusive. Those are reportable limits rather than reasons to alter thresholds after seeing test labels or discard difficult regions.

The present deliverable is a research draft for the November milestone. Written course approval has not been recorded. Final readiness is targeted for December 7, ahead of the assignment's December 8 deadline. Empirical mosquito inference remains conditional on QC, unrelated samples, polarization and documented parameters. Agreement would not identify an extinct donor or supply an accuracy estimate.

References: Durvasula and Sankararaman (2020), Science Advances, DOI 10.1126/sciadv.aax5097. Zhang et al. (2026), Science, DOI 10.1126/science.aef8874. Software and full source/version records: github.com/sriramlab/ArchIE; github.com/YulinZhang9806/trace; github.com/popgenmethods/SINGER; tskit.dev/msprime. Public empirical resource: malariagen.net.

## Figures: evidence and reproducible exports
