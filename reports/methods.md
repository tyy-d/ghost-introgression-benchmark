# Method adaptations and interpretation

## ArchIE

Use upstream `N_ton`, `distance_vector`, and `min_dist_ref` from commit
`00820bc61beca8a958f04fb17291ff66e03b9e94`. Reimplement the same S* score/recurrence
with vectorized operations and consistent physical positions after both reference
and focal-site filtering. The published scripts remove genotype columns without
consistently removing their positions; the adapter fixes that alignment. Test the
recurrence against upstream on aligned inputs and against a known filtered example.

Use a dynamically sized feature vector (2n + 8 features for n target haplotypes),
instead of the original R script's hard-coded 209th label column. Fit the same
unpenalized binomial GLM with R `glm`, standardizing nonconstant columns for numerical
conditioning. Drop constant features and replace undefined skew/kurtosis with zero.
The R adapter uses base R, removing the incidental `data.table` dependency from the
original training script. Do not substitute human-trained model weights.

Training labels come from exact simulated migration intervals, not an interpolation
between SNP labels. Use >70% donor ancestry as positive and <30% as negative;
exclude mixed windows from fitting. Balance windows across introgression versus
no-introgression scenarios, not positive versus negative haplotype labels.

These are documented adaptations of ArchIE, not an unmodified reproduction of
every implementation detail. The full empirical pipeline additionally needs
ancestral-state orientation for allele-frequency and private-derived-allele features.

## TRACE and SINGER

Use upstream TRACE commit `70cab6c6687900826eb0b2ded69c39554b1b4146` and SINGER
commit `eb8e39b1a15be4a9a4df4fdaab61847bf73515d7`. Extract genealogical statistics
on a 1-kb grid for each retained posterior ARG and average them before fitting the
upstream TRACE HMM. State 1 is archaic, following upstream output utilities.

Run SINGER with the configured mutation/recombination ratio, population size, and
polarized VCF. Convert retained posterior samples with its accompanying converter.
Use the upstream `singer_master` wrapper, which can recover from native crashes
using its documented debug/resume path; retain its logged seeds and retries.
Account separately for ARG inference and TRACE computation and report their sum.
The short smoke chain is only a functional test; inspect convergence before a
scientific performance comparison. True ARGs are separately labeled diagnostics.

## Evaluation

Project coordinates are zero-based half-open, VCF is one-based. Both predictions
and truth are integrated over the same callable bases and 50-kb grid. Classify
grid cells by their callable-weighted truth fraction; report ambiguous cells
separately. Precision–recall uses probability thresholds, and the full benchmark
freezes an F1-optimal validation threshold before test-scenario evaluation.

The structure-only scenario allows low continuous migration between ancestral
populations after target/reference divergence, with no recent pulse. Any detections
are reported as departures under this alternative history, not proof that the
detector is wrong: that model does include ancient movement of lineages.
