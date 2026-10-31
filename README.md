# Ghost Introgression Benchmark

Compare **ArchIE** and **SINGER + TRACE** for introgression from an unsampled donor,
using simulations with known ancestry and a bounded *Anopheles* empirical pilot.
This is a computational biology class project with a 2–3-week research scope.

**Status:** code and unit checks exist, but no successful full method run is claimed.
See [setup status](reports/setup-status.md) for verified work and outstanding blockers.

## Research question

How do detection accuracy, false positives, and computational cost change with
introgression age and amount when the donor is unsampled? A secondary question is
how much TRACE's performance depends on reconstructing the ARG.

The comparison is **TRACE versus ArchIE**, not a reproduction of the 2020 paper's
separate conditional site-frequency-spectrum demographic inference. Simulations
provide accuracy measurements. Empirical concordance provides biological context,
not truth labels or evidence sufficient to name an extinct donor.

## Set up

Requirements: Conda, Git, a C++17 compiler, and network access. The bootstrap script
installs only within this project, apart from normal system-provided tools. It
downloads upstream repositories into ignored `vendor/` directories, checks out
recorded commits, installs R and Python dependencies, and builds SINGER.

```bash
bash scripts/bootstrap.sh
conda activate "$PWD/.tools/environment"
pytest -q
ruff check src tests
```

The environment uses Python 3.11 and R; no human demographic parameters are used
as mosquito estimates. `RSCRIPT` and `SINGER_BINARY` can override executable paths.
The downloaded upstream methods retain their original authorship and terms.
Compatibility adaptations are listed in [methods notes](reports/methods.md).

## Run stages

Run from the repository root after activating the environment:

```bash
ghostbench simulate --config configs/smoke.yaml --output results/smoke
ghostbench train --config configs/smoke.yaml --output work/archie
ghostbench run --config configs/smoke.yaml --output results/smoke --model work/archie/model.rds
ghostbench evaluate --output results/smoke
ghostbench report --output results/smoke
```

Or run `ghostbench smoke --config configs/smoke.yaml --output results/smoke` to train
the model if necessary and execute all synthetic stages. That command runs actual
ArchIE feature extraction and an R binomial GLM, actual SINGER inference, and the
upstream TRACE HMM. It also runs TRACE on true ARGs as a separately named diagnostic.
There is no mock-method fallback.

The software-validation smoke configuration uses 10 diploid target and 10 diploid
reference samples, 1 Mb, seed 42, and 2,000 independent training windows. The
inference command does not access the truth directory; only evaluation and the
explicit true-ARG diagnostic do. Training, validation, and test seeds are separate.
The single smoke run uses a fixed 0.5 probability cutoff; it does not tune on test labels.

For the replicated study, run inside a university allocation:

```bash
bash scripts/run_hpc.sh
```

This calibrates thresholds on a held-out 10-Mb validation sequence, then freezes
them for all six scenarios and ten replicates per scenario. The runner resumes
completed replicates. Local training and full evaluation can be computationally
substantial; the short smoke MCMC chain is not suitable for claims of convergence.
The runner uses normal sequential processes and a small process pool for training;
it does not assume your university's scheduler, account, or partition.

## Inputs and outputs

- `inputs/`: phased biallelic VCF, explicit sample/haplotype order, callable intervals.
- `truth/`: simulation ARG and true donor intervals; never a method input.
- `archie/`, `trace/`, `trace_true_arg/`: predictions, scores, and fit logs.
- Experiment root: configuration, tool versions and commits, runtime, metrics,
  precision–recall figure, and generated Markdown report.

Predictions are tab-separated with columns `replicate`, `method`, `haplotype`,
`chromosome`, `start`, `end`, `probability`. Coordinates are **zero-based half-open**;
VCF positions are one-based. Missing calls mean probability zero; overlapping calls
are rejected. Metrics use a common 50-kb grid weighted by callable bases. True
ancestry fractions below 0.3 or above 0.7 define negative/positive windows; mixed
windows are counted separately. Undefined precision or recall is reported as null.

## Empirical pilot

The data target is Ag1000G phase 3, 20 unrelated quality-filtered *A. gambiae* and
20 *A. coluzzii* from Burkina Faso, chromosome 3L, internal coordinates
`[10,000,000, 12,000,000)`. Downloading is limited to the configured sample sets and
region. Local coordinates are recorded with a reversible genomic offset.

```bash
ghostbench empirical-prepare --config configs/empirical.yaml --output results/empirical
ghostbench empirical-run --config configs/empirical.yaml --output results/empirical --model work/empirical-archie/model.rds
```

Before preparation, provide `qc_metadata` in the configuration if the API metadata
does not already include documented `qc_pass` and `unrelated_keep` boolean columns.
This TSV must identify samples by `sample_id`; retain the QC and relatedness source
and selection rule in the report. Sample selection is sorted and deterministic;
unknown QC or relatedness does not silently pass.

Before inference, document mutation rate, recombination rate, effective size,
TRACE timescale, and sources in the empirical configuration. Supply an ancestrally
polarized VCF and set `ancestral_polarization_verified: true` with the supporting
source. The downloaded VCF is reference-oriented and **does not establish ancestral
state**. Retrain ArchIE for 40 target and 40 reference haplotypes using appropriate
simulations. Provide `frozen_thresholds: {archie: ..., trace: ...}` from validation.
These scientific prerequisites are explicit and currently unresolved.

Repeat preparation/inference with `reference_offset: 20` in a copied configuration
and a new output directory for comparison-population sensitivity. If too few samples
pass, the command fails rather than silently reusing related or lower-quality samples.

## Organization and publishing

Core code lives in `src/ghostbench/`; `configs/` describes experiments; `scripts/`
contains setup and HPC entry points; `tests/` verifies scientific bookkeeping;
`notebooks/` is for consuming saved outputs; `reports/` holds reviewed summaries.
Raw genomes, local settings, binaries, environments, and bulk results are ignored.

The intended public repository is `tyy-d/ghost-introgression-benchmark`.
After committing verified code and compact results, `bash scripts/publish.sh`
creates and pushes the repository and checks that remote `main` matches local HEAD.
No command stages files from the enclosing home-directory repository.

## Sources

- [Zhang et al., TRACE](https://doi.org/10.1126/science.aef8874)
- [Durvasula and Sankararaman, ghost introgression in Africans](https://doi.org/10.1126/sciadv.aax5097)
- [ArchIE software](https://github.com/sriramlab/ArchIE)
- [TRACE software](https://github.com/YulinZhang9806/trace)
- [SINGER software and parameter guidance](https://github.com/popgenmethods/SINGER)
- [Ag1000G phased data](https://www.malariagen.net/data_package/ag1000g-phase3-hap/)

## Limits

No whole-genome scan, no new extinct-species identification, and no inference of
accuracy from overlap alone. Ancestral structure is a model stress test, not a
biologically equivalent no-admixture null. Species-specific model calibration,
ancestral polarization, accessibility, reference admixture, and MCMC convergence
are essential parts of interpreting the empirical results.
