# Ghost Introgression Benchmark

A reproducible comparison of **ArchIE** and **SINGER + TRACE** for ancestry from an
unsampled donor. The project separates research evidence, diagnostic runs, and
interactive threshold previews.

**Milestones:** minimal presentable version and written professor/TA approval by
**November 14**; final readiness targeted for **December 7**. The assignment states
**December 8** as the submission deadline. See [course requirements](reports/course-requirements.md).

## Current deliverables

- [Approval proposal](reports/approval-proposal.md), ready for student review and submission.
- [November MVP report](reports/nov14-mvp/research-report-draft.pdf), with five text
  pages plus figures, 12-point type and 1.5 spacing. It is a draft, not an approval record.
- A read-only explorer with performance, ancestry-track and run-quality views.
- A portable processed-results demonstration in `reports/nov14-mvp/explorer-results/`.
- Exact callable-base scoring, validation-only FDR calibration, paired replicate
  bootstrapping, convergence screening, integrity checks, and bounded job execution.

The original [smoke report](reports/smoke/report.md) is preserved as historical
functionality evidence. Its perfect window-level scores excluded mixed regions.
**They are not the current primary evaluation.** The new diagnostic run includes
all callable bases and reveals errors at the unfrozen 0.5 threshold. Neither that
single seed nor an explorer threshold preview establishes comparative accuracy.

## Open the explorer

From this project directory, using the already installed local environment:

```bash
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
.venv/bin/ghostbench view --output results/research-v2
```

Open the printed localhost URL. For a fresh checkout without raw experiment output,
use the portable demonstration instead:

```bash
.venv/bin/ghostbench view --output reports/nov14-mvp/explorer-results
```

The explorer never launches jobs or modifies frozen thresholds. It filters saved
results, validates hashes, links genomic x-axes, and exports PNG/SVG/PDF plus TSV and
selection metadata. Start/End controls define the export range; transient mouse
zoom is not saved. Threshold previews are labeled **EXPLORATORY PREVIEW** and are
excluded from research summaries. Validation results are not displayed as test evidence.

## Install a clean environment

Requires Git, Conda, and a C++17 compiler:

```bash
bash scripts/bootstrap.sh
conda activate "$PWD/.tools/environment"
```

The bootstrap pins upstream source commits and installs Python, R and explorer
requirements inside the project. Empirical libraries are optional:
`pip install -e '.[empirical]'`. The current verified checkout instead uses `.venv/`
and `.tools/r/`; executable discovery supports both arrangements. `SINGER_BINARY`
points to the official Python wrapper `singer_master`, not the C++ core executable.
`RSCRIPT` can override R discovery. Package versions and upstream commits are recorded.

The inherited Python/readline sandbox requires this test invocation:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 MPLCONFIGDIR=.cache/matplotlib .venv/bin/python \
  -m pytest -p no:debugging -p no:capture -q
.venv/bin/ruff check src tests scripts
```

## Execute the research protocol

```bash
# Always use a new directory for a genuinely changed scientific configuration.
.venv/bin/ghostbench pilot --config configs/research.yaml --output results/research-v2
.venv/bin/ghostbench calibrate --config configs/research.yaml --output results/research-v2
.venv/bin/ghostbench benchmark --config configs/research.yaml --output results/research-v2 --stage core
.venv/bin/ghostbench summarize --output results/research-v2
.venv/bin/ghostbench export --output results/research-v2 --destination reports/research-v2
.venv/bin/ghostbench paper --output results/research-v2 --destination reports/nov14-mvp
```

The first long pilot was interrupted and remains unscreened; calibration and the
replicated study have not run. Before retrying this local campaign, archive the
existing `results/research-v2/pilot/` directory under a dated history directory,
then rerun the pilot with the same study root and budget ledger. Its source
fingerprint is preserved and must not be rewritten. See [execution status](reports/setup-status.md).

The pilot must finish with screened settings before calibration or research
inference proceeds. Unresolved convergence is reported rather than silently ignored.
The pilot runs two independent chains for baseline and null datasets, trying
50/50, 100/100 and 200/200 retained/burn-in samples, thinning five. Screening uses
classical split R-hat <=1.05 and ESS >=50 for scalar statistics; constants alone
cannot prove mixing, and a passing pilot does not prove convergence in every test run.

Calibration uses **six independent 1-Mb baseline regions**. Each method maximizes
validation recall subject to estimated base-level FDR <=5%, preferring stricter
thresholds on ties. A method with no qualifying nonempty call set abstains. Thresholds
are frozen before testing. Report achieved test FDR; a validation target is not a
guarantee under changed histories or ancestry prevalence.

Core scheduling interleaves conditions and initially targets five replicates per
condition. Only after that stage and within the remaining budget:

```bash
.venv/bin/ghostbench benchmark --output results/research-v2 --stage expansion
.venv/bin/ghostbench benchmark --output results/research-v2 --stage stress
```

Expansion adds replicates 5-9; stress tests add five replicates each of reference
admixture and demographic mismatch. The ArchIE model remains fixed. Study execution
allows at most four workers, one numerical-library thread per worker, and twelve
elapsed batch hours. The pilot is included in that budget and separately capped
at two hours. No new work is scheduled when the remaining budget is below the
estimated job duration. Keep a single study directory for the campaign so elapsed
budgets carry across batches. Changing directories must not be used to evade the budget.

## Evaluation and interpretation

Primary metrics use exact intersections of called and true ancestry at each
method's native resolution over **all callable bases**, including mixed regions.
TP, FP, FN and TN are base counts. Common 50-kb scores remain secondary diagnostics.
Undefined FDR/precision means no called bases, not perfect FDR control. The pure-null
condition is summarized by false-called callable fraction and false-called kb/Mb.
The structure scenario includes ancestral migration: report its call rate, not an
ordinary false-positive rate against a biologically unadmixed truth.

Confidence intervals resample independent simulation replicates with 2,000 bootstrap
draws and seed 20261101. Paired differences use the same replicate IDs for both
methods. With fewer than five usable pairs, results remain descriptive without
confidence intervals. Unscheduled, failed, interrupted, corrupted and diagnostic-only
runs remain visible. The displayed illustrative research cases are selected by the
largest FP/FN region with deterministic tie-breaking, not by visual appeal.

## Run integrity and data layout

Each managed run stores a manifest and fresh numbered attempt directories. Its
identity covers configuration, RNG seeds, input/model/source/software hashes and
thresholds. Only matching, completed, hash-validated results are cached. Changed
configurations and altered/missing output are rejected. Incomplete attempts are
not reused as successful results. Missing method files are never scored as zero ancestry.

- `inputs/`: phased biallelic VCF, ordered sample metadata, callable intervals.
- `truth/`: exact migration-based ancestry and true ARGs, separate from detector inputs.
- Method directories: predictions, exact scoring segments, PR curves and secondary grids.
- Run metadata: configuration, versions, wall/CPU time and memory availability notes.

Prediction TSV columns are `replicate`, `method`, `haplotype`, `chromosome`, `start`,
`end`, `probability`; intervals are **zero-based half-open**, unlike one-based VCF
positions. Genome length and generation-based times are explicitly labeled.

The explorer bundle contains processed intervals only, with a distinct fingerprint.
It cannot be resumed as an inference run. Raw genomes, ARGs, models, environments,
local settings and bulk results remain ignored. The historical smoke report is not overwritten.

## Empirical readiness and course submission

```bash
.venv/bin/ghostbench readiness --config configs/empirical.yaml --output results/empirical-readiness
```

The bounded metadata check encountered HTTP 403 responses from the documented GCS
URLs. This is an access result, not evidence that the public dataset does not exist.
Empirical inference remains gated on QC/relatedness evidence, ancestral polarization,
species parameters, a correctly sized ArchIE model, and validation-frozen thresholds.
Do not treat reference alleles as ancestral or the controlled software parameters
as fitted mosquito biology. Optional biological observations would be exploratory.

The final paper should follow the course's approximate two-page introduction,
2-3 results pages excluding figures, and 1-2 discussion pages, plus annotated
executed code. The student must obtain written approval by November 14. This
software does not send messages to instructors or claim approval was obtained.

## Sources

- [Durvasula and Sankararaman, 2020](https://doi.org/10.1126/sciadv.aax5097)
- [Zhang et al., TRACE, 2026](https://doi.org/10.1126/science.aef8874)
- [ArchIE source](https://github.com/sriramlab/ArchIE)
- [TRACE source](https://github.com/YulinZhang9806/trace)
- [SINGER source and convergence guidance](https://github.com/popgenmethods/SINGER)
- [Ag1000G public phased data](https://www.malariagen.net/data_package/ag1000g-phase3-hap/)
- [Anopheles mutation-rate study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8742016/)
