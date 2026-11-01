# Setup and verification status

## November 1 development update

The November MVP now provides exact callable-base scoring, a read-only Streamlit
explorer, static scientific exports, an approval proposal, and a course-format
draft. All 35 tests pass; the earlier 15-test record below is historical.
The corrected deadlines are November 14 for written approval and a presentable
version, December 7 for our readiness target, and December 8 in the assignment.

A fresh actual-method diagnostic includes mixed-ancestry bases and exposes errors
hidden by the historical unambiguous-window scores. It remains a single-seed,
unfrozen-threshold functionality test, not evidence of calibrated accuracy. The
latest rerun required one recorded SINGER recovery attempt; native recovery can
change the posterior, and short-chain outputs are not treated as stable accuracy estimates.

The two-chain convergence pilot was manually interrupted after 1,405 elapsed
seconds. A sandbox-denied redundant cleanup signal caused its saved failure status;
a process-table audit confirmed no inference processes remained. The runner now
preserves interruption/timeout status and records cleanup errors separately.
Calibration, replicated research experiments, and empirical inference are pending.
The failed pilot is retained with its original source fingerprint. Before retrying,
archive `results/research-v2/pilot/` outside that path, then use the same study root
so its elapsed budget remains charged. Do not overwrite the old pilot record.

The final report draft is provisional. Its empirical gate remains unresolved,
including metadata retrieval HTTP 403, QC/relatedness, polarization, and parameters.

## Historical initial setup

The isolated project repository is published publicly at
[tyy-d/ghost-introgression-benchmark](https://github.com/tyy-d/ghost-introgression-benchmark).
The initial remote commit was verified against local HEAD; the final push includes
the tested pipeline and compact synthetic results.

## Completed

- Installed project-local Python dependencies and R, compiled native Apple Silicon
  SINGER, and installed upstream TRACE with its compiled extension.
- Trained ArchIE on 2,000 independent 50-kb windows: 39,401 eligible haplotype
  examples, including 707 positives. One-time training took approximately 82 s.
- Completed both actual methods on the same seed-42, 1-Mb simulated sequences:
  10 diploid target and 10 diploid reference samples, with the donor unsampled.
- Retained four SINGER posterior ARGs after discarding four initial samples,
  then ran upstream TRACE. This functional run is deliberately short and does
  not establish MCMC convergence. The verified run required no native retries.
- Ran the separately labeled true-ARG TRACE diagnostic and exported metrics,
  precision-recall curves, callable-weighted grids, timings, tool provenance,
  and a visually checked figure to [the saved report](smoke/report.md).
- All 15 tests pass without skips. Ruff and shell syntax checks pass.
- Both actual methods classified the 392 unambiguous grid cells correctly;
  eight mixed-ancestry cells were reported separately. There were only eight
  positive cells, so perfect single-seed scores are not general accuracy estimates.

## Remaining research work

The six-scenario replicated study and the mosquito empirical pilot have not been
run. The benchmark configuration preserves longer MCMC settings for university HPC.
The empirical pipeline requires documented QC/relatedness evidence, ancestral
allele orientation, species-specific parameters, and validation-frozen thresholds.
It rejects missing prerequisites and models trained at incompatible sample counts.
Empirical overlap is agreement, not an accuracy measurement.

## Environment notes

Early network approval attempts timed out, but subsequent installation and GitHub
publication succeeded. The earlier network-disabled GitHub attempt is retained in
`github-push-attempt.txt` as diagnostic history; it is not the current status.

The inherited Conda Python readline module crashes under this sandbox when pytest
loads interactive debugging/capture support. The complete passing test command was:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 MPLCONFIGDIR=.cache/matplotlib .venv/bin/python \
  -m pytest -p no:debugging -p no:capture -q
```

The verified local environment uses `.venv/` and `.tools/r/`. A clean installation
can instead use `scripts/bootstrap.sh` and the project-local Conda environment.
Python package versions are in `requirements.lock.txt`; R session details and
upstream source commits are saved with the reports/configuration. Bulk genomes,
ARGs, environments and models remain ignored.
