# Setup and verification status

## Completed

- Created an isolated repository under the BIO course directory; verified that its
  Git root is the project directory, not the enclosing home-directory repository.
- Implemented simulation and migration-based truth extraction, upstream ArchIE
  feature adapters and R GLM training, a SINGER/TRACE runner, callable-weighted
  evaluation, threshold calibration, reporting, and bounded empirical preparation.
- Added fixed smoke and replicated-study configurations, a project-local Conda
  installer, HPC entry point, source commit identifiers, and public-publish script.
- Ran available unit checks with the existing scientific Python environment:
  **11 passed, 2 skipped**. The two skipped checks require `msprime` and cover
  simulation/VCF sample ordering and exact ancestry-label extraction.
- Compiled Python sources for syntax checking and validated shell-script syntax.
- Built upstream SINGER's C++ source into a native Apple Silicon executable and
  completed its official wrapper on a deterministic toy VCF, obtaining two
  posterior samples in approximately 1.34 seconds. The wrapper used one of its
  built-in recovery attempts after the native program initially crashed. This
  component-only check is separate from the required two-method simulation run.

## Not completed or verified

- **No successful ArchIE + SINGER + TRACE smoke run.** Required dependencies could
  not be installed, so method execution, training, exact simulated truth extraction,
  SINGER conversion, and end-to-end figures remain unverified.
- No replicated benchmark, empirical ancestry calls, empirical overlap, or scientific
  accuracy claims. The generated report command will create figures only after a
  successful actual run; there are no fabricated results or mock fallback methods.
- Empirical preparation still requires release-specific QC/relatedness evidence;
  inference also requires species-specific parameters, ancestral polarization, a
  correctly sized trained model, and validation-frozen thresholds.

## Environment blockers

The sandbox cannot resolve external package hosts. Network-enabled Python and R
installation commands were submitted, including a narrower wheels-only download.
Automatic approval review timed out before those installers ran, including after
the user explicitly approved project-local installation. These were review-service
timeouts, not findings that the dependencies were unsafe.

An upstream ArchIE checkout and TRACE checkout were downloaded successfully.
SINGER was retrieved with a source-only sparse checkout and built locally. All
downloads live in ignored directories, not in committed project code.

The existing base Python's pytest debugging plugin segfaulted while importing
`readline`. The available unit checks succeeded with:

```bash
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=src /opt/anaconda3/bin/python \
  -m pytest -p no:debugging -q
```

This is an environment workaround, not evidence that the new isolated environment
has been installed or tested. Ruff is declared but has not been run locally.

## Reproduce the remaining setup

```bash
bash scripts/bootstrap.sh
conda activate "$PWD/.tools/environment"
pytest -q
ruff check src tests
ghostbench smoke --config configs/smoke.yaml --output results/smoke
```

After a successful real run, copy its compact report, metrics, figure and provenance
into a tracked report directory; keep bulk intermediate files ignored.

## GitHub

The intended remote is public `tyy-d/ghost-introgression-benchmark`. The GitHub CLI
recognizes account `tyy-d`, but its authentication check reports an invalid token.
Publication status and exact attempt output are retained in `github-push-attempt.txt`.
Repair authentication with `gh auth login --hostname github.com`, then run
`bash scripts/publish.sh` and verify the printed remote commit identifier.
