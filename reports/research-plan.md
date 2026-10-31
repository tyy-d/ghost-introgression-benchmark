# Two-to-three-week research schedule

| Period | Deliverable | Acceptance |
|---|---|---|
| Days 1–4 | Environment, isolated repository, synthetic run | Both real methods finish; labels and coordinates verified; public push verified |
| Days 5–10 | Six-scenario simulation study | Ten independent replicates each, frozen validation thresholds, accuracy and runtime figures |
| Days 11–15 | Mosquito pilot and report | Documented QC, unrelated samples, polarization, biological parameters; empirical overlap and reference sensitivity |

Stop expanding biological scope if environment or ARG inference is slower than
expected. Preserve method execution and scientific bookkeeping as the priorities.
The scientific report must distinguish completed experiments from configuration-only
work and must not turn a short single-seed smoke run into a ranking of methods.

The available compute assumption is laptop development plus a university allocation.
No paid cloud resource is authorized. Scheduler, account and cluster connection are
not assumed; use `scripts/run_hpc.sh` inside an allocated job.

Empirical prerequisites still to establish from primary sources: mutation and
recombination rates, an appropriate effective size and timescale, ancestral allele
orientation, and release-specific QC/relatedness evidence. Population structure and
comparison-panel admixture are part of the interpretation, not nuisance filters
that can be assumed away.
