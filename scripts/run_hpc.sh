#!/usr/bin/env bash
# Scheduler-independent entry point. Invoke inside an allocated university job;
# no account, partition, scheduler, or memory allocation is silently assumed.
set -euo pipefail
cd "$(dirname "$0")/.."
export MPLCONFIGDIR="$PWD/.cache/matplotlib"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
ghostbench benchmark --config configs/benchmark.yaml --output results/benchmark --model work/archie/model.rds
