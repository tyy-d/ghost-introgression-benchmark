"""Shared configuration, coordinates, and result formats."""

import json
import subprocess
import time
from pathlib import Path

import pandas as pd
import yaml

PREDICTION_COLUMNS = [
    "replicate", "method", "haplotype", "chromosome", "start", "end", "probability"
]


def read_config(path):
    path = Path(path)
    config = yaml.safe_load(path.read_text())
    if "base" in config:
        config = read_config(path.parent / config["base"]) | config
    return config


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_command(command, log, cwd=None):
    """Execute a real tool and retain its complete output and wall time."""
    log = Path(log).resolve()
    log.parent.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with log.open("w") as stream:
        process = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, cwd=cwd)
    elapsed = time.monotonic() - start
    write_json(log.with_suffix(".command.json"), {
        "command": [str(x) for x in command], "cwd": str(cwd) if cwd else None,
        "seconds": elapsed, "returncode": process.returncode,
    })
    if process.returncode:
        raise RuntimeError(f"Command failed ({process.returncode}); see {log}")
    return elapsed


def write_predictions(rows, path):
    frame = pd.DataFrame(rows, columns=PREDICTION_COLUMNS)
    if len(frame):
        if not ((frame.start >= 0) & (frame.end > frame.start)).all():
            raise ValueError("Predictions must use nonempty zero-based half-open intervals")
        if not frame.probability.between(0, 1).all():
            raise ValueError("Ancestry probabilities must lie in [0, 1]")
    frame.to_csv(path, sep="\t", index=False)


def vcf_position(zero_based):
    """VCF is one-based; all internal intervals are zero-based half-open."""
    return int(zero_based) + 1
