"""Shared configuration, coordinates, and result formats."""

import json
import os
import signal
import subprocess
import tempfile
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
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
        temporary = stream.name
    os.replace(temporary, path)


def run_command(command, log, cwd=None, timeout=None):
    """Execute a real tool and retain its complete output and wall time."""
    log = Path(log).resolve()
    log.parent.mkdir(parents=True, exist_ok=True)
    from .process import execute

    seconds, returncode, telemetry = execute(command, log, cwd, timeout)
    write_json(log.with_suffix(".command.json"), telemetry)
    if returncode:
        raise RuntimeError(f"Command failed ({returncode}); see {log}")
    return seconds


def terminate_group(process):
    """Terminate the complete job, including native grandchildren."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        pass
    # The leader can exit before a grandchild: always kill remaining group members.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


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
