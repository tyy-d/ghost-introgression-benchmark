"""Experiment orchestration and recorded provenance."""

import hashlib
import importlib.metadata
import json
import platform
import resource
import subprocess
import time
from pathlib import Path

import pandas as pd

from . import archie, trace_runner
from .evaluation import (
    base_curve,
    base_metrics,
    base_segments,
    metrics,
    precision_recall,
    score_grid,
)
from .io import write_json
from .simulation import generate

ROOT = Path(__file__).resolve().parents[2]


def provenance(output):
    commits = {}
    for name in ["ArchIE", "trace", "SINGER", "SINGER-source"]:
        path = ROOT / "vendor" / name
        if (path / ".git").exists():
            result = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                                    capture_output=True, text=True)
            if result.returncode == 0:
                commits[name] = result.stdout.strip()
    versions = {package.metadata["Name"]: package.version
                for package in importlib.metadata.distributions()}
    source_hash = hashlib.sha256()
    for path in sorted((ROOT / "src").rglob("*.py")):
        source_hash.update(str(path.relative_to(ROOT)).encode())
        source_hash.update(path.read_bytes())
    native = Path(trace_runner.singer_binary()).parent / "singer"
    r_version = subprocess.run([archie.rscript(), "--version"], capture_output=True, text=True)
    write_json(Path(output) / "provenance.json", {
        "python": platform.python_version(), "platform": platform.platform(),
        "upstream_commits": commits, "packages": versions,
        "project_source_sha256": source_hash.hexdigest(),
        "singer_binary_sha256": hashlib.sha256(native.read_bytes()).hexdigest(),
        "r_version": (r_version.stdout + r_version.stderr).strip(),
    })


def evaluate(output, config, thresholds=None):
    output = Path(output)
    metadata = pd.read_csv(output / "inputs/samples.tsv", sep="\t")
    truth = pd.read_csv(output / "truth/ancestry.tsv", sep="\t")
    callable_frame = pd.read_csv(output / "inputs/callable.tsv", sep="\t")
    callable_intervals = list(zip(callable_frame.start, callable_frame.end))
    targets = metadata.query("population == 'target'").haplotype.tolist()
    results = []
    primary = []
    for method in ["archie", "trace", "trace_true_arg"]:
        path = output / method / "predictions.tsv"
        if not path.exists():
            raise FileNotFoundError(f"Missing method output cannot be scored as no ancestry: {path}")
        predictions = pd.read_csv(path, sep="\t")
        grid = score_grid(predictions, truth, targets, config["sequence_length"],
                          config["window_size"], callable_intervals)
        grid.to_csv(output / method / "grid.tsv", sep="\t", index=False)
        precision_recall(grid).to_csv(output / method / "precision_recall.tsv", sep="\t", index=False)
        threshold = (thresholds or {}).get(method, 0.5)
        segments = base_segments(predictions, truth, targets, callable_frame, require_dense=True)
        segments.to_csv(output / method / "segments.tsv", sep="\t", index=False)
        base_curve(segments).to_csv(output / method / "base_curve.tsv", sep="\t", index=False)
        interpretation = "alternative-history call rate" if config.get("ancestral_structure") else "simulation accuracy"
        row = {"method": method, "interpretation": interpretation,
               "threshold_status": "frozen" if thresholds is not None else "diagnostic-unfrozen",
               **base_metrics(segments, threshold)}
        if config.get("ancestral_structure"):
            for key in ["recall", "precision", "fdr", "f1", "false_called_fraction", "false_called_kb_per_mb"]:
                row[key] = None
        primary.append(row)
        result = {"method": method, **metrics(grid, threshold if threshold is not None else 2.0)}
        results.append(result)
    write_json(output / "metrics.json", results)
    write_json(output / "base_metrics.json", primary)
    return primary


def _execute_run(config, output, model, include_true_arg=True, thresholds=None):
    before = resource.getrusage(resource.RUSAGE_SELF)
    children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    generate(config, output)
    provenance(output)
    metadata = pd.read_csv(output / "inputs/samples.tsv", sep="\t")
    runtimes = {}
    training_summary = json.loads((Path(model).parent / "training_summary.json").read_text())
    runtimes["archie_training_seconds_one_time"] = training_summary["total_seconds"]
    started = time.monotonic()
    archie.predict(output / "inputs/sequences.vcf", metadata, config, model,
                   output / "archie", replicate=str(config["seed"]))
    runtimes["archie_inference_seconds"] = time.monotonic() - started
    trees, runtimes["singer_seconds"] = trace_runner.infer_arg(
        output / "inputs/sequences.vcf", config, output / "singer"
    )
    runtimes["trace_hmm_seconds"] = trace_runner.predict(
        trees, metadata, config, output / "trace", replicate=str(config["seed"])
    )
    runtimes["trace_end_to_end_seconds"] = (
        runtimes["singer_seconds"] + runtimes["trace_hmm_seconds"]
    )
    if include_true_arg:
        runtimes["trace_true_arg_seconds"] = trace_runner.predict(
            [output / "truth/ancestry.trees"], metadata, config,
            output / "trace_true_arg", replicate=str(config["seed"]), method="trace_true_arg",
        )
    after = resource.getrusage(resource.RUSAGE_SELF)
    children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    scale = 1 if platform.system() == "Darwin" else 1024
    runtimes.update(worker_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                    children_cpu_seconds=children_after.ru_utime+children_after.ru_stime-children_before.ru_utime-children_before.ru_stime,
                    worker_peak_rss_bytes=after.ru_maxrss*scale,
                    largest_child_peak_rss_bytes=children_after.ru_maxrss*scale)
    write_json(output / "runtime.json", runtimes)
    evaluate(output, config, thresholds)
    from .report import create_report
    create_report(output)


def run(config, output, model, include_true_arg=True, thresholds=None):
    from .registry import Run

    if not include_true_arg:
        raise ValueError("Research protocol requires the separately labeled true-ARG diagnostic")
    with Run(output, config, model, thresholds,
             diagnostic=config.get("role", "diagnostic") == "diagnostic") as attempt:
        if not attempt.cached:
            _execute_run(config, attempt.data, model, include_true_arg, thresholds)
    return attempt.data


def benchmark(config, output, model):
    """Calibrate once on held-out validation; then freeze across all scenarios."""
    from .study import run_study

    return run_study(config, output, model)
