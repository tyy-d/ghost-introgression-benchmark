"""Experiment orchestration and recorded provenance."""

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time
from pathlib import Path

import pandas as pd

from . import archie, trace_runner
from .evaluation import metrics, precision_recall, score_grid
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
    for method in ["archie", "trace", "trace_true_arg"]:
        path = output / method / "predictions.tsv"
        if not path.exists():
            continue
        predictions = pd.read_csv(path, sep="\t")
        grid = score_grid(predictions, truth, targets, config["sequence_length"],
                          config["window_size"], callable_intervals)
        grid.to_csv(output / method / "grid.tsv", sep="\t", index=False)
        precision_recall(grid).to_csv(output / method / "precision_recall.tsv", sep="\t", index=False)
        threshold = (thresholds or {}).get(method, 0.5)
        result = {"method": method, **metrics(grid, threshold)}
        results.append(result)
    write_json(output / "metrics.json", results)
    return results


def run(config, output, model, include_true_arg=True, thresholds=None):
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
    write_json(output / "runtime.json", runtimes)
    evaluate(output, config, thresholds)
    from .report import create_report
    create_report(output)


def benchmark(config, output, model):
    """Calibrate once on held-out validation; then freeze across all scenarios."""
    from .evaluation import select_threshold

    output = Path(output)
    validation = config | {"seed": config["validation_seed"],
                           "sequence_length": config["validation_windows"] * config["window_size"]}
    validation_path = output / "validation"
    if not (validation_path / "metrics.json").exists():
        run(validation, validation_path, model)
    thresholds = {method: select_threshold(pd.read_csv(
        validation_path / method / "grid.tsv", sep="\t"
    )) for method in ["archie", "trace", "trace_true_arg"]}
    write_json(output / "frozen_thresholds.json", {
        "validation_seed": config["validation_seed"], "thresholds": thresholds,
        "rule": "callable-weighted F1; stricter threshold wins ties",
    })
    rows = []
    for scenario_index, scenario in enumerate(config["scenarios"]):
        for replicate in range(config["replicates"]):
            experiment = config | scenario | {
                "seed": config["seed"] + 1000 * scenario_index + 2 * replicate
            }
            path = output / scenario["name"] / str(replicate)
            if not (path / "metrics.json").exists():
                run(experiment, path, model, thresholds=thresholds)
            for row in json.loads((path / "metrics.json").read_text()):
                runtime = json.loads((path / "runtime.json").read_text())
                rows.append({"scenario": scenario["name"], "replicate": replicate, **row,
                             "interpretation": "alternative-history call rate" if scenario.get("ancestral_structure")
                             else "simulation accuracy", **runtime})
    pd.DataFrame(rows).to_csv(output / "benchmark.tsv", sep="\t", index=False)
