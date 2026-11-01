"""Frozen validation calibration and a bounded local worker queue."""

import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from .evaluation import base_curve, fdr_threshold
from .io import terminate_group, write_json
from .registry import METHODS, digest, identity, validate_complete

COMPUTE_KEYS = ["sequence_length", "target_diploids", "reference_diploids", "effective_size",
                "donor_split", "reference_split", "admixture_time", "admixture_fraction",
                "mutation_rate", "recombination_rate", "window_size", "trace_time", "trace_resolution",
                "singer_samples", "singer_burnin", "singer_thin", "ancestral_structure",
                "reference_admixture_fraction", "simulation_effective_size", "convergence_status",
                "pilot_fingerprint"]


def experiment_config(config, role, seed, scenario="baseline", replicate=0, changes=None):
    clean = {key: config[key] for key in COMPUTE_KEYS if key in config}
    clean.update(changes or {})
    clean.update(role=role, seed=seed, scenario=scenario, replicate=replicate)
    return clean


class Budget:
    """Persist elapsed batch time, not the sum of parallel job durations."""
    def __init__(self, root, seconds):
        self.path = Path(root) / "budget.json"
        self.limit = seconds

    def __enter__(self):
        self.previous = json.loads(self.path.read_text()) if self.path.exists() else {"elapsed_seconds": 0}
        if self.previous.get("running"):
            raise ValueError("Previous batch was not reconciled; inspect active jobs before resuming")
        self.started = time.monotonic()
        self.deadline = self.started + max(0, self.limit-self.previous["elapsed_seconds"])
        write_json(self.path, self.previous | {"running": True, "started_at": time.time(), "limit_seconds": self.limit})
        return self

    def remaining(self):
        return max(0, self.deadline-time.monotonic())

    def __exit__(self, kind, error, trace):
        write_json(self.path, {"running": False, "elapsed_seconds": self.previous["elapsed_seconds"]+time.monotonic()-self.started,
                               "limit_seconds": self.limit, "ended_at": time.time()})


def batch(jobs, config, output, budget):
    """Each worker has its own process group; no more than four run at once."""
    output = Path(output)
    specs = output / ".jobs"
    specs.mkdir(parents=True, exist_ok=True)
    pending, active, finished = list(jobs), [], []
    planned_path = output / "planned_jobs.json"
    planned = json.loads(planned_path.read_text()) if planned_path.exists() else {}
    planned.update({job["output"]: job for job in jobs})
    write_json(planned_path, planned)
    workers = min(4, max(1, int(config.get("workers", 4))))
    estimate = float(config.get("estimated_job_seconds", 600))
    try:
        while pending or active:
            while pending and len(active) < workers and budget.remaining() >= estimate:
                job = pending.pop(0)
                name = hashlib.sha256(str(job["output"]).encode()).hexdigest()[:16]
                spec = specs / f"{name}.json"
                write_json(spec, job | {"timeout": budget.remaining()})
                stream = (specs / f"{name}.stdout.log").open("w")
                env = os.environ | {key: "1" for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]}
                process = subprocess.Popen([sys.executable, "-m", "ghostbench.worker", str(spec)],
                                           stdout=stream, stderr=subprocess.STDOUT, env=env, start_new_session=True)
                active.append((process, stream, job))
            for item in list(active):
                process, stream, job = item
                if process.poll() is not None:
                    stream.close()
                    finished.append({"output": job["output"], "returncode": process.returncode})
                    active.remove(item)
            if not active and pending and budget.remaining() < estimate:
                break
            if budget.remaining() <= 0:
                break
            if active:
                time.sleep(0.2)
    finally:
        for process, stream, job in active:
            terminate_group(process)
            stream.close()
            finished.append({"output": job["output"], "returncode": process.poll(), "interrupted": True})
    write_json(output / "last_batch.json", {"finished": finished, "unscheduled": [job["output"] for job in pending]})
    return finished


def calibration(config, output, model, budget):
    validate_seed_ranges(config, model)
    output = Path(output)
    settings = experiment_config(config, "validation", config["validation_seed"])
    fingerprint, _ = identity(settings, model, None)
    freeze_path = output / "frozen_thresholds.json"
    if freeze_path.exists():
        frozen = json.loads(freeze_path.read_text())
        if frozen["identity"] != fingerprint or frozen["target_fdr"] != config.get("target_fdr", 0.05):
            raise ValueError("Frozen calibration is stale; select a new study directory")
        for record in frozen["validation_runs"]:
            _, manifest = validate_complete(record["path"])
            if manifest["fingerprint"] != record["fingerprint"]:
                raise ValueError("Validation provenance changed")
        return frozen
    jobs = []
    for replicate in range(6):
        parameters = experiment_config(config, "validation", config["validation_seed"]+2*replicate, replicate=replicate)
        jobs.append({"config": parameters, "output": str((output / "validation" / str(replicate)).resolve()),
                     "model": str(Path(model).resolve()), "thresholds": None})
    batch(jobs, config, output, budget)
    pools = {method: [] for method in METHODS}
    records = []
    for job in jobs:
        data, manifest = validate_complete(job["output"])
        records.append({"path": job["output"], "fingerprint": manifest["fingerprint"], "seed": job["config"]["seed"]})
        for method in METHODS:
            pools[method].append(pd.read_csv(data / method / "segments.tsv", sep="\t"))
    selected = {}
    for method, frames in pools.items():
        segments = pd.concat(frames, ignore_index=True)
        selected[method] = fdr_threshold(segments, config.get("target_fdr", 0.05))
        base_curve(segments).to_csv(output / f"validation_{method}_curve.tsv", sep="\t", index=False)
    frozen = {"schema": 2, "identity": fingerprint, "source_role": "validation", "target_fdr": config.get("target_fdr", 0.05),
              "validation_runs": records, "operating_points": selected,
              "thresholds": {method: value["threshold"] for method, value in selected.items()},
              "note": "Validation target, not a guarantee of test FDR; thresholds are never selected on test results"}
    write_json(freeze_path, frozen)
    return frozen


def validate_seed_ranges(config, model):
    """Include mutation seeds (+1), not just ancestry seeds, in leakage checks."""
    training = json.loads((Path(model).parent / "training_summary.json").read_text())
    train = set(range(training["seed"], training["seed"]+2*training["windows"]))
    validation = {config["validation_seed"]+offset for offset in range(12)}
    test = {config["seed"]+offset+1000*condition+2*replicate+mutation
            for offset, conditions in [(0, len(config["scenarios"])), (20000, 2)]
            for condition in range(conditions) for replicate in range(10) for mutation in [0, 1]}
    if train & validation or train & test or validation & test:
        raise ValueError("Training, validation and test RNG seeds overlap")


def configured(config, output):
    config = dict(config)
    pilot_path = Path(output) / "pilot" / "pilot.json"
    if not pilot_path.exists():
        raise ValueError("Run the convergence pilot before calibration or research inference")
    report = json.loads(pilot_path.read_text())
    from .convergence import pilot_identity
    if report["fingerprint"] != pilot_identity(config):
        raise ValueError("Pilot source/configuration changed; select a new study directory")
    if report["status"] != "screened":
        raise ValueError("Convergence is unresolved: retain diagnostic-only results, do not silently run the research benchmark")
    config.update(report["selected_settings"])
    config.update(convergence_status="screened", pilot_fingerprint=digest(pilot_path))
    return config


def run_study(config, output, model, stage="core"):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = configured(config, output)
    with Budget(output, config.get("budget_seconds", 43200)) as budget:
        frozen = calibration(config, output, model, budget)
        scenarios = list(config["scenarios"])
        counts = range(5) if stage == "core" else range(5, 10)
        if stage == "stress":
            scenarios = [{"name": "reference_admixed", "reference_admixture_fraction": 0.05},
                         {"name": "demography_mismatch", "simulation_effective_size": 2*config["effective_size"]}]
            counts = range(5)
        jobs = []
        for replicate in counts:
            for index, scenario in enumerate(scenarios):
                offset = 20000 if stage == "stress" else 0
                seed = config["seed"] + offset + 1000*index + 2*replicate
                changes = {key: value for key, value in scenario.items() if key != "name"}
                parameters = experiment_config(config, "test", seed, scenario["name"], replicate, changes)
                jobs.append({"config": parameters,
                             "output": str(output / "test" / scenario["name"] / str(replicate)),
                             "model": str(Path(model).resolve()), "thresholds": frozen["thresholds"]})
        result = batch(jobs, config, output, budget)
    from .analysis import summarize
    summarize(output)
    return result


def calibrate_only(config, output, model):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = configured(config, output)
    with Budget(output, config.get("budget_seconds", 43200)) as budget:
        return calibration(config, output, model, budget)
