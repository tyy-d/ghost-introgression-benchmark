"""Content-addressed run manifests, immutable completed attempts, and integrity checks."""

import hashlib
import importlib.metadata
import json
import time
from pathlib import Path

from .io import write_json

ROOT = Path(__file__).resolve().parents[2]
CORE_FILES = ["simulation.py", "archie.py", "trace_runner.py", "evaluation.py", "pipeline.py", "io.py", "process.py", "registry.py"]
METHODS = ["archie", "trace", "trace_true_arg"]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            value.update(block)
    return value.hexdigest()


def identity(config, model, thresholds):
    import subprocess

    from .archie import rscript
    from .trace_runner import singer_binary

    versions = {}
    for package in ["numpy", "scipy", "msprime", "tskit", "pandas", "scikit-learn", "tracehmm", "psutil"]:
        versions[package] = importlib.metadata.version(package)
    r = subprocess.run([rscript(), "--version"], capture_output=True, text=True, check=True)
    wrapper = Path(singer_binary())
    source = {name: digest(ROOT / "src/ghostbench" / name) for name in CORE_FILES}
    source["archie_glm.R"] = digest(ROOT / "scripts/archie_glm.R")
    payload = {"schema": 2, "config": config, "thresholds": thresholds, "software": versions,
               "source": source, "model_sha256": digest(model),
               "model_metadata_sha256": digest(Path(model).parent / "training_summary.json"),
               "singer_sha256": digest(wrapper.parent / "singer"), "wrapper_sha256": digest(wrapper),
               "r_version": (r.stdout+r.stderr).strip()}
    fingerprint = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return fingerprint, payload


def active_dir(root):
    root = Path(root)
    manifest = root / "manifest.json"
    if not manifest.exists():
        return root
    return root / json.loads(manifest.read_text())["artifact_dir"]


def validate_complete(root):
    root = Path(root)
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest["status"] not in {"complete", "diagnostic-only"}:
        raise ValueError(f"Run is not complete: {root} ({manifest['status']})")
    data = root / manifest["artifact_dir"]
    if not data.resolve().is_relative_to(root.resolve()):
        raise ValueError("Run artifact directory escapes its root")
    for filename, expected in manifest["artifact_hashes"].items():
        path = data / filename
        if not path.resolve().is_relative_to(data.resolve()):
            raise ValueError("Artifact path escapes the run directory")
        if not path.is_file() or digest(path) != expected:
            raise ValueError(f"Missing or altered result artifact: {path}")
    return data, manifest


class Run:
    def __init__(self, root, config, model, thresholds, diagnostic=False):
        self.root = Path(root).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.path = self.root / "manifest.json"
        self.fingerprint, self.payload = identity(config, model, thresholds)
        self.diagnostic = diagnostic
        self.cached = False

    def __enter__(self):
        previous = json.loads(self.path.read_text()) if self.path.exists() else None
        if previous:
            if previous["fingerprint"] != self.fingerprint:
                raise ValueError(f"Stale run fingerprint at {self.root}; use a new output directory")
            if previous["status"] in {"complete", "diagnostic-only"}:
                self.data, self.manifest = validate_complete(self.root)
                self.cached = True
                return self
        elif any(self.root.iterdir()):
            raise ValueError(f"Unmanaged output directory cannot be resumed: {self.root}")
        attempt = previous.get("attempt", 0) + 1 if previous else 1
        self.data = self.root / "attempts" / f"{attempt:04d}"
        self.data.mkdir(parents=True)
        self.manifest = {"schema": 2, "fingerprint": self.fingerprint, "identity": self.payload,
                         "status": "running", "attempt": attempt,
                         "artifact_dir": str(self.data.relative_to(self.root)),
                         "started_at": time.time(), "role": self.payload["config"].get("role", "diagnostic"),
                         "scenario": self.payload["config"].get("scenario", "baseline"),
                         "replicate": self.payload["config"].get("replicate", 0)}
        write_json(self.path, self.manifest)
        self.started = time.monotonic()
        return self

    def __exit__(self, kind, error, traceback):
        if self.cached:
            return False
        self.manifest["elapsed_seconds"] = time.monotonic() - self.started
        self.manifest["ended_at"] = time.time()
        if kind is not None:
            self.manifest["status"] = "interrupted" if issubclass(kind, (KeyboardInterrupt, SystemExit, TimeoutError)) else "failed"
            self.manifest["error"] = f"{kind.__name__}: {error}"
        else:
            required = ["config.json", "runtime.json", "base_metrics.json", "provenance.json",
                        "inputs/sequences.vcf", "inputs/samples.tsv", "inputs/callable.tsv",
                        "truth/ancestry.tsv", "truth/ancestry.trees"]
            required += [f"{method}/{name}" for method in METHODS for name in ["predictions.tsv", "segments.tsv", "grid.tsv", "base_curve.tsv"]]
            try:
                self.manifest["artifact_hashes"] = {name: digest(self.data / name) for name in required}
                self.manifest["input_hashes"] = {name: self.manifest["artifact_hashes"][name] for name in required if name.startswith("inputs/")}
                self.manifest["status"] = "diagnostic-only" if self.diagnostic else "complete"
            except Exception as failure:
                self.manifest["status"] = "failed"
                self.manifest["error"] = str(failure)
                write_json(self.path, self.manifest)
                raise
        write_json(self.path, self.manifest)
        return False
