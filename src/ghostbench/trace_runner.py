"""Run upstream SINGER and TRACE without access to simulated ancestry labels."""

import os
import sys
import time
from pathlib import Path

import numpy as np
import tskit

from .io import run_command, write_json, write_predictions

ROOT = Path(__file__).resolve().parents[2]


def singer_binary():
    if os.environ.get("SINGER_BINARY"):
        return os.environ["SINGER_BINARY"]
    matches = list((ROOT / ".tools").rglob("singer_master"))
    if not matches:
        raise FileNotFoundError("SINGER missing. Run scripts/bootstrap.sh or set SINGER_BINARY")
    return str(matches[0].resolve())


def infer_arg(vcf, config, output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    prefix = output / "singer"
    vcf = Path(vcf).resolve()
    command = [sys.executable, singer_binary(), "-vcf", str(vcf).removesuffix(".vcf"),
               "-output", str(prefix), "-m", str(config["mutation_rate"]),
               "-Ne", str(config["effective_size"]),
               "-ratio", str(config["recombination_rate"] / config["mutation_rate"]),
               "-polar", "0.99",
               "-start", "0", "-end", str(config["sequence_length"]),
               "-n", str(config["singer_samples"] + config["singer_burnin"]),
               "-thin", str(config["singer_thin"]), "-seed", str(config.get("singer_seed", config["seed"]))]
    # SINGER owns <prefix>.log for recovery state. Capturing stdout there would
    # corrupt its checkpoints, especially when the official wrapper retries.
    seconds = run_command(command, output / "singer.stdout.log", timeout=config.get("job_timeout_seconds"))
    # Different SINGER releases expose tree conversion as a separate tool.
    converter = Path(singer_binary()).parent / "convert_to_tskit"
    # Fresh attempt directories prevent accidental reuse of an older posterior.
    trees = []
    if converter.exists():
        seconds += run_command([
            sys.executable, str(converter), "-input", str(prefix), "-output", str(prefix),
            "-start", str(config["singer_burnin"]),
            "-end", str(config["singer_samples"] + config["singer_burnin"]),
            "-step", "1",
        ], output / "convert.log", timeout=config.get("job_timeout_seconds"))
        trees = sorted(output.glob("*.trees"))
    if len(trees) != config["singer_samples"]:
        raise RuntimeError(f"Expected {config['singer_samples']} retained ARGs, found {len(trees)} at {output}")
    stdout = (output / "singer.stdout.log").read_text(errors="replace")
    write_json(output / "quality.json", {"retained_samples": len(trees), "burnin_samples": config["singer_burnin"],
                                         "thinning": config["singer_thin"], "seed": config.get("singer_seed", config["seed"]),
                                         "recovery_attempts": stdout.count("Auto-debug iteration:"),
                                         "convergence": config.get("convergence_status", "unresolved")})
    return trees, seconds


def predict(trees, metadata, config, output, replicate="0", method="trace", chromosome="1",
            callable_intervals=None):
    from tracehmm import TRACE
    from tracehmm.extract_cli import get_data

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    targets = metadata[metadata.population == "target"]
    nodes = targets.node.astype(int).tolist()
    arrays = []
    start_time = time.monotonic()
    resolution = config.get("trace_resolution", 1000)
    spans = None
    for path in trees:
        ts = tskit.load(path)
        if list(ts.samples()) != metadata.node.astype(int).tolist():
            raise ValueError("ARG leaf order does not match VCF haplotype order")
        ncoal, _, _, _, spans, _, _ = get_data(
            ts, nodes, config["trace_time"], resolution
        )
        arrays.append(ncoal)
    ncoal = np.mean(arrays, axis=0)
    length = int(config["sequence_length"])
    spans = np.array([(left, min(length, left + resolution))
                      for left in range(0, length, resolution)])
    from .evaluation import overlap_length

    accessible = np.ones(len(spans))
    if callable_intervals is not None:
        accessible = np.array([overlap_length(left, right, callable_intervals) > 0
                               for left, right in spans], dtype=float)
    rows, parameters = [], []
    # TRACE defaults to 1e-8/bp. Scale physical spans into its expected map units
    # to use the configured uniform rate while retaining physical output coordinates.
    rate_scale = config["recombination_rate"] / 1e-8
    for index, target in enumerate(targets.itertuples()):
        hmm = TRACE()
        observed = ncoal[index].copy()
        observed[accessible == 0] = 0
        hmm.init_hmm(observed, spans * rate_scale, include_regions=accessible,
                     seed=config["seed"])
        fitted = hmm.train(seed=config["seed"])
        posterior, _, _ = hmm.decode(seed=config["seed"])
        probabilities = np.exp(posterior)[1]
        if not np.isfinite(probabilities).all():
            raise RuntimeError(f"Nonfinite TRACE posterior for {target.haplotype}")
        parameters.append({"haplotype": target.haplotype, "fit": str(fitted)})
        for (left, right), probability, include in zip(spans, probabilities, accessible, strict=True):
            if include:
                rows.append({"replicate": replicate, "method": method,
                             "haplotype": target.haplotype, "chromosome": chromosome,
                             "start": int(left), "end": int(right),
                             "probability": float(probability)})
    write_predictions(rows, output / "predictions.tsv")
    write_json(output / "fits.json", parameters)
    return time.monotonic() - start_time
