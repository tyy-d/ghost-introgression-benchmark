"""Use ArchIE's published features and unpenalized binomial GLM.

    Compatibility fixes: dynamic feature count, finite constant-site statistics,
    and consistent SNP positions after reference/focal filtering for S*.
"""

import importlib.util
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .evaluation import overlap_length
from .io import run_command, write_json, write_predictions

ROOT = Path(__file__).resolve().parents[2]


def upstream():
    path = ROOT / "vendor/ArchIE/data/calc_stats_window_data.py"
    if not path.exists():
        raise FileNotFoundError("ArchIE not installed. Run scripts/bootstrap.sh")
    spec = importlib.util.spec_from_file_location("archie_features", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def s_star(haplotypes, positions, focal):
    """Vectorized ArchIE recurrence, preserving physical positions after filtering."""
    keep = haplotypes[focal] == 1
    haplotypes, positions = haplotypes[:, keep], positions[keep]
    count = len(positions)
    if count < 2:
        return 0.0
    best = np.zeros(count)
    frequencies = haplotypes.sum(axis=0)
    for i in range(1, count):
        distance = np.abs(haplotypes[:, :i] - haplotypes[:, [i]]).sum(axis=0)
        gap = positions[i] - positions[:i]
        valid = (gap >= 10) & (gap <= 50000) & (frequencies[:i] > 1) & (frequencies[i] > 1)
        score = np.where(distance == 0, gap + 5000, np.where(distance <= 5, -10000, -np.inf))
        best[i] = max(0.0, float(np.max(np.where(valid, best[:i] + score, -np.inf))))
    return float(best.max())


def features(genotypes, positions, target_haplotypes, window_start, window_end):
    """One feature row per focal haplotype; target and reference columns stay ordered."""
    module = upstream()
    keep = (positions >= window_start) & (positions < window_end)
    genotype = genotypes[keep]
    positions = positions[keep]
    target = genotype[:, :target_haplotypes]
    reference = genotype[:, target_haplotypes:]
    if reference.shape[1] == 0:
        raise ValueError("ArchIE requires a comparison population")
    private = ~reference.any(axis=1)
    private_haplotypes = target[private].T
    private_positions = positions[private]
    result = []
    for focal in range(target_haplotypes):
        if not len(target):
            result.append(np.zeros(2 * target_haplotypes + 8))
            continue
        target_list = target.tolist()
        row = module.N_ton(target_list, target_haplotypes, focal)
        row += module.distance_vector(target_list, focal)
        row += [module.min_dist_ref(target_list, reference.tolist(), focal)]
        row += [s_star(private_haplotypes, private_positions, focal)]
        row += [float(private_haplotypes[focal].sum())]
        result.append(np.nan_to_num(row, nan=0.0, posinf=0.0, neginf=0.0))
    return np.asarray(result)


def training_window(task):
    from .simulation import donor_intervals, simulate

    config, index, offset = task
    seed = config["training_seed"] + offset + 2 * index
    fraction = config["admixture_fraction"] if index % 2 else 0.0
    ts = simulate(config, seed, length=config["window_size"], fraction=fraction)
    count = 2 * config["target_diploids"]
    matrix = features(ts.genotype_matrix(), np.array(ts.tables.sites.position), count,
                      0, config["window_size"])
    truth = donor_intervals(ts, config["admixture_time"], range(count))
    labels = []
    for node in range(count):
        proportion = overlap_length(0, config["window_size"], truth[node]) / config["window_size"]
        labels.append(1 if proportion > 0.7 else 0 if proportion < 0.3 else -1)
    return matrix, labels


def rscript():
    local = ROOT / ".tools/r/bin/Rscript"
    return os.environ.get("RSCRIPT", str(local) if local.exists() else "Rscript")


def train(config, output):
    started = time.monotonic()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    matrices, labels = [], []
    tasks = [(config, i, 0) for i in range(config["training_windows"])]
    workers = min(config.get("training_jobs", 4), os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for index, (matrix, label) in enumerate(pool.map(training_window, tasks, chunksize=5)):
            matrices.append(matrix)
            labels.extend(label)
            if index % 100 == 0:
                print(f"ArchIE training windows: {index}/{len(tasks)}", flush=True)
    matrix, labels = np.vstack(matrices), np.asarray(labels)
    eligible = labels != -1
    frame = pd.DataFrame(matrix[eligible], columns=[f"x{i}" for i in range(matrix.shape[1])])
    frame["label"] = labels[eligible]
    frame.to_csv(output / "training.tsv", sep="\t", index=False)
    seconds = run_command([
        rscript(), str(ROOT / "scripts/archie_glm.R"), "train",
        str(output / "training.tsv"), str(output / "model.rds"),
    ], output / "train.log")
    write_json(output / "training_summary.json", {
        "windows": len(tasks), "labeled_haplotypes": int(eligible.sum()),
        "positive_haplotypes": int((labels == 1).sum()), "seed": config["training_seed"],
        "features": matrix.shape[1], "glm_seconds": seconds,
        "total_seconds": time.monotonic() - started,
        "target_haplotypes": 2 * config["target_diploids"],
        "reference_haplotypes": 2 * config["reference_diploids"],
        "window_size": config["window_size"],
    })


def predict(vcf, metadata, config, model, output, replicate="0", chromosome="1",
            start=0, end=None):
    from .vcf import read_vcf

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    positions, genotypes, names = read_vcf(vcf)
    expected = metadata.drop_duplicates("individual").individual.tolist()
    if names != expected:
        raise ValueError("VCF individual order does not match sample metadata")
    count = int((metadata.population == "target").sum())
    training = json.loads((Path(model).parent / "training_summary.json").read_text())
    reference_count = int((metadata.population == "reference").sum())
    if (training.get("target_haplotypes") != count
            or training.get("reference_haplotypes") != reference_count
            or training.get("window_size") != config["window_size"]):
        raise ValueError("ArchIE model sample counts/window size do not match input; retrain the model")
    end = end or config["sequence_length"]
    matrix, intervals = [], []
    for left in range(start, end, config["window_size"]):
        right = min(end, left + config["window_size"])
        matrix.append(features(genotypes, positions, count, left, right))
        intervals.extend([(f"h{hap}", left, right) for hap in range(count)])
    matrix = np.vstack(matrix)
    pd.DataFrame(matrix, columns=[f"x{i}" for i in range(matrix.shape[1])]).to_csv(
        output / "features.tsv", sep="\t", index=False
    )
    seconds = run_command([
        rscript(), str(ROOT / "scripts/archie_glm.R"), "predict",
        str(output / "features.tsv"), str(Path(model).resolve()),
        str((output / "probabilities.tsv").resolve()),
    ], output / "predict.log")
    probabilities = pd.read_csv(output / "probabilities.tsv", sep="\t").probability
    rows = [{"replicate": replicate, "method": "archie", "haplotype": hap,
             "chromosome": chromosome, "start": left, "end": right, "probability": float(p)}
            for (hap, left, right), p in zip(intervals, probabilities, strict=True)]
    write_predictions(rows, output / "predictions.tsv")
    return seconds
