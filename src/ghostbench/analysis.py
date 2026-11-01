"""One analysis layer for CLI summaries, the explorer, and figure exports."""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .evaluation import base_metrics
from .io import write_json
from .registry import validate_complete

COUNT_KEYS = ["tp_bases", "fp_bases", "fn_bases", "tn_bases"]
PERFORMANCE = ["recall", "precision", "fdr", "false_called_fraction", "false_called_kb_per_mb", "call_rate"]


def resampled_values(frame, indices):
    counts = frame[COUNT_KEYS].to_numpy(dtype=float)[indices].sum(axis=1)
    tp, fp, fn, tn = counts.T
    total = counts.sum(axis=1)

    def ratio(a, b):
        return np.divide(a, b, out=np.full_like(a, np.nan), where=b != 0)

    return {"recall": ratio(tp, tp+fn), "precision": ratio(tp, tp+fp), "fdr": ratio(fp, tp+fp),
            "false_called_fraction": ratio(fp, total), "false_called_kb_per_mb": ratio(1000*fp, total),
            "call_rate": ratio(tp+fp, total)}


def pooled(frame):
    tp, fp, fn, tn = frame[COUNT_KEYS].sum().to_numpy(dtype=float)
    total = tp+fp+fn+tn
    return {"tp_bases": tp, "fp_bases": fp, "fn_bases": fn, "tn_bases": tn,
            "recall": tp/(tp+fn) if tp+fn else None,
            "precision": tp/(tp+fp) if tp+fp else None, "fdr": fp/(tp+fp) if tp+fp else None,
            "false_called_fraction": fp/total if total else None,
            "false_called_kb_per_mb": 1000*fp/total if total else None,
            "call_rate": (tp+fp)/total if total else None,
            "callable_bases": total, "called_bases": tp+fp, "true_bases": tp+fn}


def collect(root):
    """Corruption and incomplete runs become visible status rows, never zero calls."""
    root = Path(root)
    runs, metrics = [], []
    for path in sorted(root.rglob("manifest.json")):
        manifest = json.loads(path.read_text())
        row = {"run_id": str(path.parent.relative_to(root)), "path": str(path.parent.resolve()),
               "status": manifest["status"], "role": manifest.get("role"),
               "scenario": manifest.get("scenario"), "replicate": manifest.get("replicate"),
               "error": manifest.get("error", ""), "fingerprint": manifest["fingerprint"]}
        configuration = manifest.get("identity", {}).get("config", {})
        row["eligible_test"] = (row["role"] == "test" and row["status"] == "complete"
                                and configuration.get("convergence_status") == "screened"
                                and manifest.get("identity", {}).get("thresholds") is not None)
        if manifest["status"] in {"complete", "diagnostic-only"}:
            try:
                data, checked = validate_complete(path.parent)
                row["artifact_dir"] = str(data)
                runtime = json.loads((data / "runtime.json").read_text())
                for result in json.loads((data / "base_metrics.json").read_text()):
                    metrics.append(row | result | runtime)
            except (ValueError, OSError, KeyError) as error:
                row.update(status="invalid", error=str(error))
        runs.append(row)
    planned = root / "planned_jobs.json"
    if planned.exists():
        existing = {row["path"] for row in runs}
        for job in json.loads(planned.read_text()).values():
            if job["output"] not in existing:
                runs.append({"run_id": str(Path(job["output"]).relative_to(root.resolve())),
                             "path": job["output"], "status": "not-run", "role": job["config"]["role"],
                             "scenario": job["config"]["scenario"], "replicate": job["config"]["replicate"], "error": ""})
    return pd.DataFrame(runs), pd.DataFrame(metrics)


def bootstrap_summary(frame, draws=2000, seed=20261101):
    rows = []
    if frame.empty:
        return pd.DataFrame()
    for (scenario, method), group in frame.groupby(["scenario", "method"], sort=True):
        if group.replicate.duplicated().any():
            raise ValueError("Duplicate replicate IDs in a scenario/method")
        group = group.sort_values("replicate").reset_index(drop=True)
        estimate = pooled(group)
        alternative = (group.interpretation == "alternative-history call rate").any()
        n = len(group)
        random = np.random.default_rng(seed)
        distributions = {key: [] for key in PERFORMANCE}
        if n >= 5:
            distributions = resampled_values(group, random.integers(0, n, size=(draws, n)))
        for metric in PERFORMANCE:
            meaningful = not alternative or metric == "call_rate"
            values = np.asarray(distributions[metric] if meaningful else [])
            values = values[np.isfinite(values)]
            interval = np.quantile(values, [0.025, 0.975]) if len(values) else [None, None]
            rows.append({"scenario": scenario, "method": method, "metric": metric,
                         "estimate": estimate[metric] if meaningful else None,
                         "ci_low": interval[0], "ci_high": interval[1],
                         "bootstrap_defined_draws": len(values), "bootstrap_draws": draws if n >= 5 else 0,
                         "replicates": n, "qualification": "replicate-bootstrap" if n >= 5 else "descriptive: fewer than five replicates",
                         "interpretation": "alternative-history call rate" if alternative else "simulation accuracy",
                         **{key: estimate[key] for key in ["callable_bases", "called_bases", "true_bases"]}})
    return pd.DataFrame(rows)


def paired_differences(frame, draws=2000, seed=20261101):
    rows = []
    if frame.empty:
        return pd.DataFrame()
    for scenario, group in frame.groupby("scenario", sort=True):
        first = group[group.method == "archie"].set_index("replicate")
        second = group[group.method == "trace"].set_index("replicate")
        ids = sorted(set(first.index) & set(second.index))
        if not ids:
            continue
        first, second = first.loc[ids], second.loc[ids]
        alternative = (group.interpretation == "alternative-history call rate").any()
        keys = ["call_rate"] if alternative else ["recall", "fdr", "precision"]
        estimates = [pooled(first), pooled(second)]
        random = np.random.default_rng(seed)
        samples = {key: [] for key in keys}
        if len(ids) >= 5:
            indices = random.integers(0, len(ids), (draws, len(ids)))
            a, b = resampled_values(first, indices), resampled_values(second, indices)
            samples = {key: (b[key]-a[key])[np.isfinite(b[key]-a[key])] for key in keys}
        for key in keys:
            a, b = estimates[0][key], estimates[1][key]
            interval = np.quantile(samples[key], [0.025, 0.975]) if len(samples[key]) else [None, None]
            rows.append({"scenario": scenario, "metric": key, "difference_trace_minus_archie": b-a if a is not None and b is not None else None,
                         "ci_low": interval[0], "ci_high": interval[1], "paired_replicates": len(ids),
                         "paired_ids": ",".join(map(str, ids)), "bootstrap_defined_draws": len(samples[key]),
                         "qualification": "paired replicate bootstrap" if len(ids) >= 5 else "descriptive: fewer than five paired replicates"})
    return pd.DataFrame(rows)


def error_regions(segments, threshold):
    """Merge contiguous FP/FN pieces without crossing masks, haplotypes or chromosomes."""
    frame = segments.copy()
    called = frame.probability.ge(threshold) if threshold is not None else pd.Series(False, index=frame.index)
    frame["error"] = np.where(called & ~frame.truth.astype(bool), "FP",
                              np.where(~called & frame.truth.astype(bool), "FN", ""))
    rows = []
    for (hap, chrom), group in frame.groupby(["haplotype", "chromosome"], sort=True):
        current = None
        for piece in group.sort_values("start").itertuples():
            if not piece.error:
                if current:
                    rows.append(current)
                    current = None
                continue
            if current and current["end"] == piece.start and current["error"] == piece.error:
                current["end"] = piece.end
                current["bases"] += piece.bases
            else:
                if current:
                    rows.append(current)
                current = {"haplotype": hap, "chromosome": str(chrom), "start": piece.start,
                           "end": piece.end, "bases": piece.bases, "error": piece.error}
        if current:
            rows.append(current)
    return pd.DataFrame(rows, columns=["haplotype", "chromosome", "start", "end", "bases", "error"])


def select_errors(frame):
    candidates = []
    for result in frame.itertuples():
        if result.interpretation != "simulation accuracy" or result.method == "trace_true_arg":
            continue
        data, _ = validate_complete(result.path)
        segments = pd.read_csv(data / result.method / "segments.tsv", sep="\t")
        threshold = None if pd.isna(result.threshold) else float(result.threshold)
        regions = error_regions(segments, threshold)
        for row in regions.to_dict("records"):
            candidates.append(row | {"run_id": result.run_id, "path": result.path, "scenario": result.scenario,
                                      "replicate": result.replicate, "method": result.method})
    if not candidates:
        return []
    table = pd.DataFrame(candidates).sort_values(
        ["bases", "scenario", "replicate", "method", "haplotype", "chromosome", "start"],
        ascending=[False, True, True, True, True, True, True], kind="stable")
    return table.groupby("error", sort=True).head(1).to_dict("records")


def summarize(root, draws=2000):
    root = Path(root)
    destination = root / "summary"
    destination.mkdir(parents=True, exist_ok=True)
    runs, frame = collect(root)
    runs.to_csv(destination / "runs.tsv", sep="\t", index=False)
    frame.to_csv(destination / "replicate_metrics.tsv", sep="\t", index=False)
    # Validation and diagnostic-only runs are never scientific test evidence.
    if not frame.empty:
        frame = frame[frame.eligible_test.eq(True)].copy()
        # Paired comparison uses complete pairs for all group summaries as well.
        paired = frame.groupby(["scenario", "replicate"]).method.transform(lambda x: {"archie", "trace"}.issubset(set(x)))
        frame = frame[paired]
    summary = bootstrap_summary(frame, draws)
    differences = paired_differences(frame, draws)
    summary.to_csv(destination / "performance.tsv", sep="\t", index=False)
    differences.to_csv(destination / "paired_differences.tsv", sep="\t", index=False)
    write_json(destination / "selected_errors.json", select_errors(frame) if not frame.empty else [])
    write_json(destination / "analysis.json", {"version": 2, "bootstrap_unit": "independent simulation replicate",
                                               "bootstrap_seed": 20261101, "bootstrap_draws": draws,
                                               "test_rows": len(frame), "status_counts": runs.status.value_counts().to_dict() if len(runs) else {},
                                               "note": "No ranking is inferred from smoke or preview results"})
    return runs, frame, summary, differences


def preview(segments, threshold):
    return base_metrics(segments, threshold) | {"operating_point": "preview", "threshold_status": "EXPLORATORY PREVIEW - not a frozen research result"}
