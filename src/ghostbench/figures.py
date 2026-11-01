"""Shared view data, interactive figures, and headless publication exports."""

import io
import json
import shutil
import zipfile
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

from .analysis import bootstrap_summary, error_regions
from .registry import validate_complete

matplotlib.use("Agg")

COLORS = {"archie": "#0072B2", "trace": "#D55E00", "trace_true_arg": "#777777", "truth": "#009E73"}
LABELS = {"archie": "ArchIE", "trace": "SINGER + TRACE", "trace_true_arg": "TRACE: true ARG (diagnostic)"}


def track_data(root, haplotype, start=0, end=None, preview_threshold=None):
    data, manifest = validate_complete(root)
    config = json.loads((data / "config.json").read_text())
    end = min(config["sequence_length"], end if end is not None else config["sequence_length"])
    if not 0 <= start < end:
        raise ValueError("Track range must lie inside the sequence")
    metrics = {row["method"]: row for row in json.loads((data / "base_metrics.json").read_text())}
    frames, regions = [], []
    for method in ["archie", "trace", "trace_true_arg"]:
        segments = pd.read_csv(data / method / "segments.tsv", sep="\t")
        segments = segments[(segments.haplotype == haplotype) & (segments.end > start) & (segments.start < end)].copy()
        segments["start"] = segments.start.clip(lower=start)
        segments["end"] = segments.end.clip(upper=end)
        segments["bases"] = segments.end-segments.start
        threshold = preview_threshold if preview_threshold is not None else metrics[method]["threshold"]
        segments["called"] = segments.probability.ge(threshold) if threshold is not None else False
        segments["method"] = method
        frames.append(segments)
        errors = error_regions(segments, threshold)
        errors["method"] = method
        regions.append(errors)
    frame = pd.concat(frames, ignore_index=True)
    if frame.empty:
        raise ValueError("Selected range/haplotype has no callable bases")
    mask = pd.read_csv(data / "inputs/callable.tsv", sep="\t")
    view = {"run_id": str(root), "fingerprint": manifest["fingerprint"], "haplotype": haplotype,
            "chromosome": "1", "start": start, "end": end, "units": "zero-based half-open bp",
            "threshold_status": "EXPLORATORY PREVIEW" if preview_threshold is not None else metrics["archie"]["threshold_status"],
            "preview_threshold": preview_threshold, "stored_thresholds": {key: row["threshold"] for key, row in metrics.items()},
            "role": manifest["role"], "scenario": manifest["scenario"]}
    return frame, pd.concat(regions, ignore_index=True), mask, view


def plotly_tracks(frame, errors, mask, view):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    methods = list(frame.method.unique())
    figure = make_subplots(rows=4, cols=1, shared_xaxes=True, vertical_spacing=0.07,
                          subplot_titles=["True donor ancestry / accessible sequence", *[LABELS[x] for x in methods]])
    truth = frame[frame.method == methods[0]]
    for row_index, (name, pieces) in enumerate([("truth", truth), *[(m, frame[frame.method == m]) for m in methods]], 1):
        for state in [False, True]:
            selected = pieces[pieces.truth == state] if name == "truth" else pieces[pieces.called == state]
            figure.add_trace(go.Bar(x=selected.end-selected.start, base=selected.start,
                y=[0.5]*len(selected), orientation="h", width=0.7,
                marker_color=COLORS.get(name) if state else "#dddddd", opacity=1 if state else 0.5,
                name=f"{name}: {'donor' if state else 'other'}", showlegend=False,
                customdata=selected[["start", "end", "probability"]].to_numpy(),
                hovertemplate="[%{customdata[0]}, %{customdata[1]}) bp<br>Probability: %{customdata[2]:.4f}<extra></extra>"),
                row=row_index, col=1)
        figure.update_yaxes(range=[0, 1], showticklabels=False, row=row_index, col=1)
    for error in errors[errors.method != "trace_true_arg"].itertuples():
        row = methods.index(error.method)+2
        figure.add_vrect(x0=error.start, x1=error.end, line_width=1,
                         line_color="#CC79A7" if error.error == "FN" else "#000000", fillcolor="rgba(0,0,0,0)", row=row, col=1)
    figure.update_xaxes(range=[view["start"], view["end"]])
    figure.update_xaxes(title_text="Position (bp; zero-based)", row=4, col=1)
    figure.update_layout(height=540, barmode="overlay", template="plotly_white", title=f"{view['haplotype']} — {view['threshold_status']}",
                         margin=dict(l=30, r=20, t=80, b=30))
    return figure


def static_tracks(frame, errors, mask, view):
    import matplotlib.pyplot as plt

    methods = list(frame.method.unique())
    fig, axes = plt.subplots(4, 1, figsize=(10, 5.5), sharex=True, constrained_layout=True)
    truth = frame[frame.method == methods[0]]
    for axis, name, pieces in zip(axes, ["truth", *methods], [truth, *[frame[frame.method == m] for m in methods]], strict=True):
        chosen = pieces[pieces.truth.astype(bool)] if name == "truth" else pieces[pieces.called]
        axis.broken_barh([(max(view["start"], a), min(view["end"], b)-max(view["start"], a))
                          for a, b in zip(mask.start, mask.end) if b > view["start"] and a < view["end"]], (0.2, 0.6), facecolors="#eeeeee")
        axis.broken_barh(list(zip(chosen.start, chosen.end-chosen.start)), (0.2, 0.6), facecolors=COLORS[name])
        for error in errors[errors.method == name].itertuples():
            axis.axvspan(error.start, error.end, fill=False, edgecolor="#CC79A7" if error.error == "FN" else "black", linewidth=1)
        axis.set(ylim=(0, 1), yticks=[], ylabel="")
        axis.text(0.005, 0.96, "Truth" if name == "truth" else LABELS[name],
                  transform=axis.transAxes, ha="left", va="top", fontsize=8,
                  bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1})
        axis.set_xlim(view["start"], view["end"])
    axes[0].set_title(f"{view['haplotype']} | {view['scenario']} | {view['threshold_status']}")
    axes[-1].set_xlabel("Position (bp; zero-based half-open intervals)")
    return fig


def static_performance(frame, summary=None):
    import matplotlib.pyplot as plt

    summary = bootstrap_summary(frame) if summary is None else summary
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    if summary.empty:
        for axis in axes:
            axis.text(0.5, 0.5, "No eligible test results", ha="center", transform=axis.transAxes)
        return fig
    scenarios = list(summary.scenario.unique())
    methods = list(summary.method.unique())
    for axis, metric in zip(axes, ["recall", "fdr"], strict=True):
        for index, method in enumerate(methods):
            points = summary[(summary.method == method) & (summary.metric == metric)]
            x = np.array([scenarios.index(s) for s in points.scenario]) + (index-1)*0.16
            axis.scatter(x, points.estimate, label=LABELS[method], color=COLORS[method], s=22)
            for position, point in zip(x, points.itertuples(), strict=True):
                if pd.notna(point.ci_low):
                    axis.vlines(position, point.ci_low, point.ci_high, color=COLORS[method])
        axis.set_xticks(range(len(scenarios)), scenarios, rotation=25, ha="right")
        axis.set_ylim(-0.02, 1.03)
        axis.set_ylabel("Recall (all callable bases)" if metric == "recall" else "Achieved test FDR")
    axes[1].axhline(0.05, color="black", linestyle="--", linewidth=0.8, label="Validation target: 5%")
    axes[0].legend(fontsize=7)
    return fig


def static_cost(frame):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    if frame.empty:
        for axis in axes:
            axis.text(0.5, 0.5, "No eligible test results", ha="center", transform=axis.transAxes)
        return fig
    unique = frame.drop_duplicates(["scenario", "replicate"])
    a = unique.archie_inference_seconds
    b = unique.trace_end_to_end_seconds
    axes[0].boxplot([a, b], tick_labels=["ArchIE", "SINGER + TRACE"])
    axes[0].set(ylabel="Inference wall time (s)", yscale="log")
    stress = frame[frame.scenario.isin(["baseline", "reference_admixed", "demography_mismatch"])]
    stress_summary = bootstrap_summary(stress)
    if len(stress_summary):
        for method, group in stress_summary[stress_summary.metric == "recall"].groupby("method"):
            axes[1].plot(group.scenario, group.estimate, "o-", color=COLORS[method], label=LABELS[method])
    axes[1].set_ylabel("Recall: baseline and robustness conditions")
    axes[1].tick_params(axis="x", rotation=20)
    if len(stress):
        axes[1].legend(fontsize=7)
    return fig


def figure_bytes(figure, format="png"):
    output = io.BytesIO()
    figure.savefig(output, format=format, dpi=180, bbox_inches="tight")
    return output.getvalue()


def export_bundle(figure, table, selection):
    """Identical data/selection contract for the explorer download and CLI export."""
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for format in ["png", "svg", "pdf"]:
            archive.writestr(f"figure.{format}", figure_bytes(figure, format))
        archive.writestr("data.tsv", table.to_csv(sep="\t", index=False))
        archive.writestr("selection.json", json.dumps(selection, indent=2, allow_nan=False))
    return output.getvalue()


def export_study(root, destination):
    import matplotlib.pyplot as plt

    from .analysis import select_errors, summarize

    root, destination = Path(root), Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    runs, frame, summary, differences = summarize(root)
    figures = {"performance": static_performance(frame, summary), "cost_robustness": static_cost(frame)}
    selected = select_errors(frame) if not frame.empty else []
    if selected:
        selected.sort(key=lambda row: (-row["bases"], row["scenario"], row["replicate"], row["method"], row["start"]))
        chosen = selected[0]
        data, manifest = validate_complete(chosen["path"])
        config = json.loads((data / "config.json").read_text())
        margin = 50000
        args = track_data(chosen["path"], chosen["haplotype"], max(0, chosen["start"]-margin), min(config["sequence_length"], chosen["end"]+margin))
        figures["localization"] = static_tracks(*args)
        args[0].to_csv(destination / "localization_data.tsv", sep="\t", index=False)
        (destination / "localization_selection.json").write_text(json.dumps(args[3], indent=2))
    else:
        fig, ax = plt.subplots(figsize=(9, 3))
        ax.text(0.5, 0.5, "No eligible FP/FN case available; no illustrative case selected", ha="center", transform=ax.transAxes)
        ax.axis("off")
        figures["localization"] = fig
    for name, figure in figures.items():
        for format in ["png", "svg", "pdf"]:
            (destination / f"{name}.{format}").write_bytes(figure_bytes(figure, format))
        plt.close(figure)
    for name, table in [("performance", summary), ("paired_differences", differences), ("replicate_metrics", frame), ("run_status", runs)]:
        table.to_csv(destination / f"{name}.tsv", sep="\t", index=False)
    selection = {"threshold_status": "frozen research results only", "preview": False,
                 "root": str(root.resolve()), "replicate_bootstrap_draws": 2000,
                 "run_fingerprints": sorted(frame.fingerprint.unique().tolist()) if len(frame) else []}
    (destination / "selection.json").write_text(json.dumps(selection, indent=2))
    return frame, summary, differences, runs


def pack_view(root, destination):
    """Portable processed-result bundle: no raw genomes, ARGs, models, or resume."""
    import hashlib

    from .io import write_json
    from .registry import digest

    data, manifest = validate_complete(root)
    destination = Path(destination)
    existing = destination / "manifest.json"
    if existing.exists() and json.loads(existing.read_text()).get("bundle_kind") != "processed-view-only":
        raise ValueError("Refusing to overwrite a real inference run with a viewer bundle")
    target = destination / "processed"
    names = ["config.json", "runtime.json", "base_metrics.json", "provenance.json", "inputs/samples.tsv", "inputs/callable.tsv"]
    names += [f"{method}/segments.tsv" for method in ["archie", "trace", "trace_true_arg"]]
    if (data / "singer/quality.json").exists():
        names.append("singer/quality.json")
    for name in names:
        (target / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(data / name, target / name)
    packed = {key: manifest[key] for key in ["role", "scenario", "replicate", "identity", "status"]}
    packed.update(schema=2, bundle_kind="processed-view-only", artifact_dir="processed",
                  source_fingerprint=manifest["fingerprint"],
                  fingerprint=hashlib.sha256(("view-only:"+manifest["fingerprint"]).encode()).hexdigest(),
                  artifact_hashes={name: digest(target / name) for name in names},
                  note="Processed probabilities and ancestry only. Raw inputs and ARGs are not distributed. Fingerprint differs from an inference run, so this bundle cannot be resumed as one.")
    write_json(existing, packed)
    return destination
