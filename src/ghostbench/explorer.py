"""Read-only local Streamlit explorer. No research jobs are launched here."""

import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

from ghostbench.analysis import bootstrap_summary, collect, paired_differences, preview
from ghostbench.figures import (
    LABELS,
    export_bundle,
    figure_bytes,
    plotly_tracks,
    static_performance,
    static_tracks,
    track_data,
)
from ghostbench.registry import validate_complete

st.set_page_config(page_title="Ghostbench research explorer", layout="wide")


@st.cache_data(show_spinner=False)
def load(root, signature):
    return collect(root)


def signature(root):
    """Include artifact metadata so changed data invalidate the validated cache."""
    records = []
    for path in sorted(Path(root).rglob("manifest.json")):
        records.append((str(path), path.stat().st_mtime_ns))
        try:
            manifest = json.loads(path.read_text())
            for name in manifest.get("artifact_hashes", {}):
                artifact = path.parent / manifest["artifact_dir"] / name
                records.append((str(artifact), artifact.stat().st_mtime_ns, artifact.stat().st_size))
        except (OSError, ValueError, KeyError):
            records.append((str(path), "invalid"))
    return tuple(records)


def main():
    st.title("Ghost introgression: methods and evidence")
    st.caption("Read-only results explorer • Exact callable-base scoring • No jobs or threshold files are modified")
    root = st.sidebar.text_input("Study directory", os.environ.get("GHOSTBENCH_RESULTS", "results/research-v2"))
    if st.sidebar.button("Refresh validated results"):
        load.clear()
    runs, frame = load(root, signature(root))
    if runs.empty:
        st.info("No run manifests found. Run the diagnostic example or convergence pilot from the CLI.")
        return
    performance_tab, tracks_tab, quality_tab = st.tabs(["Performance", "Ancestry tracks", "Run quality"])
    with quality_tab:
        st.subheader("Every scheduled result, including failures")
        st.dataframe(runs, hide_index=True)
        for row in runs.to_dict("records"):
            if row["status"] not in {"complete", "diagnostic-only"}:
                continue
            data, manifest = validate_complete(row["path"])
            with st.expander(f"Configuration and software: {row['run_id']}"):
                st.json(manifest.get("identity", {}))
                quality = data / "singer/quality.json"
                if quality.exists():
                    st.json(json.loads(quality.read_text()))
        for filename in ["pilot/pilot.json", "frozen_thresholds.json", "budget.json", "last_batch.json"]:
            path = Path(root) / filename
            if path.exists():
                with st.expander(filename):
                    st.json(json.loads(path.read_text()))
        st.caption("Pilot-screened settings are not proof of convergence of every test ARG posterior.")
    if frame.empty:
        with performance_tab:
            st.warning("No complete, integrity-checked method outputs. Missing outputs are not no-ancestry calls.")
        return
    available = frame[frame.role != "validation"].copy()
    if available.empty:
        st.info("Only validation results are present. They are not test evidence.")
        return
    diagnostic_default = not available.eligible_test.eq(True).any()
    show_diagnostics = st.sidebar.checkbox("Show diagnostic-only runs", value=diagnostic_default)
    if not show_diagnostics:
        available = available[available.eligible_test.eq(True)]
    scenarios = sorted(available.scenario.unique())
    chosen = st.sidebar.multiselect("Conditions", scenarios, default=scenarios)
    methods = st.sidebar.multiselect("Methods", sorted(available.method.unique()), default=sorted(available.method.unique()), format_func=lambda key: LABELS[key])
    selected = available[available.scenario.isin(chosen) & available.method.isin(methods)]
    with performance_tab:
        if selected.empty:
            st.info("Choose at least one condition and method.")
        else:
            if not selected.eligible_test.eq(True).all():
                st.warning("DIAGNOSTIC / UNFROZEN RESULTS: exploratory functionality evidence, not the calibrated research comparison.")
            summary = bootstrap_summary(selected)
            figure = static_performance(selected, summary)
            st.image(figure_bytes(figure))
            st.dataframe(summary, hide_index=True)
            st.subheader("Paired differences: TRACE minus ArchIE")
            st.dataframe(paired_differences(selected), hide_index=True)
            st.subheader("Per-replicate denominators and cost")
            st.dataframe(selected, hide_index=True)
            selection = {"conditions": chosen, "methods": methods, "units": "callable bases / seconds",
                         "threshold_status": "diagnostic-unfrozen" if not selected.eligible_test.eq(True).all() else "frozen",
                         "run_fingerprints": sorted(selected.fingerprint.unique())}
            st.download_button("Export this performance view (PNG/SVG/PDF + data + metadata)",
                               export_bundle(figure, summary, selection), "performance-view.zip", mime="application/zip")
            import matplotlib.pyplot as plt
            plt.close(figure)
    with tracks_tab:
        options = available.drop_duplicates("path").sort_values(["scenario", "replicate"])
        path = st.selectbox("Run", options.path.tolist(), format_func=lambda value: options.set_index("path").loc[value, "run_id"])
        data, manifest = validate_complete(path)
        config = json.loads((data / "config.json").read_text())
        samples = pd.read_csv(data / "inputs/samples.tsv", sep="\t")
        haplotype = st.selectbox("Target haplotype", samples[samples.population == "target"].haplotype.tolist())
        left, right = st.columns(2)
        start = int(left.number_input("Start (zero-based bp)", 0, config["sequence_length"]-1, 0, key=f"start-{path}"))
        end = int(right.number_input("End (exclusive bp)", start+1, config["sequence_length"], config["sequence_length"], key=f"end-{path}"))
        exploratory = st.checkbox("Preview a different probability threshold (does not change frozen results)")
        threshold = st.slider("Exploratory probability threshold", 0.0, 1.0, 0.5, 0.01) if exploratory else None
        if exploratory:
            st.warning("EXPLORATORY PREVIEW — excluded from frozen research results and report conclusions.")
        frame, errors, mask, view = track_data(path, haplotype, start, end, threshold)
        st.plotly_chart(plotly_tracks(frame, errors, mask, view), width="stretch")
        st.caption("All tracks share the x-axis. White gaps are inaccessible. Black outlines: FP; purple: FN. Export uses the Start/End controls; transient mouse zoom is not saved.")
        st.dataframe(errors, hide_index=True)
        if exploratory:
            st.dataframe(pd.DataFrame([preview(group, threshold) | {"method": method}
                                      for method, group in frame.groupby("method")]))
        figure = static_tracks(frame, errors, mask, view)
        st.download_button("Export this track view (PNG/SVG/PDF + data + metadata)",
                           export_bundle(figure, frame, view), "ancestry-track-view.zip", mime="application/zip")
        import matplotlib.pyplot as plt
        plt.close(figure)


if __name__ == "__main__":
    main()
