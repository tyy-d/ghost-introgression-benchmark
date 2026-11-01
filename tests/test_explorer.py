import io
import json
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from ghostbench.analysis import collect
from ghostbench.evaluation import base_metrics, base_segments
from ghostbench.figures import export_bundle, plotly_tracks, static_tracks, track_data
from ghostbench.io import PREDICTION_COLUMNS, write_json
from ghostbench.registry import digest


@pytest.fixture
def saved_run(tmp_path):
    """Toy saved-result fixture only; never exported as research evidence."""
    run = tmp_path / "fixture"
    data = run / "attempt"
    (data / "inputs").mkdir(parents=True)
    config = {"sequence_length": 100, "role": "diagnostic"}
    write_json(data / "config.json", config)
    mask = pd.DataFrame([{"chromosome": "1", "start": 0, "end": 100}])
    mask.to_csv(data / "inputs/callable.tsv", sep="\t", index=False)
    pd.DataFrame([{"haplotype": "h0", "population": "target"}]).to_csv(data / "inputs/samples.tsv", sep="\t", index=False)
    truth = pd.DataFrame([{"haplotype": "h0", "start": 25, "end": 75}])
    metrics = []
    for method in ["archie", "trace", "trace_true_arg"]:
        (data / method).mkdir()
        calls = pd.DataFrame([["0", method, "h0", "1", 0, 50, 0.8],
                              ["0", method, "h0", "1", 50, 100, 0.2]], columns=PREDICTION_COLUMNS)
        segments = base_segments(calls, truth, ["h0"], mask, require_dense=True)
        segments.to_csv(data / method / "segments.tsv", sep="\t", index=False)
        metrics.append({"method": method, "threshold_status": "diagnostic-unfrozen", "interpretation": "simulation accuracy", **base_metrics(segments)})
    write_json(data / "base_metrics.json", metrics)
    write_json(data / "runtime.json", {"archie_inference_seconds": 1, "trace_end_to_end_seconds": 3})
    files = [path for path in data.rglob("*") if path.is_file()]
    write_json(run / "manifest.json", {"status": "diagnostic-only", "artifact_dir": "attempt", "role": "diagnostic",
                                       "scenario": "fixture", "replicate": 0, "fingerprint": "unit-fixture",
                                       "identity": {"config": config, "thresholds": None},
                                       "artifact_hashes": {str(path.relative_to(data)): digest(path) for path in files}})
    return run


def test_preview_and_download_do_not_modify_saved_results(saved_run):
    before = {str(path): digest(path) for path in saved_run.rglob("*") if path.is_file()}
    frame, errors, mask, view = track_data(saved_run, "h0", 10, 90, 0.9)
    assert not frame.called.any()
    assert view["threshold_status"] == "EXPLORATORY PREVIEW"
    assert (frame.start >= 10).all() and (frame.end <= 90).all()
    figure = static_tracks(frame, errors, mask, view)
    payload = export_bundle(figure, frame, view)
    import matplotlib.pyplot as plt
    plt.close(figure)
    with zipfile.ZipFile(io.BytesIO(payload)) as archive:
        selection = json.loads(archive.read("selection.json"))
        table = pd.read_csv(io.BytesIO(archive.read("data.tsv")), sep="\t")
        assert selection["preview_threshold"] == 0.9
        assert selection["start"] == 10 and selection["end"] == 90
        assert table.bases.sum() == frame.bases.sum()
        assert archive.read("figure.pdf").startswith(b"%PDF")
        assert "figure.svg" in archive.namelist()
    after = {str(path): digest(path) for path in saved_run.rglob("*") if path.is_file()}
    assert before == after


def test_linked_plot_ranges_match_export_selection(saved_run):
    frame, errors, mask, view = track_data(saved_run, "h0", 10, 90)
    figure = plotly_tracks(frame, errors, mask, view)
    for index in range(1, 5):
        name = "xaxis" if index == 1 else f"xaxis{index}"
        assert list(figure.layout[name].range) == [10, 90]
    assert figure.layout.xaxis.matches == "x4"


def test_missing_file_is_a_visible_invalid_run_not_zero_calls(saved_run):
    (saved_run / "attempt/trace/segments.tsv").unlink()
    runs, metrics = collect(saved_run.parent)
    assert runs.status.tolist() == ["invalid"]
    assert metrics.empty


def test_streamlit_views_and_preview_are_operational(saved_run, monkeypatch):
    from streamlit.testing.v1 import AppTest

    monkeypatch.setenv("GHOSTBENCH_RESULTS", str(saved_run.parent))
    app = AppTest.from_file(str(Path("src/ghostbench/explorer.py").resolve()), default_timeout=30).run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == ["Performance", "Ancestry tracks", "Run quality"]
    previews = [box for box in app.checkbox if box.label.startswith("Preview")]
    assert len(previews) == 1
    previews[0].set_value(True).run()
    assert not app.exception
    assert any("EXPLORATORY PREVIEW" in message.value for message in app.warning)
