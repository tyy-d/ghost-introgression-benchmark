import json
import sys
import time

import numpy as np
import pandas as pd
import psutil
import pytest

from ghostbench.analysis import bootstrap_summary, error_regions, paired_differences
from ghostbench.convergence import chain_diagnostics, diagnose_tables
from ghostbench.evaluation import base_curve, base_metrics, base_segments, fdr_threshold
from ghostbench.io import PREDICTION_COLUMNS, run_command
from ghostbench.registry import digest, validate_complete
from ghostbench.study import Budget, experiment_config


def example():
    predictions = pd.DataFrame([["0", "test", "h0", "1", 0, 100, 0.8]], columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame([{"haplotype": "h0", "start": 0, "end": 40}])
    mask = pd.DataFrame([{"chromosome": "1", "start": 0, "end": 25},
                         {"chromosome": "1", "start": 50, "end": 100}])
    return predictions, truth, mask


def test_exact_mixed_bases_are_not_discarded():
    calls, truth, _ = example()
    mask = pd.DataFrame([{"chromosome": "1", "start": 0, "end": 100}])
    segments = base_segments(calls, truth, ["h0"], mask, require_dense=True)
    result = base_metrics(segments)
    assert result["tp_bases"] == 40
    assert result["fp_bases"] == 60
    assert result["fdr"] == 0.6
    assert result["recall"] == 1
    assert result["callable_bases"] == 100


def test_masks_and_adjacent_boundaries_count_exactly_once():
    calls, truth, mask = example()
    segments = base_segments(calls, truth, ["h0"], mask)
    result = base_metrics(segments)
    assert result["tp_bases"] == 25
    assert result["fp_bases"] == 50
    assert result["callable_bases"] == 75
    assert result["false_called_kb_per_mb"] == pytest.approx(1000*50/75)


def test_empty_calls_are_distinct_from_partial_dense_output():
    calls, truth, mask = example()
    empty = pd.DataFrame(columns=PREDICTION_COLUMNS)
    result = base_metrics(base_segments(empty, truth, ["h0"], mask))
    assert result["precision"] is None
    assert result["fdr"] is None
    assert result["recall"] == 0
    with pytest.raises(ValueError, match="Incomplete"):
        base_segments(empty, truth, ["h0"], mask, require_dense=True)
    calls.loc[0, "end"] = 20
    with pytest.raises(ValueError, match="Incomplete"):
        base_segments(calls, truth, ["h0"], mask, require_dense=True)


def test_fdr_calibration_abstains_instead_of_claiming_empty_success():
    calls, truth, mask = example()
    segments = base_segments(calls, truth, ["h0"], mask)
    point = fdr_threshold(segments)
    assert point["status"] == "no qualifying operating point"
    assert point["threshold"] is None
    assert base_metrics(segments, None)["fdr"] is None


def test_fdr_threshold_ties_choose_stricter_nonempty_point():
    segments = pd.DataFrame({"probability": [0.9, 0.8, 0.7], "truth": [True, False, False], "bases": [100, 1, 100]})
    result = fdr_threshold(segments)
    assert result["threshold"] == 0.9
    assert result["validation_recall"] == 1
    assert base_curve(segments).iloc[-1].fp_bases == 101


def metric_rows(n=5):
    rows = []
    for replicate in range(n):
        for method in ["archie", "trace"]:
            rows.append({"scenario": "baseline", "replicate": replicate, "method": method,
                         "interpretation": "simulation accuracy", "tp_bases": 10+replicate,
                         "fp_bases": replicate, "fn_bases": 10, "tn_bases": 100})
    return pd.DataFrame(rows)


def test_bootstrap_is_paired_at_independent_replicate_level():
    differences = paired_differences(metric_rows(), draws=100)
    assert (differences.difference_trace_minus_archie == 0).all()
    assert (differences.ci_low == 0).all() and (differences.ci_high == 0).all()
    assert (differences.paired_replicates == 5).all()
    assert (differences.paired_ids == "0,1,2,3,4").all()
    short = bootstrap_summary(metric_rows(4), draws=100)
    assert short.ci_low.isna().all()
    assert (short.bootstrap_draws == 0).all()


def test_structure_is_not_mislabeled_false_positive_evidence():
    frame = metric_rows()
    frame["interpretation"] = "alternative-history call rate"
    result = bootstrap_summary(frame, draws=100)
    assert result[result.metric != "call_rate"].estimate.isna().all()
    assert result[result.metric == "call_rate"].estimate.notna().all()


def test_validation_test_and_training_seeds_are_disjoint(config):
    validation = {200000+2*i for i in range(6)}
    testing = {300000+1000*i+2*j for i in range(6) for j in range(10)}
    training = {100000+2*i for i in range(2000)}
    assert not validation & testing and not validation & training and not testing & training
    first = experiment_config(config | {"workers": 1}, "test", 300000)
    second = experiment_config(config | {"workers": 4}, "test", 300000)
    assert first == second


def test_completed_output_hashes_reject_missing_or_modified_results(tmp_path):
    data = tmp_path / "attempt"
    data.mkdir()
    file = data / "predictions.tsv"
    file.write_text("real saved prediction bytes")
    manifest = {"status": "complete", "artifact_dir": "attempt", "artifact_hashes": {"predictions.tsv": digest(file)}}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    assert validate_complete(tmp_path)[0] == data
    file.write_text("changed prediction bytes")
    with pytest.raises(ValueError, match="altered"):
        validate_complete(tmp_path)
    file.unlink()
    with pytest.raises(ValueError, match="Missing"):
        validate_complete(tmp_path)


def test_incomplete_status_is_never_eligible(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"status": "interrupted"}))
    with pytest.raises(ValueError, match="not complete"):
        validate_complete(tmp_path)


def test_budget_accounts_elapsed_batches_and_rejects_unreconciled_state(tmp_path):
    with Budget(tmp_path, 10) as budget:
        assert 0 < budget.remaining() <= 10
    assert json.loads((tmp_path / "budget.json").read_text())["running"] is False
    (tmp_path / "budget.json").write_text(json.dumps({"running": True, "elapsed_seconds": 1}))
    with pytest.raises(ValueError, match="reconciled"):
        with Budget(tmp_path, 10):
            pass


def test_convergence_distinguishes_constants_and_drifting_chains():
    rng = np.random.default_rng(42)
    good = chain_diagnostics(rng.normal(size=(2, 400)))
    assert good["status"] == "pass"
    bad = chain_diagnostics(np.array([np.arange(400), np.arange(400)+400]))
    assert bad["status"] == "nonmixing"
    assert chain_diagnostics(np.zeros((2, 100)))["status"] == "constant-uninformative"
    tables = [pd.DataFrame({name: np.zeros(100) for name in ["diversity_fit", "incompatibilities", "ancestry_fraction_h0"]}) for _ in range(2)]
    assert diagnose_tables(tables)["status"] == "unresolved"


def test_error_selection_does_not_merge_across_mask_gaps():
    calls, truth, mask = example()
    regions = error_regions(base_segments(calls, truth, ["h0"], mask), 0.5)
    assert regions.error.tolist() == ["FP"]
    assert regions.start.tolist() == [50]
    assert regions.end.tolist() == [100]


def test_timeout_kills_stubborn_grandchild_and_retains_status(tmp_path):
    pidfile = tmp_path / "child.pid"
    child = "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)"
    parent = ("import subprocess,sys,time; from pathlib import Path; "
              f"p=subprocess.Popen([sys.executable,'-c',{child!r}]); "
              f"Path({str(pidfile)!r}).write_text(str(p.pid)); time.sleep(60)")
    with pytest.raises(TimeoutError):
        run_command([sys.executable, "-c", parent], tmp_path / "command.log", timeout=0.6)
    pid = int(pidfile.read_text())
    time.sleep(0.1)
    try:
        assert psutil.Process(pid).status() == psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        pass
    assert json.loads((tmp_path / "command.command.json").read_text())["status"] == "timeout"


@pytest.fixture
def config():
    from ghostbench.io import read_config
    return read_config("configs/smoke.yaml")


def test_cleanup_denial_does_not_hide_timeout(tmp_path, monkeypatch):
    import ghostbench.io as io

    original = io.terminate_group

    def denied_after_cleanup(process):
        original(process)
        raise PermissionError("sandbox denied redundant group signal")

    monkeypatch.setattr(io, "terminate_group", denied_after_cleanup)
    with pytest.raises(TimeoutError):
        run_command([sys.executable, "-c", "import time; time.sleep(60)"],
                    tmp_path / "denied.log", timeout=0.2)
    telemetry = json.loads((tmp_path / "denied.command.json").read_text())
    assert telemetry["status"] == "timeout"
    assert "PermissionError" in telemetry["cleanup_error"]


def test_run_rejects_changed_identity_and_retries_in_fresh_attempt(tmp_path, monkeypatch):
    import ghostbench.registry as registry

    monkeypatch.setattr(registry, "identity", lambda config, model, thresholds:
                        (config["key"], {"config": config}))
    with pytest.raises(KeyboardInterrupt):
        with registry.Run(tmp_path, {"key": "first"}, None, None):
            raise KeyboardInterrupt()
    assert json.loads((tmp_path / "manifest.json").read_text())["status"] == "interrupted"
    with pytest.raises(ValueError, match="Stale run"):
        with registry.Run(tmp_path, {"key": "changed"}, None, None):
            pass
    with pytest.raises(KeyboardInterrupt):
        with registry.Run(tmp_path, {"key": "first"}, None, None) as run:
            assert run.data.name == "0002"
            assert (tmp_path / "attempts/0001").is_dir()
            raise KeyboardInterrupt()
