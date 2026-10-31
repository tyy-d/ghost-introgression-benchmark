import numpy as np
import pandas as pd
import pytest

from ghostbench.evaluation import metrics, overlap_length, score_grid
from ghostbench.io import PREDICTION_COLUMNS, read_config, vcf_position
from ghostbench.vcf import read_vcf


@pytest.fixture
def config():
    return read_config("configs/smoke.yaml")


def test_intervals_do_not_double_count():
    assert overlap_length(0, 20, [(0, 10), (5, 15)]) == 15
    assert overlap_length(10, 20, [(0, 10)]) == 0


def test_callable_weighting_and_ambiguous_windows():
    predictions = pd.DataFrame([
        ["0", "test", "h0", "1", 0, 100, 0.8],
    ], columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame([{"haplotype": "h0", "start": 0, "end": 40}])
    grid = score_grid(predictions, truth, ["h0"], 100, 50, [(0, 25), (50, 100)])
    assert grid.callable_bases.tolist() == [25, 50]
    result = metrics(grid)
    assert result["tp_bases"] == 25
    assert result["fp_bases"] == 50
    assert result["precision"] == pytest.approx(1 / 3)


def test_empty_calls_and_null_truth():
    predictions = pd.DataFrame(columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame(columns=["haplotype", "start", "end"])
    grid = score_grid(predictions, truth, ["h0"], 100, 50, [(0, 100)])
    result = metrics(grid)
    assert result["false_positive_callable_fraction"] == 0
    assert result["recall"] is None
    assert result["precision"] is None


def test_mixed_ancestry_is_separate():
    predictions = pd.DataFrame(columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame([{"haplotype": "h0", "start": 0, "end": 50}])
    grid = score_grid(predictions, truth, ["h0"], 100, 100, [(0, 100)])
    assert metrics(grid)["ambiguous_windows"] == 1
    assert metrics(grid)["evaluated_windows"] == 0


def test_sample_order_and_vcf_roundtrip(config, tmp_path):
    pytest.importorskip("msprime")
    from ghostbench.simulation import export_vcf, sample_metadata, simulate

    ts = simulate(config, 71, length=1000)
    metadata = sample_metadata(ts, config["target_diploids"])
    assert metadata.population.tolist() == ["target"] * 20 + ["reference"] * 20
    export_vcf(ts, tmp_path / "input.vcf")
    positions, matrix, names = read_vcf(tmp_path / "input.vcf")
    assert names == [f"sample{i}" for i in range(20)]
    assert np.array_equal(positions, ts.tables.sites.position)
    assert np.array_equal(matrix, ts.genotype_matrix())
    assert vcf_position(0) == 1


def test_full_and_zero_introgression_labels(config):
    pytest.importorskip("msprime")
    from ghostbench.simulation import donor_intervals, simulate

    full = simulate(config, 83, length=1000, fraction=1.0)
    truth = donor_intervals(full, config["admixture_time"], range(20))
    assert all(overlap_length(0, 1000, value) == 1000 for value in truth.values())
    null = simulate(config, 85, length=1000, fraction=0.0)
    assert all(not value for value in donor_intervals(null, config["admixture_time"], range(20)).values())


def test_overlap_predictions_rejected():
    predictions = pd.DataFrame([
        ["0", "test", "h0", "1", 0, 60, 0.8],
        ["0", "test", "h0", "1", 50, 100, 0.8],
    ], columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame(columns=["haplotype", "start", "end"])
    with pytest.raises(ValueError, match="Overlapping"):
        score_grid(predictions, truth, ["h0"], 100, 50, [(0, 100)])


def test_fully_masked_region_has_no_metric_weight():
    predictions = pd.DataFrame(columns=PREDICTION_COLUMNS)
    truth = pd.DataFrame(columns=["haplotype", "start", "end"])
    grid = score_grid(predictions, truth, ["h0"], 100, 50, [])
    assert metrics(grid)["evaluated_windows"] == 0
    assert metrics(grid)["false_positive_callable_fraction"] is None
