from pathlib import Path

import numpy as np
import pytest

from ghostbench.archie import features, s_star, upstream
from ghostbench.empirical import mask_intervals, select_samples
from ghostbench.vcf import read_vcf


def test_s_star_matches_upstream_when_positions_are_aligned():
    if not Path("vendor/ArchIE/data/calc_stats_window_data.py").exists():
        pytest.skip("Upstream ArchIE checkout is required")
    haplotypes = np.array([[1, 1, 1], [1, 1, 1], [0, 0, 0], [0, 0, 0]])
    positions = np.array([10, 2000, 4000])
    assert s_star(haplotypes, positions, 0) == upstream().s_star_ind(haplotypes, positions, 0)


def test_s_star_preserves_positions_when_focal_has_missing_derived_sites():
    haplotypes = np.array([[1, 0, 1], [1, 0, 1], [0, 1, 0], [0, 1, 0]])
    assert s_star(haplotypes, np.array([100, 200, 1100]), 0) == 6000


def test_no_snps_produces_finite_sample_dependent_features():
    if not Path("vendor/ArchIE/data/calc_stats_window_data.py").exists():
        pytest.skip("Upstream ArchIE checkout is required")
    matrix = features(np.zeros((0, 8), dtype=int), np.array([]), 4, 0, 50000)
    assert matrix.shape == (4, 16)
    assert np.isfinite(matrix).all()


def test_dense_mask_conversion():
    assert mask_intervals([False, True, True, False, True]) == [(1, 3), (4, 5)]
    assert mask_intervals([]) == []


def test_empirical_selection_needs_quality_and_relatedness():
    import pandas as pd
    with pytest.raises(ValueError, match="QC/relatedness"):
        select_samples(pd.DataFrame({"sample_id": ["a"]}), {})


def test_vcf_rejects_unphased_input(tmp_path):
    path = tmp_path / "bad.vcf"
    path.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
                    "1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
    with pytest.raises(ValueError, match="phased"):
        read_vcf(path)


def test_model_with_different_sample_count_is_rejected_before_inference(tmp_path):
    import json

    import pandas as pd

    from ghostbench.archie import predict

    path = tmp_path / "input.vcf"
    path.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n"
                    "1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t0|0\n")
    metadata = pd.DataFrame({"individual": ["a", "a", "b", "b"],
                             "population": ["target", "target", "reference", "reference"]})
    (tmp_path / "training_summary.json").write_text(json.dumps({
        "target_haplotypes": 20, "reference_haplotypes": 20, "window_size": 50000,
    }))
    with pytest.raises(ValueError, match="retrain"):
        predict(path, metadata, {"window_size": 50000}, tmp_path / "model.rds", tmp_path / "out")
