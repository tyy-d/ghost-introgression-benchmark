"""Callable-base-weighted evaluation on a common genomic grid."""

import numpy as np
import pandas as pd


def merged_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def overlap_length(start, end, intervals):
    return sum(max(0, min(end, right) - max(start, left))
               for left, right in merged_intervals(intervals))


def intersect_intervals(first, second):
    return merged_intervals([(max(a, c), min(b, d)) for a, b in first for c, d in second
                             if min(b, d) > max(a, c)])


def score_grid(predictions, truth, haplotypes, length, window_size, callable_intervals):
    """Missing predictions mean zero probability; overlapping calls are rejected.

    Truth fractions and scores are averaged only over callable bases. A region
    with no callable bases contributes no observations or metric weight.
    """
    if length <= 0 or window_size <= 0:
        raise ValueError("Sequence length and window size must be positive")
    rows = []
    for haplotype in haplotypes:
        calls = predictions[predictions.haplotype == haplotype].sort_values("start")
        if len(calls) > 1 and (calls.start.to_numpy()[1:] < calls.end.to_numpy()[:-1]).any():
            raise ValueError("Overlapping prediction intervals would double-count bases")
        actual = truth[truth.haplotype == haplotype]
        true_callable = intersect_intervals(
            list(zip(actual.start, actual.end)), callable_intervals
        )
        for start in range(0, length, window_size):
            end = min(length, start + window_size)
            weight = overlap_length(start, end, callable_intervals)
            if not weight:
                continue
            probability_sum = sum(
                float(call.probability) * overlap_length(
                    max(start, int(call.start)), min(end, int(call.end)), callable_intervals
                ) for call in calls.itertuples() if call.end > start and call.start < end
            )
            fraction = overlap_length(start, end, true_callable) / weight
            label = 1 if fraction > 0.7 else 0 if fraction < 0.3 else -1
            rows.append({"haplotype": haplotype, "start": start, "end": end,
                         "callable_bases": weight, "truth_fraction": fraction,
                         "label": label, "probability": probability_sum / weight})
    return pd.DataFrame(rows, columns=["haplotype", "start", "end", "callable_bases",
                                       "truth_fraction", "label", "probability"])


def metrics(grid, threshold=0.5):
    eligible = grid[grid.label != -1]
    positive = eligible.probability.to_numpy() >= threshold
    actual = eligible.label.to_numpy() == 1
    weights = eligible.callable_bases.to_numpy()
    tp = float(weights[positive & actual].sum())
    fp = float(weights[positive & ~actual].sum())
    fn = float(weights[~positive & actual].sum())
    tn = float(weights[~positive & ~actual].sum())

    def divide(numerator, denominator):
        return numerator / denominator if denominator else None

    return {
        "threshold": float(threshold), "tp_bases": tp, "fp_bases": fp,
        "fn_bases": fn, "tn_bases": tn,
        "precision": divide(tp, tp + fp), "recall": divide(tp, tp + fn),
        "f1": divide(2 * tp, 2 * tp + fp + fn),
        "false_positive_callable_fraction": divide(fp, fp + tn),
        "ambiguous_windows": int((grid.label == -1).sum()),
        "ambiguous_callable_bases": int(grid.loc[grid.label == -1, "callable_bases"].sum()),
        "evaluated_windows": int(len(eligible)),
    }


def precision_recall(grid):
    thresholds = np.unique(np.concatenate(([0.0], grid.probability.to_numpy(), [1.0])))
    return pd.DataFrame([metrics(grid, threshold) for threshold in thresholds])


def select_threshold(grid):
    """Freeze a validation-only F1 threshold; prefer stricter thresholds on ties."""
    curve = precision_recall(grid)
    valid = curve.dropna(subset=["f1"])
    if valid.empty:
        return 0.5
    return float(valid.sort_values(["f1", "threshold"], ascending=False).iloc[0].threshold)
