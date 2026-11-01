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


SEGMENT_COLUMNS = ["haplotype", "chromosome", "start", "end", "bases", "truth", "probability"]


def base_segments(predictions, truth, haplotypes, callable_frame, require_dense=False):
    """Partition callable sequence at every truth/prediction boundary, without rounding.

    Missing *files* are rejected by the result loader. An explicitly empty call
    table can be scored as no calls; dense detector output must cover every
    callable base when require_dense=True. No mixed-ancestry region is discarded.
    """
    predictions = predictions.copy()
    truth = truth.copy()
    callable_frame = callable_frame.copy()
    for frame in [predictions, truth, callable_frame]:
        if "chromosome" not in frame:
            frame["chromosome"] = "1"
        frame["chromosome"] = frame.chromosome.astype(str)
        if len(frame) and (not np.isfinite(frame[["start", "end"]].to_numpy()).all()
                           or (frame.start < 0).any() or (frame.end <= frame.start).any()):
            raise ValueError("Invalid genomic intervals")
    if len(predictions) and (not predictions.probability.between(0, 1).all()
                             or not set(predictions.haplotype).issubset(set(haplotypes))):
        raise ValueError("Invalid probabilities or unknown haplotypes")
    rows = []
    for chromosome, mask in callable_frame.groupby("chromosome", sort=True):
        masks = merged_intervals(list(zip(mask.start, mask.end)))
        for haplotype in haplotypes:
            calls = predictions[(predictions.haplotype == haplotype)
                                & (predictions.chromosome == chromosome)].sort_values("start")
            if len(calls) > 1 and (calls.start.to_numpy()[1:] < calls.end.to_numpy()[:-1]).any():
                raise ValueError("Overlapping prediction intervals")
            actual = truth[(truth.haplotype == haplotype) & (truth.chromosome == chromosome)]
            positives = merged_intervals(list(zip(actual.start, actual.end)))
            call_rows = list(calls.itertuples())
            for left, right in masks:
                cuts = {left, right}
                for a, b in positives:
                    if a < right and b > left:
                        cuts.update([max(a, left), min(b, right)])
                local_calls = [call for call in call_rows if call.start < right and call.end > left]
                for call in local_calls:
                    cuts.update([max(call.start, left), min(call.end, right)])
                boundaries = sorted(cuts)
                p_index = t_index = 0
                for start, end in zip(boundaries[:-1], boundaries[1:], strict=True):
                    while p_index < len(local_calls) and local_calls[p_index].end <= start:
                        p_index += 1
                    covered = p_index < len(local_calls) and local_calls[p_index].start <= start
                    if require_dense and not covered:
                        raise ValueError(f"Incomplete prediction coverage: {haplotype}:{chromosome}:{start}-{end}")
                    probability = float(local_calls[p_index].probability) if covered else 0.0
                    while t_index < len(positives) and positives[t_index][1] <= start:
                        t_index += 1
                    is_true = t_index < len(positives) and positives[t_index][0] <= start
                    rows.append([haplotype, chromosome, start, end, end-start, is_true, probability])
    return pd.DataFrame(rows, columns=SEGMENT_COLUMNS)


def base_metrics(segments, threshold=0.5):
    weights = segments.bases.to_numpy(dtype=float)
    actual = segments.truth.to_numpy(dtype=bool)
    called = (segments.probability.to_numpy(dtype=float) >= threshold
              if threshold is not None else np.zeros(len(segments), dtype=bool))
    tp = float(weights[called & actual].sum())
    fp = float(weights[called & ~actual].sum())
    fn = float(weights[~called & actual].sum())
    tn = float(weights[~called & ~actual].sum())
    total = tp + fp + fn + tn

    def ratio(a, b):
        return a / b if b else None

    return {"threshold": threshold, "tp_bases": tp, "fp_bases": fp, "fn_bases": fn,
            "tn_bases": tn, "callable_bases": total, "true_bases": tp + fn,
            "called_bases": tp + fp, "precision": ratio(tp, tp+fp),
            "recall": ratio(tp, tp+fn), "fdr": ratio(fp, tp+fp),
            "f1": ratio(2*tp, 2*tp+fp+fn),
            "false_called_fraction": ratio(fp, total), "false_called_kb_per_mb": ratio(fp*1000, total),
            "call_rate": ratio(tp+fp, total), "operating_point": "abstain" if threshold is None else "frozen"}


def base_curve(segments):
    """Exact weighted PR curve in O(n log n), pooling all validation replicates."""
    if segments.empty:
        return pd.DataFrame(columns=["threshold", "tp_bases", "fp_bases", "called_bases", "precision", "recall", "fdr"])
    frame = segments.assign(
        tp=segments.bases * segments.truth.astype(int),
        fp=segments.bases * (~segments.truth.astype(bool)).astype(int),
    ).groupby("probability", sort=True)[["tp", "fp"]].sum().sort_index(ascending=False).cumsum()
    frame = frame.rename(columns={"tp": "tp_bases", "fp": "fp_bases"})
    frame["called_bases"] = frame.tp_bases + frame.fp_bases
    frame["precision"] = frame.tp_bases / frame.called_bases
    frame["fdr"] = frame.fp_bases / frame.called_bases
    total_true = float(segments.loc[segments.truth.astype(bool), "bases"].sum())
    frame["recall"] = frame.tp_bases / total_true if total_true else np.nan
    return frame.reset_index(names="threshold")


def fdr_threshold(segments, target=0.05):
    if not 0 < target < 1:
        raise ValueError("FDR target must be between zero and one")
    curve = base_curve(segments)
    eligible = curve[(curve.called_bases > 0) & (curve.fdr <= target)].dropna(subset=["recall"])
    if eligible.empty:
        return {"threshold": None, "status": "no qualifying operating point", "target_fdr": target}
    best = eligible.sort_values(["recall", "threshold"], ascending=False).iloc[0]
    return {"threshold": float(best.threshold), "status": "qualified", "target_fdr": target,
            "validation_recall": float(best.recall), "validation_fdr": float(best.fdr),
            "validation_called_bases": float(best.called_bases)}
