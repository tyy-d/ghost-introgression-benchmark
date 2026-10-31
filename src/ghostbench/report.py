"""Export compact figures and an honest research-status report."""

import json
import os
from pathlib import Path


def create_report(output):
    output = Path(output)
    os.environ.setdefault("MPLCONFIGDIR", str(output / ".mplcache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd

    metrics = json.loads((output / "metrics.json").read_text())
    runtime = json.loads((output / "runtime.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for row in metrics:
        curve = pd.read_csv(output / row["method"] / "precision_recall.tsv", sep="\t")
        axes[0].plot(curve.recall, curve.precision, label=row["method"])
    axes[0].set(xlabel="Recall", ylabel="Precision", xlim=(0, 1), ylim=(0, 1.02))
    axes[0].legend()
    keys = ["archie_inference_seconds", "trace_end_to_end_seconds"]
    axes[1].bar(["ArchIE", "SINGER + TRACE"], [runtime[key] for key in keys])
    axes[1].set(ylabel="Wall time (seconds)")
    fig.savefig(output / "benchmark.png", dpi=160)
    plt.close(fig)
    lines = ["# Synthetic functionality run", "",
             "This single run checks functionality; it does not establish comparative accuracy.", "",
             "Parameters are software-validation settings, not mosquito estimates. "
             "All reported metrics use callable-base weighting on a common 50-kb grid. "
             "Windows with 30–70% true donor ancestry are reported separately.", "",
             "| Method | Precision | Recall | F1 | Ambiguous windows |",
             "|---|---:|---:|---:|---:|"]
    for row in metrics:
        def display(value):
            return "undefined" if value is None else f"{value:.3f}"
        lines.append(f"| {row['method']} | {display(row['precision'])} | "
                     f"{display(row['recall'])} | {display(row['f1'])} | {row['ambiguous_windows']} |")
    lines.extend(["", "![Precision–recall and runtime](benchmark.png)", "",
                  "## Limitations", "",
                  "One seed is not a performance estimate. The short SINGER chain is a smoke test, "
                  "not a convergence demonstration. ArchIE uses a comparison population and "
                  "simulation-trained parameters; TRACE uses inferred genealogies. Their differing "
                  "information requirements must be considered when interpreting performance. "
                  "The true-ARG TRACE result is diagnostic and is not an empirical inference result.",
                  "", "## Empirical status", "",
                  "The mosquito pilot has not been run by this synthetic command. "
                  "Empirical overlap is not an accuracy measurement.", ""])
    (output / "report.md").write_text("\n".join(lines))
