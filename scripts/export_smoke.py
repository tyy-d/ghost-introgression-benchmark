"""Publish only compact, completed synthetic artifacts; keep raw data ignored."""

import argparse
import json
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="results/smoke-verified")
    parser.add_argument("--output", default="reports/smoke")
    parser.add_argument("--training", default="work/archie/training_summary.json")
    args = parser.parse_args()
    source, destination = Path(args.input), Path(args.output)
    methods = {row["method"] for row in json.loads((source / "metrics.json").read_text())}
    if methods != {"archie", "trace", "trace_true_arg"}:
        raise ValueError("Export requires both actual methods and the true-ARG diagnostic")
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ["metrics.json", "runtime.json", "config.json", "provenance.json",
                     "report.md", "benchmark.png"]:
        shutil.copy2(source / filename, destination / filename)
    shutil.copy2(args.training, destination / "training_summary.json")
    for method in sorted(methods):
        (destination / method).mkdir(exist_ok=True)
        for filename in ["grid.tsv", "precision_recall.tsv"]:
            shutil.copy2(source / method / filename, destination / method / filename)


if __name__ == "__main__":
    main()
