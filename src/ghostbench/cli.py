"""Command-line entry points; every long stage can be run independently."""

import argparse
import json
import time
from pathlib import Path

from .io import read_config, write_json


def main():
    parser = argparse.ArgumentParser(prog="ghostbench")
    parser.add_argument("command", choices=["simulate", "train", "run", "smoke", "benchmark",
                                           "evaluate", "report", "empirical-prepare", "empirical-run"])
    parser.add_argument("--config", default="configs/smoke.yaml")
    parser.add_argument("--output", default="results/smoke")
    parser.add_argument("--model", default="work/archie/model.rds")
    args = parser.parse_args()
    config = read_config(args.config)
    output = Path(args.output)
    if args.command == "simulate":
        from .simulation import generate
        generate(config, output)
    elif args.command == "train":
        from .archie import train
        started = time.monotonic()
        train(config, output)
        write_json(output / "total_runtime.json", {"seconds": time.monotonic() - started})
    elif args.command in {"run", "smoke", "benchmark"}:
        from .archie import train
        from .pipeline import benchmark, run
        if args.command == "smoke" and not Path(args.model).exists():
            train(config, Path(args.model).parent)
        if args.command == "benchmark":
            benchmark(config, output, args.model)
        else:
            run(config, output, args.model)
    elif args.command == "evaluate":
        from .pipeline import evaluate
        config = json.loads((output / "config.json").read_text())
        evaluate(output, config)
    elif args.command == "report":
        from .report import create_report
        create_report(output)
    else:
        from .empirical import prepare, run_empirical
        if args.command == "empirical-prepare":
            prepare(config, output)
        else:
            run_empirical(config, output, args.model)


if __name__ == "__main__":
    main()
