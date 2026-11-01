"""Command-line entry points; every long stage can be run independently."""

import argparse
import json
import time
from pathlib import Path

from .io import read_config, write_json


def main():
    parser = argparse.ArgumentParser(prog="ghostbench")
    parser.add_argument("command", choices=["simulate", "train", "run", "smoke", "benchmark",
                                           "evaluate", "report", "empirical-prepare", "empirical-run",
                                           "pilot", "calibrate", "summarize", "view", "export", "paper", "readiness"])
    parser.add_argument("--config")
    parser.add_argument("--output")
    parser.add_argument("--model", default="work/archie/model.rds")
    parser.add_argument("--stage", choices=["core", "expansion", "stress"], default="core")
    parser.add_argument("--destination", default="reports/research-v2")
    args = parser.parse_args()
    research_commands = {"pilot", "calibrate", "benchmark", "summarize", "view", "export", "paper"}
    empirical_commands = {"readiness", "empirical-prepare", "empirical-run"}
    default_config = "configs/research.yaml" if args.command in research_commands else "configs/empirical.yaml" if args.command in empirical_commands else "configs/smoke.yaml"
    default_output = "results/research-v2" if args.command in research_commands else "results/empirical-readiness" if args.command == "readiness" else "results/empirical" if args.command in empirical_commands else "results/smoke"
    config = read_config(args.config or default_config)
    output = Path(args.output or default_output)
    if args.command == "view" and args.output is None and not output.exists():
        output = Path("reports/nov14-mvp/explorer-results")
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
        from .pipeline import run
        if args.command == "smoke" and not Path(args.model).exists():
            train(config, Path(args.model).parent)
        if args.command == "benchmark":
            from .study import run_study
            run_study(config, output, args.model, args.stage)
        else:
            run(config, output, args.model)
    elif args.command == "evaluate":
        from .pipeline import evaluate
        from .registry import active_dir
        output = active_dir(output)
        if (output.parent.parent / "manifest.json").exists():
            raise ValueError("Completed managed runs are immutable; inspect or export saved scores instead")
        config = json.loads((output / "config.json").read_text())
        evaluate(output, config)
    elif args.command == "report":
        from .registry import active_dir
        from .report import create_report
        output = active_dir(output)
        create_report(output)
    elif args.command == "pilot":
        from .convergence import pilot
        from .study import Budget
        output.mkdir(parents=True, exist_ok=True)
        with Budget(output, config.get("budget_seconds", 43200)) as budget:
            pilot(config | {"pilot_budget_seconds": min(config.get("pilot_budget_seconds", 7200), budget.remaining())}, output / "pilot")
    elif args.command == "calibrate":
        from .study import calibrate_only
        calibrate_only(config, output, args.model)
    elif args.command == "summarize":
        from .analysis import summarize
        summarize(output)
    elif args.command == "export":
        from .figures import export_study
        export_study(output, args.destination)
    elif args.command == "paper":
        from .paper import write_paper
        write_paper(output, args.destination)
    elif args.command == "readiness":
        from .readiness import audit
        audit(config, output)
    elif args.command == "view":
        import os
        import subprocess
        import sys

        from .registry import ROOT
        subprocess.run([sys.executable, "-m", "streamlit", "run", str(ROOT / "src/ghostbench/explorer.py"),
                        "--server.address", "127.0.0.1", "--server.headless", "true",
                        "--browser.gatherUsageStats", "false"],
                       env=os.environ | {"GHOSTBENCH_RESULTS": str(output.resolve())}, check=True)
    else:
        from .empirical import prepare, run_empirical
        if args.command == "empirical-prepare":
            prepare(config, output)
        else:
            run_empirical(config, output, args.model)


if __name__ == "__main__":
    main()
