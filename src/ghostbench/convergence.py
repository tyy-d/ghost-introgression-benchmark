"""Split-chain diagnostics and a bounded, real SINGER convergence pilot."""

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .io import write_json


def chain_diagnostics(chains):
    """Classical split R-hat and Geyer initial-positive/monotone ESS.

    These scalar summaries screen mixing, not prove ARG posterior convergence.
    Identical constants are explicitly uninformative; differing constants fail.
    """
    chains = np.asarray(chains, dtype=float)
    if chains.ndim != 2 or chains.shape[0] < 2 or chains.shape[1] < 8 or not np.isfinite(chains).all():
        return {"rhat": None, "ess": None, "status": "insufficient"}
    half = chains.shape[1] // 2
    split = np.concatenate((chains[:, :half], chains[:, -half:]), axis=0)
    m, n = split.shape
    variances = split.var(axis=1, ddof=1)
    within = float(variances.mean())
    if within == 0:
        identical = bool(np.all(split == split[0, 0]))
        return {"rhat": None, "ess": None, "status": "constant-uninformative" if identical else "nonmixing"}
    between = float(n * split.mean(axis=1).var(ddof=1))
    variance = (n-1)/n*within + between/n
    rhat = float(np.sqrt(variance / within))
    centered = split - split.mean(axis=1, keepdims=True)
    covariances = np.array([np.mean(np.sum(centered[:, :n-lag] * centered[:, lag:], axis=1)/n)
                            for lag in range(n)])
    rho = 1 - (within - covariances) / variance
    rho[0] = 1.0
    pairs = []
    for lag in range(0, n-1, 2):
        pair = float(rho[lag] + rho[lag+1])
        if pair < 0:
            break
        pairs.append(min(pair, pairs[-1]) if pairs else pair)
    tau = max(1.0, -1 + 2*sum(pairs))
    ess = float(min(m*n, m*n/tau))
    return {"rhat": rhat, "ess": ess, "status": "pass" if rhat <= 1.05 and ess >= 50 else "nonmixing"}


def diagnose_tables(tables):
    diagnostics = {name: chain_diagnostics([table[name].to_numpy() for table in tables])
                   for name in ["diversity_fit", "incompatibilities", "ancestry_fraction_h0"]}
    informative = [value for value in diagnostics.values() if value["status"] != "constant-uninformative"]
    passed = len(informative) >= 2 and all(value["status"] == "pass" for value in informative)
    return {"status": "screened" if passed else "unresolved", "diagnostics": diagnostics,
            "note": "At least two informative statistics must pass; constant statistics cannot establish mixing"}


def posterior_statistics(trees, config, destination, deadline=None):
    """Use real per-posterior TRACE calls for the predeclared focal haplotype h0."""
    import tskit

    from . import trace_runner
    from .simulation import sample_metadata

    rows = []
    for index, path in enumerate(trees):
        if deadline is not None and time.monotonic() >= deadline:
            raise TimeoutError("Pilot budget exhausted during posterior diagnostics")
        ts = tskit.load(path)
        edges = np.arange(0, ts.sequence_length, config["window_size"], dtype=float)
        edges = np.append(edges, ts.sequence_length)
        fit = float(np.mean((ts.diversity(windows=edges, mode="site")
                             - config["mutation_rate"] * ts.diversity(windows=edges, mode="branch"))**2))
        incompatible = 0
        for tree in ts.trees():
            for site in tree.sites():
                if len(site.mutations) > 2 or (len(site.mutations) == 2
                        and all(mutation.node != tree.root for mutation in site.mutations)):
                    incompatible += 1
        metadata = sample_metadata(ts, config["target_diploids"])
        # Only h0 is focal; retain every other leaf as a comparison sample.
        metadata.loc[metadata.haplotype != "h0", "population"] = "reference"
        result_dir = Path(destination) / f"posterior-{index}"
        trace_runner.predict([path], metadata, config, result_dir)
        predictions = pd.read_csv(result_dir / "predictions.tsv", sep="\t")
        ancestry = float(((predictions.end-predictions.start) * (predictions.probability >= 0.5)).sum() / ts.sequence_length)
        rows.append({"sample": index, "diversity_fit": fit, "incompatibilities": incompatible,
                     "ancestry_fraction_h0": ancestry})
    return pd.DataFrame(rows)


def pilot_identity(config):
    from .registry import CORE_FILES, ROOT, digest
    from .study import COMPUTE_KEYS

    request = {"config": {key: config[key] for key in COMPUTE_KEYS if key in config},
               "ladder": config.get("pilot_ladder", [50, 100, 200]),
               "source": {name: digest(ROOT / "src/ghostbench" / name) for name in CORE_FILES},
               "diagnostics_source": digest(Path(__file__)), "kind": "convergence-pilot-v2"}
    return hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()


def pilot(config, output):
    from .simulation import generate
    from .trace_runner import infer_arg

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    limit = config.get("pilot_budget_seconds", 7200)
    fingerprint = pilot_identity(config)
    report_path = output / "pilot.json"
    if report_path.exists():
        old = json.loads(report_path.read_text())
        if old["fingerprint"] != fingerprint:
            raise ValueError("Stale pilot configuration: choose a new output directory")
        if old.get("finished"):
            return old
        raise ValueError("Interrupted pilot requires a fresh output directory")
    report = {"fingerprint": fingerprint, "status": "unresolved", "finished": False,
              "attempts": [], "statistic_definition": "50-kb diversity MSE, unmapped sites, actual TRACE ancestry fraction for h0"}
    write_json(report_path, report)
    try:
        for retained in config.get("pilot_ladder", [50, 100, 200]):
            conditions = []
            setting = {"singer_samples": retained, "singer_burnin": retained, "singer_thin": 5}
            for scenario_index, (name, fraction) in enumerate([("baseline", 0.05), ("null", 0.0)]):
                parameters = config | setting | {"seed": 900000+1000*scenario_index, "admixture_fraction": fraction}
                base = output / name / str(retained)
                generate(parameters, base)
                tables = []
                for chain in range(2):
                    remaining = limit - (time.monotonic()-started)
                    if remaining <= 0:
                        raise TimeoutError("Convergence pilot budget exhausted")
                    settings = parameters | {"singer_seed": 950000+1000*scenario_index+chain,
                                             "job_timeout_seconds": remaining}
                    chain_dir = base / f"chain-{chain}"
                    trees, seconds = infer_arg(base / "inputs/sequences.vcf", settings, chain_dir)
                    table = posterior_statistics(trees, settings, chain_dir / "trace_diagnostics", deadline=started+limit)
                    table.to_csv(chain_dir / "diagnostics.tsv", sep="\t", index=False)
                    tables.append(table)
                    if time.monotonic()-started >= limit:
                        raise TimeoutError("Convergence pilot budget exhausted")
                result = {"scenario": name, **diagnose_tables(tables)}
                conditions.append(result)
            report["attempts"].append({"settings": setting, "conditions": conditions})
            write_json(report_path, report)
            if all(item["status"] == "screened" for item in conditions):
                report.update(status="screened", selected_settings=setting)
                break
    except TimeoutError as error:
        report["error"] = str(error)
    except KeyboardInterrupt as error:
        report.update(status="interrupted", error=f"KeyboardInterrupt: {error}")
        raise
    except BaseException as error:
        report.update(status="failed", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        report.update(finished=True, elapsed_seconds=time.monotonic()-started)
        write_json(report_path, report)
    return report
