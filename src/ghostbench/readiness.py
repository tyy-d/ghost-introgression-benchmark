"""Bounded public-data readiness audit; never substitutes reference for ancestral alleles."""

import io
import time
import urllib.request
from pathlib import Path

import pandas as pd

from .io import write_json


def audit(config, output):
    started = time.monotonic()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    report = {"status": "not ready for inference", "data": [], "sources": [
        "https://www.malariagen.net/data_package/ag1000g-phase3-hap/",
        "https://malariagen.github.io/vector-data/ag3/cloud.html",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC8742016/"],
        "biological_note": "A. coluzzii has a reported mutation-rate estimate near 1e-9/site/generation, with wide uncertainty; this is not a fitted A. gambiae demographic model",
        "gates": {"qc_relatedness_evidence": bool(config.get("qc_metadata")),
                  "ancestral_polarization": bool(config.get("ancestral_polarization_verified")),
                  "species_parameters": all(config.get(key) for key in ["mutation_rate", "recombination_rate", "effective_size", "trace_time", "parameter_source"]),
                  "frozen_thresholds": bool(config.get("frozen_thresholds"))}}
    for sample_set in config["sample_sets"]:
        url = f"https://storage.googleapis.com/vo_agam_release_master_us_central1/v3/metadata/general/{sample_set}/samples.meta.csv"
        try:
            with urllib.request.urlopen(url, timeout=15) as response:
                payload = response.read(5*1024*1024)
            table = pd.read_csv(io.BytesIO(payload))
            report["data"].append({"sample_set": sample_set, "url": url, "metadata_available": True,
                                   "samples": len(table), "columns": table.columns.tolist(),
                                   "selected_country_samples": int((table.country == config["country"]).sum())})
            table.to_csv(output / f"{sample_set}.metadata.csv", index=False)
        except Exception as error:
            report["data"].append({"sample_set": sample_set, "url": url, "metadata_available": False,
                                   "error": f"{type(error).__name__}: {error}"})
        if time.monotonic()-started > 3600:
            break
    report["elapsed_seconds"] = time.monotonic()-started
    report["note"] = "Metadata availability does not establish QC, unrelatedness, polarization, or ancestry correctness"
    write_json(output / "readiness.json", report)
    return report
