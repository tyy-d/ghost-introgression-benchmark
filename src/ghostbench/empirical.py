"""Prepare a bounded Ag1000G pilot with explicit QC and relatedness evidence.

The public API provides phased genotypes, geographic metadata and accessibility.
QC/relatedness tables must be provided when the selected API release does not
include the required evidence. Unknown quality or relatedness is never passed.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from .io import write_json


def select_samples(metadata, config):
    required = {"sample_id", "country", "taxon", "qc_pass", "unrelated_keep"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Missing QC/relatedness metadata columns: {sorted(missing)}. "
                         "Provide qc_metadata with documented qc_pass and unrelated_keep flags.")
    if metadata.sample_id.duplicated().any():
        raise ValueError("Sample identifiers must be unique")
    eligible = metadata[(metadata.country == config["country"])
                        & metadata.qc_pass.eq(True) & metadata.unrelated_keep.eq(True)]
    selected = []
    for species in [config["target_species"], config["reference_species"]]:
        candidates = eligible[eligible.taxon == species].sort_values("sample_id")
        offset = config.get("reference_offset", 0) if species == config["reference_species"] else 0
        candidates = candidates.iloc[offset:offset + config["samples_per_population"]]
        if len(candidates) < config["samples_per_population"]:
            raise ValueError(f"Insufficient unrelated quality-filtered samples for {species}")
        selected.extend(candidates.sample_id.tolist())
    return selected


def mask_intervals(mask):
    """Convert a zero-based dense accessibility mask to half-open intervals."""
    padded = np.concatenate(([False], np.asarray(mask, dtype=bool), [False]))
    changes = np.diff(padded.astype(int))
    return list(zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)))


def prepare(config, output):
    import malariagen_data

    output = Path(output)
    inputs = output / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    api = malariagen_data.Ag3()
    metadata = api.sample_metadata(sample_sets=config["sample_sets"])
    if config.get("qc_metadata"):
        evidence = pd.read_csv(config["qc_metadata"], sep="\t")
        metadata = metadata.merge(evidence, on="sample_id", validate="one_to_one")
    selected = select_samples(metadata, config)
    # Public coordinates are one-based inclusive; internal coordinates are
    # zero-based half-open [10,000,000, 12,000,000).
    region = f"{config['chromosome']}:{config['start'] + 1}-{config['end']}"
    dataset = api.haplotypes(region=region, analysis=config["analysis"],
                             sample_sets=config["sample_sets"],
                             sample_query=f"sample_id in {selected!r}")
    sample_ids = dataset.sample_id.values.astype(str).tolist()
    order = [sample_ids.index(sample) for sample in selected]
    genotypes = dataset.call_genotype.values[:, order, :]
    positions = dataset.variant_position.values.astype(int) - 1
    alleles = dataset.variant_allele.values.astype(str)
    accessible = np.asarray(api.is_accessible(region=region, site_mask=config["site_mask"]))
    if len(accessible) != config["end"] - config["start"]:
        raise ValueError("Accessibility mask length does not match requested interval")
    local = positions - config["start"]
    valid = (genotypes >= 0).all(axis=(1, 2)) & (genotypes <= 1).all(axis=(1, 2))
    valid &= np.array([len(row[0]) == 1 and len(row[1]) == 1 for row in alleles])
    # Excluded variant positions are not silently counted as invariant callable bases.
    accessible[local[~valid]] = False
    keep = accessible[local] & valid
    keep &= (genotypes == 0).any(axis=(1, 2)) & (genotypes == 1).any(axis=(1, 2))
    with (inputs / "sequences.vcf").open("w") as stream:
        stream.write("##fileformat=VCFv4.2\n")
        stream.write(f"##contig=<ID=1,length={len(accessible)}>\n")
        stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        stream.write("\t".join(selected) + "\n")
        for index in np.flatnonzero(keep):
            calls = [f"{a}|{b}" for a, b in genotypes[index]]
            fields = ["1", str(local[index] + 1), ".", alleles[index, 0], alleles[index, 1],
                      ".", "PASS", ".", "GT", *calls]
            stream.write("\t".join(fields) + "\n")
    rows = []
    for index, sample in enumerate(selected):
        for phase in range(2):
            rows.append({"node": 2 * index + phase, "haplotype": f"h{2 * index + phase}",
                         "individual": sample, "phase": phase,
                         "population": "target" if index < config["samples_per_population"] else "reference"})
    pd.DataFrame(rows).to_csv(inputs / "samples.tsv", sep="\t", index=False)
    intervals = mask_intervals(accessible)
    pd.DataFrame([{"chromosome": "1", "start": int(a), "end": int(b)} for a, b in intervals]).to_csv(
        inputs / "callable.tsv", sep="\t", index=False
    )
    metadata.set_index("sample_id").loc[selected].to_csv(output / "selected_samples.tsv", sep="\t")
    write_json(output / "coordinate_map.json", {
        "chromosome": config["chromosome"], "offset": config["start"],
        "rule": "genomic_coordinate = local_coordinate + offset; intervals are half-open",
        "ancestral_polarization": "Reference alleles are not ancestral alleles; see empirical-run guard",
    })
    write_json(output / "config.json", config)


def run_empirical(config, output, model):
    from . import archie, trace_runner
    from .evaluation import overlap_length
    from .pipeline import provenance

    required = ["mutation_rate", "recombination_rate", "effective_size", "trace_time", "parameter_source"]
    if any(config.get(key) is None for key in required):
        raise ValueError("Document species-specific rates, Ne, TRACE timescale, and parameter_source first")
    if not config.get("ancestral_polarization_verified"):
        raise ValueError("ArchIE requires ancestral/derived allele orientation. Supply an ancestrally "
                         "polarized input VCF and set ancestral_polarization_verified with a documented source.")
    output = Path(output)
    metadata = pd.read_csv(output / "inputs/samples.tsv", sep="\t")
    mask = pd.read_csv(output / "inputs/callable.tsv", sep="\t")
    intervals = list(zip(mask.start, mask.end))
    settings = config | {"sequence_length": config["end"] - config["start"],
                         "seed": 42, "singer_samples": 50, "singer_burnin": 50, "singer_thin": 20}
    archie.predict(output / "inputs/sequences.vcf", metadata, settings, model, output / "archie")
    trees, singer_time = trace_runner.infer_arg(output / "inputs/sequences.vcf", settings, output / "singer")
    trace_time = trace_runner.predict(trees, metadata, settings, output / "trace", callable_intervals=intervals)
    write_json(output / "runtime.json", {"singer_seconds": singer_time, "trace_hmm_seconds": trace_time})
    first = pd.read_csv(output / "archie/predictions.tsv", sep="\t")
    second = pd.read_csv(output / "trace/predictions.tsv", sep="\t")
    thresholds = config.get("frozen_thresholds")
    if thresholds is None:
        raise ValueError("Provide validation-frozen thresholds before summarizing empirical overlap")
    intersection = union = 0
    for haplotype in metadata.query("population == 'target'").haplotype:
        a = first[(first.haplotype == haplotype) & (first.probability >= thresholds["archie"])]
        b = second[(second.haplotype == haplotype) & (second.probability >= thresholds["trace"])]
        a_intervals, b_intervals = list(zip(a.start, a.end)), list(zip(b.start, b.end))
        from .evaluation import intersect_intervals
        shared = intersect_intervals(a_intervals, b_intervals)
        intersection += sum(overlap_length(left, right, intervals) for left, right in shared)
        from .evaluation import merged_intervals
        union += sum(overlap_length(left, right, intervals)
                     for left, right in merged_intervals(a_intervals + b_intervals))
    write_json(output / "overlap.json", {"shared_callable_bases": intersection,
                                         "union_callable_bases": union,
                                         "jaccard": intersection / union if union else None,
                                         "interpretation": "Agreement is not accuracy"})
    provenance(output)
