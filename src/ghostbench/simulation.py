"""Simulations with migration-record ancestry labels, isolated from detector input."""

from pathlib import Path

import msprime
import numpy as np
import pandas as pd

from .io import vcf_position, write_json


def simulate(config, seed, length=None, fraction=None):
    demography = msprime.Demography()
    for name in ["target", "reference", "donor", "modern", "ancestor"]:
        demography.add_population(name=name, initial_size=config.get("simulation_effective_size", config["effective_size"]))
    proportion = config["admixture_fraction"] if fraction is None else fraction
    if proportion:
        # Backwards in time, target lineages move into the unsampled donor.
        demography.add_mass_migration(
            time=config["admixture_time"], source="target", dest="donor", proportion=proportion
        )
    if config.get("reference_admixture_fraction", 0):
        demography.add_mass_migration(time=config["admixture_time"], source="reference", dest="donor",
                                      proportion=config["reference_admixture_fraction"])
    if config.get("ancestral_structure"):
        # Bidirectional ancestral migration without a discrete pulse. This is
        # a structured-history stress test, not a biologically unadmixed null.
        demography.add_migration_rate_change(
            time=config["reference_split"] + 1, source="modern", dest="donor", rate=1e-6
        )
        demography.add_migration_rate_change(
            time=config["reference_split"] + 1, source="donor", dest="modern", rate=1e-6
        )
    demography.add_population_split(
        time=config["reference_split"], derived=["target", "reference"], ancestral="modern"
    )
    demography.add_population_split(
        time=config["donor_split"], derived=["modern", "donor"], ancestral="ancestor"
    )
    demography.sort_events()
    ancestry = msprime.sim_ancestry(
        samples={"target": config["target_diploids"], "reference": config["reference_diploids"]},
        demography=demography, sequence_length=length or config["sequence_length"],
        recombination_rate=config["recombination_rate"], random_seed=int(seed),
        record_migrations=True,
    )
    return msprime.sim_mutations(
        ancestry, rate=config["mutation_rate"], random_seed=int(seed) + 1,
        model=msprime.MatrixMutationModel(
            alleles=["A", "G"], root_distribution=[1.0, 0.0],
            transition_matrix=[[0.0, 1.0], [1.0, 0.0]],
        ),
    )


def donor_intervals(ts, admixture_time, target_nodes, donor_population=2):
    """Identify descendants of each target-to-donor pulse migration.

    A migration node can be an ancestor of several sampled haplotypes. Traverse
    every overlapping marginal tree rather than assigning ancestry only to tips.
    Merge overlapping pieces so each base is counted at most once.
    """
    intervals = {int(node): [] for node in target_nodes}
    for migration in ts.migrations():
        if migration.dest != donor_population or not np.isclose(migration.time, admixture_time):
            continue
        for tree in ts.trees():
            left = max(tree.interval.left, migration.left)
            right = min(tree.interval.right, migration.right)
            if right <= left:
                continue
            for node in target_nodes:
                if tree.is_descendant(int(node), migration.node):
                    intervals[int(node)].append((int(left), int(right)))
    merged = {}
    for node, pieces in intervals.items():
        result = []
        for left, right in sorted(pieces):
            if result and left <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], right))
            else:
                result.append((left, right))
        merged[node] = result
    return merged


def sample_metadata(ts, target_count):
    rows = []
    for index, node in enumerate(ts.samples()):
        population = "target" if index < target_count * 2 else "reference"
        rows.append({"node": int(node), "haplotype": f"h{index}",
                     "individual": f"sample{index // 2}", "phase": index % 2,
                     "population": population})
    return pd.DataFrame(rows)


def export_vcf(ts, path):
    names = [f"sample{i}" for i in range(ts.num_individuals)]
    with Path(path).open("w") as stream:
        ts.write_vcf(stream, contig_id="1", individual_names=names,
                     position_transform=lambda positions: np.array([vcf_position(p) for p in positions]))


def generate(config, output):
    output = Path(output)
    inputs, truth = output / "inputs", output / "truth"
    inputs.mkdir(parents=True, exist_ok=True)
    truth.mkdir(parents=True, exist_ok=True)
    ts = simulate(config, config["seed"])
    metadata = sample_metadata(ts, config["target_diploids"])
    metadata.to_csv(inputs / "samples.tsv", sep="\t", index=False)
    export_vcf(ts, inputs / "sequences.vcf")
    # True ARGs and migration labels are never stored alongside inference inputs.
    ts.dump(truth / "ancestry.trees")
    targets = metadata.query("population == 'target'")
    intervals = donor_intervals(ts, config["admixture_time"], targets.node)
    rows = [{"haplotype": f"h{node}", "start": left, "end": right}
            for node, pieces in intervals.items() for left, right in pieces]
    pd.DataFrame(rows, columns=["haplotype", "start", "end"]).to_csv(
        truth / "ancestry.tsv", sep="\t", index=False
    )
    pd.DataFrame([{"chromosome": "1", "start": 0, "end": int(ts.sequence_length)}]).to_csv(
        inputs / "callable.tsv", sep="\t", index=False
    )
    write_json(output / "config.json", config)
    return ts
