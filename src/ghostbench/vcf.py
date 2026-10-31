"""Minimal phased biallelic VCF reader, with strict sample and position checks."""

import gzip

import numpy as np


def read_vcf(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    positions, rows, names = [], [], None
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith("#CHROM"):
                names = line.rstrip().split("\t")[9:]
            elif not line.startswith("#"):
                fields = line.rstrip().split("\t")
                if len(fields[3]) != 1 or len(fields[4]) != 1 or "," in fields[4]:
                    raise ValueError("Only biallelic SNPs are supported")
                formats = fields[8].split(":")
                gt_index = formats.index("GT")
                row = []
                for sample in fields[9:]:
                    genotype = sample.split(":")[gt_index]
                    if genotype not in {"0|0", "0|1", "1|0", "1|1"}:
                        raise ValueError("Expected completely called, phased diploid genotypes")
                    row.extend(int(allele) for allele in genotype.split("|"))
                positions.append(int(fields[1]) - 1)
                rows.append(row)
    if names is None:
        raise ValueError("Missing VCF sample header")
    positions = np.asarray(positions, dtype=int)
    if len(positions) > 1 and np.any(np.diff(positions) <= 0):
        raise ValueError("VCF positions must be strictly increasing within the input region")
    return positions, np.asarray(rows, dtype=np.int8).reshape(-1, len(names) * 2), names
