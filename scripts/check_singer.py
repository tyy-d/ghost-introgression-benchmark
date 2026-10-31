"""Small executable check, explicitly NOT the full simulation benchmark.

Uses a fixed toy VCF to verify that the locally built upstream C++ program runs.
It does not provide ancestry labels, invoke TRACE, or satisfy the two-method run.
"""

import argparse
import random
import sys
from pathlib import Path

from ghostbench.io import run_command, write_json


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default=".tools/singer/singer_master")
    parser.add_argument("--output", default="results/singer-build-check")
    args = parser.parse_args()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    vcf = output / "toy.vcf"
    with vcf.open("w") as stream:
        stream.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n")
        stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n')
        stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\n")
        random_generator = random.Random(42)
        for index in range(1000):
            count = random_generator.choice([1, 2, 3, 4])
            derived = set(random_generator.sample(range(8), count))
            pattern = [phase in derived for phase in range(8)]
            genotypes = [f"{int(pattern[2*i])}|{int(pattern[2*i+1])}" for i in range(4)]
            stream.write("\t".join(["1", str(50 + index * 90), ".", "A", "G",
                                    ".", "PASS", ".", "GT", *genotypes]) + "\n")
    wrapper = Path(args.binary).name == "singer_master"
    command = ([sys.executable] if wrapper else []) + [str(Path(args.binary).resolve()),
               "-vcf" if wrapper else "-input", str(vcf.with_suffix("")),
               "-output", str(output / "toy"), "-m", "1e-8",
               "-ratio" if wrapper else "-r", "1" if wrapper else "1e-8",
               "-Ne", "100000", "-start", "0", "-end", "100000", "-n", "2",
               "-thin", "1", "-seed", "42", "-polar", "0.99"]
    seconds = run_command(command, output / "singer.log")
    nodes = list(output.glob("toy_nodes_*.txt"))
    if len(nodes) != 2:
        raise RuntimeError("SINGER did not produce the requested two samples")
    write_json(output / "summary.json", {"kind": "SINGER executable-only check", "seconds": seconds,
                                         "posterior_samples": len(nodes),
                                         "full_benchmark_passed": False})


if __name__ == "__main__":
    main()
