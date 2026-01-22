#!/usr/bin/env python3

import random
from typing import List, Tuple

# ----------------------------
# Basic sequence utilities
# ----------------------------

def rand_dna(n: int, seed: int | None = None) -> str:
    if seed is not None:
        random.seed(seed)
    return "".join(random.choice("ACGT") for _ in range(n))


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]


# ----------------------------
# Concatemer generators
# ----------------------------

def forward_concatemer(seq: str, n_copies: int) -> str:
    return "".join(seq for _ in range(n_copies))


def inverted_concatemer(seq: str) -> str:
    return seq + revcomp(seq)


def template_switch(
    seq: str,
    cut1: int,
    cut2: int,
    rc: bool = False,
) -> str:
    left = seq[:cut1]
    right = seq[cut2:]
    if rc:
        right = revcomp(right)
    return left + right


def multi_switch_concatemer(
    seq: str,
    switches: List[Tuple[int, int, str]],
) -> str:
    """
    switches: list of (start, end, orientation)
    orientation is "+" or "-"
    """
    out = []
    for s, e, orient in switches:
        frag = seq[s:e]
        if orient == "-":
            frag = revcomp(frag)
        out.append(frag)
    return "".join(out)


# ----------------------------
# Main driver
# ----------------------------

def main():
    # Base template (think "true molecule")
    template = rand_dna(4000, seed=42)

    reads = {}
    truth = []

    # 1. Simple forward concatemer
    r = forward_concatemer(template[:1000], 3)
    reads["forward_concatemer"] = r
    truth.append(("forward_concatemer", "forward", "2", "1000,2000"))

    # 2. Inverted concatemer (foldback)
    r = inverted_concatemer(template[:800])
    reads["inverted_concatemer"] = r
    truth.append(("inverted_concatemer", "inverted", "1", "800"))

    # 3. Template switch (strand switch)
    r = template_switch(template, cut1=1200, cut2=300, rc=True)
    reads["template_switch"] = r
    truth.append(("template_switch", "switch", "1", "1200"))

    # 4. Multi-switch MDA-style read
    switches = [
        (0, 1200, "+"),
        (400, 1600, "-"),
        (1600, 2600, "+"),
    ]
    r = multi_switch_concatemer(template, switches)
    reads["multi_switch"] = r
    truth.append(("multi_switch", "multi", "2", "1200,2400"))

    # ----------------------------
    # Write FASTA
    # ----------------------------

    with open("test_reads.fa", "w") as f:
        for name, seq in reads.items():
            f.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

    # ----------------------------
    # Write ground-truth TSV
    # ----------------------------

    with open("truth.tsv", "w") as f:
        f.write("read_id\ttype\tn_breakpoints\tapprox_positions\n")
        for row in truth:
            f.write("\t".join(row) + "\n")

    print("Wrote:")
    print("  - test_reads.fa")
    print("  - truth.tsv")


if __name__ == "__main__":
    main()

