#!/usr/bin/env python3
"""Generate a synthetic reference genome with controlled, ground-truthed
inverted-repeat (foldback-like) loci baked directly into the sequence.

The point: real inverted repeats living in a genome should look, to a
per-read foldback detector, exactly like a genuine "palindrome" -- because
every read that spans the locus contains the same structure, regardless of
which molecule/amplification event it came from. This script places such
loci at controlled arm length / gap / identity / local-coverage-tier, so we
can later test whether mdax's cross-read support-count heuristic correctly
recognises them as "real" across a range of local sequencing depths.

Also injects a smaller number of *tandem* (non-inverted) repeats and general
unique sequence so the assembly graph isn't trivial.

Outputs (into --outdir):
  genome.fa              -- the synthetic reference (single contig, id="chr1")
  repeats_truth.bed      -- locus_id, start, end, arm_len, gap_len, identity, cov_tier
  coverage_profile.tsv   -- per-window relative sampling weight (for simulate_reads.py)
"""
import argparse
import json
import random
import sys

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def rand_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choices(BASES, k=n))


def mutate(rng: random.Random, s: str, identity: float) -> str:
    """Introduce substitutions so that expected identity to the original ~= identity."""
    if identity >= 1.0:
        return s
    sub_rate = 1.0 - identity
    out = list(s)
    for i in range(len(out)):
        if rng.random() < sub_rate:
            choices = [b for b in BASES if b != out[i]]
            out[i] = rng.choice(choices)
    return "".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genome-len", type=int, default=2_000_000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--flank-min", type=int, default=5000,
                     help="minimum unique bp between inserted features")
    ap.add_argument(
        "--arm-lens", type=int, nargs="+", default=[800, 1500, 3000, 6000],
        help="candidate arm lengths (bp) for real inverted-repeat loci",
    )
    ap.add_argument(
        "--gap-lens", type=int, nargs="+", default=[0, 0, 50, 300],
        help="candidate spacer lengths (bp) between the two arms",
    )
    ap.add_argument(
        "--identities", type=float, nargs="+", default=[1.0, 0.99, 0.97, 0.93],
        help="candidate arm-to-arm identities (some real IRs are diverged, not perfect)",
    )
    ap.add_argument(
        "--cov-tiers", type=float, nargs="+", default=[0.15, 0.4, 1.0, 2.5],
        help="candidate relative local-coverage multipliers for repeat loci "
             "(mimics uneven MDA amplification -- some real repeats end up "
             "in low-coverage patches purely by chance)",
    )
    ap.add_argument("--n-loci-per-combo", type=int, default=3,
                     help="how many independent genomic loci per (arm_len x cov_tier) combo")
    ap.add_argument("--tandem-repeats", type=int, default=6,
                     help="number of tandem (non-inverted) duplications to also insert, "
                          "as an assembly-difficulty control unrelated to foldback detection")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    # ---- build a list of features to insert, then lay them out along the genome ----
    features = []  # each: dict(kind, arm_len, gap_len, identity, cov_tier)
    locus_id = 0
    for arm_len in args.arm_lens:
        for cov_tier in args.cov_tiers:
            for _ in range(args.n_loci_per_combo):
                gap_len = rng.choice(args.gap_lens)
                identity = rng.choice(args.identities)
                features.append(dict(
                    kind="inverted_repeat",
                    locus_id=f"ir{locus_id:03d}",
                    arm_len=arm_len, gap_len=gap_len,
                    identity=identity, cov_tier=cov_tier,
                ))
                locus_id += 1
    rng.shuffle(features)

    for i in range(args.tandem_repeats):
        unit_len = rng.choice(args.arm_lens)
        features.append(dict(
            kind="tandem_repeat",
            locus_id=f"td{i:03d}",
            arm_len=unit_len, gap_len=0,
            identity=1.0, cov_tier=1.0,
        ))

    rng.shuffle(features)

    # ---- lay out genome: unique flank, feature, unique flank, feature, ... ----
    seq_parts = []
    bed_rows = []
    cov_windows = []  # (start, end, weight)
    pos = 0

    def add_unique(n):
        nonlocal pos
        s = rand_dna(rng, n)
        seq_parts.append(s)
        cov_windows.append((pos, pos + n, 1.0))
        pos += n

    for feat in features:
        add_unique(rng.randint(args.flank_min, args.flank_min * 2))

        arm = rand_dna(rng, feat["arm_len"])
        gap = rand_dna(rng, feat["gap_len"]) if feat["gap_len"] else ""

        if feat["kind"] == "inverted_repeat":
            arm2 = mutate(rng, revcomp(arm), feat["identity"])
            block = arm + gap + arm2
        else:  # tandem_repeat: A + A (mutated copy), same orientation
            arm2 = mutate(rng, arm, feat["identity"])
            block = arm + gap + arm2

        start = pos
        seq_parts.append(block)
        cov_windows.append((pos, pos + len(block), feat["cov_tier"]))
        pos += len(block)
        end = pos

        bed_rows.append(dict(
            locus_id=feat["locus_id"], kind=feat["kind"],
            start=start, end=end,
            arm_len=feat["arm_len"], gap_len=feat["gap_len"],
            identity=feat["identity"], cov_tier=feat["cov_tier"],
        ))

    remaining = args.genome_len - pos
    if remaining > 0:
        add_unique(remaining)

    genome = "".join(seq_parts)

    # ---- write outputs ----
    import os
    os.makedirs(args.outdir, exist_ok=True)

    fa_path = os.path.join(args.outdir, "genome.fa")
    with open(fa_path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(genome), 80):
            f.write(genome[i:i + 80] + "\n")

    bed_path = os.path.join(args.outdir, "repeats_truth.bed")
    with open(bed_path, "w") as f:
        f.write("locus_id\tkind\tstart\tend\tarm_len\tgap_len\tidentity\tcov_tier\n")
        for r in sorted(bed_rows, key=lambda r: r["start"]):
            f.write(f"{r['locus_id']}\t{r['kind']}\t{r['start']}\t{r['end']}\t"
                    f"{r['arm_len']}\t{r['gap_len']}\t{r['identity']:.3f}\t{r['cov_tier']:.3f}\n")

    cov_path = os.path.join(args.outdir, "coverage_profile.tsv")
    with open(cov_path, "w") as f:
        f.write("start\tend\tweight\n")
        for s, e, w in cov_windows:
            f.write(f"{s}\t{e}\t{w}\n")

    meta_path = os.path.join(args.outdir, "genome_meta.json")
    with open(meta_path, "w") as f:
        json.dump(dict(genome_len=len(genome), n_features=len(features),
                        seed=args.seed, args=vars(args)), f, indent=2)

    n_ir = sum(1 for r in bed_rows if r["kind"] == "inverted_repeat")
    n_td = sum(1 for r in bed_rows if r["kind"] == "tandem_repeat")
    print(f"[make_genome] wrote {fa_path} ({len(genome):,} bp)", file=sys.stderr)
    print(f"[make_genome] {n_ir} inverted-repeat loci, {n_td} tandem-repeat loci", file=sys.stderr)
    print(f"[make_genome] truth bed: {bed_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
