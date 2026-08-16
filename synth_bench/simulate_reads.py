#!/usr/bin/env python3
"""Simulate a HiFi-like read set over a synth_bench genome.

Two read classes are produced:

  clean   -- a contiguous extraction from the genome (possibly spanning a
             real inverted-repeat locus baked in by make_genome.py). Start
             positions are sampled with per-locus coverage-tier weighting,
             so real repeat loci get uneven depth just like real MDA
             libraries -- this is the thing that stresses mdax's
             cross-read support-count heuristic.

  chimera -- a synthetic MDA-style foldback artefact: [lead][A][gap][RC(A)][trail],
             built independently of the genome's real repeat loci (source
             sequence is unique, single-copy), with identity/gap drawn from
             classes modelled on real observed MDA artefacts (see
             data/README.md in the main repo: art_hi/md/lo/vlo, li_*).

Outputs (into --outdir):
  reads.fa
  reads_truth.tsv
"""
import argparse
import json
import os
import random
import sys

import numpy as np

BASES = "ACGT"
BASES_ARR = np.frombuffer(BASES.encode(), dtype=np.uint8)
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def load_fasta_single(path):
    seq_chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_chunks.append(line.strip())
    return "".join(seq_chunks)


def load_bed(path):
    rows = []
    with open(path) as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(
                locus_id=parts[0], kind=parts[1],
                start=int(parts[2]), end=int(parts[3]),
                arm_len=int(parts[4]), gap_len=int(parts[5]),
                identity=float(parts[6]), cov_tier=float(parts[7]),
            ))
    return rows


def load_coverage_profile(path, genome_len):
    w = np.ones(genome_len, dtype=np.float64)
    with open(path) as f:
        header = f.readline()
        for line in f:
            s, e, weight = line.rstrip("\n").split("\t")
            w[int(s):int(e)] = float(weight)
    return w


def apply_errors(seq: str, rate: float, rng: np.random.Generator) -> str:
    """Vectorised uniform-random substitution error model (HiFi: indels negligible)."""
    if rate <= 0 or not seq:
        return seq
    arr = np.frombuffer(seq.encode(), dtype=np.uint8).copy()
    n = len(arr)
    mask = rng.random(n) < rate
    n_sub = int(mask.sum())
    if n_sub == 0:
        return seq
    # pick a replacement base index 0..2 that skips the current base, then map
    repl_choice = rng.integers(0, 3, size=n_sub)
    idx = np.nonzero(mask)[0]
    for k, i in enumerate(idx):
        cur = arr[i]
        opts = [b for b in BASES_ARR if b != cur]
        arr[i] = opts[repl_choice[k]]
    return arr.tobytes().decode()


def mutate_to_identity(s: str, identity: float, rng: random.Random) -> str:
    if identity >= 1.0:
        return s
    sub_rate = 1.0 - identity
    out = list(s)
    for i in range(len(out)):
        if rng.random() < sub_rate:
            choices = [b for b in BASES if b != out[i]]
            out[i] = rng.choice(choices)
    return "".join(out)


CHIMERA_CLASSES = [
    # (name, identity_range, gap_range, weight) -- modelled on data/README.md
    ("hi", (0.995, 1.0), (0, 0), 0.30),
    ("md", (0.985, 0.995), (0, 0), 0.25),
    ("lo", (0.950, 0.985), (0, 0), 0.20),
    ("vlo", (0.900, 0.950), (0, 0), 0.10),
    ("li", (0.75, 0.82), (1600, 4000), 0.15),
]


def pick_chimera_class(rng: random.Random):
    names = [c[0] for c in CHIMERA_CLASSES]
    weights = [c[3] for c in CHIMERA_CLASSES]
    name = rng.choices(names, weights=weights, k=1)[0]
    cls = next(c for c in CHIMERA_CLASSES if c[0] == name)
    identity = rng.uniform(*cls[1])
    gap = rng.randint(*cls[2]) if cls[2][1] > 0 else 0
    return name, identity, gap


def weighted_starts(weight_arr: np.ndarray, n: int, rng: np.random.Generator) -> np.ndarray:
    cum = np.cumsum(weight_arr)
    total = cum[-1]
    draws = rng.random(n) * total
    return np.searchsorted(cum, draws, side="right")


def overlapping_loci(bed_rows, start, end):
    hits = []
    for r in bed_rows:
        if not (end <= r["start"] or start >= r["end"]):
            hits.append(r)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--coverage-profile", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--target-coverage", type=float, default=30.0)
    ap.add_argument("--mean-read-len", type=int, default=12000)
    ap.add_argument("--min-read-len", type=int, default=3000)
    ap.add_argument("--max-read-len", type=int, default=30000)
    ap.add_argument("--error-rate", type=float, default=0.001,
                     help="uniform per-base substitution rate (HiFi-like)")
    ap.add_argument("--chimera-frac", type=float, default=0.20,
                     help="fraction of read 'slots' that become MDA foldback chimeras")
    ap.add_argument("--chimera-arm-lens", type=int, nargs="+",
                     default=[1000, 2000, 3000, 4500, 6000])
    ap.add_argument("--chimera-lead-max", type=int, default=1500)
    ap.add_argument("--chimera-trail-max", type=int, default=6000)
    ap.add_argument("--repeat-buffer", type=int, default=3000,
                     help="keep chimera source sequence at least this far from any real "
                          "repeat locus, so ground truth stays unambiguous")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng_py = random.Random(args.seed)
    rng_np = np.random.default_rng(args.seed)

    genome = load_fasta_single(args.genome)
    L = len(genome)
    bed_rows = load_bed(args.bed)
    weight_arr = load_coverage_profile(args.coverage_profile, L)

    # forbid chimera source starts near any real repeat locus
    forbidden = np.zeros(L, dtype=bool)
    for r in bed_rows:
        s = max(0, r["start"] - args.repeat_buffer)
        e = min(L, r["end"] + args.repeat_buffer)
        forbidden[s:e] = True
    allowed_starts = np.nonzero(~forbidden)[0]

    n_reads = int(args.target_coverage * L / args.mean_read_len)
    print(f"[simulate_reads] genome={L:,}bp target_cov={args.target_coverage}x "
          f"-> ~{n_reads:,} read slots", file=sys.stderr)

    # read lengths: lognormal clipped
    mu = np.log(args.mean_read_len)
    lens = rng_np.lognormal(mean=mu, sigma=0.35, size=n_reads)
    lens = np.clip(lens, args.min_read_len, args.max_read_len).astype(int)

    is_chimera = rng_np.random(n_reads) < args.chimera_frac
    starts = weighted_starts(weight_arr, n_reads, rng_np)
    strands = rng_np.integers(0, 2, size=n_reads)

    os.makedirs(args.outdir, exist_ok=True)
    fa_path = os.path.join(args.outdir, "reads.fa")
    truth_path = os.path.join(args.outdir, "reads_truth.tsv")

    n_written = 0
    n_chim_written = 0
    with open(fa_path, "w") as fa, open(truth_path, "w") as tt:
        tt.write("read_id\ttype\tgenome_start\tgenome_end\tstrand\t"
                  "overlap_loci\tfully_spans_locus\tchimera_class\t"
                  "true_split_in_read\ttrue_identity\ttrue_gap_len\tread_len\n")

        for i in range(n_reads):
            rlen = int(lens[i])
            strand = "+" if strands[i] == 0 else "-"

            if not is_chimera[i]:
                start = int(starts[i])
                if start + rlen > L:
                    start = max(0, L - rlen)
                end = start + rlen
                seq = genome[start:end]
                if strand == "-":
                    seq = revcomp(seq)
                seq = apply_errors(seq, args.error_rate, rng_np)

                hits = overlapping_loci(bed_rows, start, end)
                overlap_ids = ",".join(h["locus_id"] for h in hits) if hits else "NA"
                fully = any(start <= h["start"] and end >= h["end"] for h in hits)

                read_id = f"clean_{i}"
                tt.write(f"{read_id}\tclean\t{start}\t{end}\t{strand}\t"
                         f"{overlap_ids}\t{int(fully)}\tNA\tNA\tNA\tNA\t{len(seq)}\n")
            else:
                arm_len = rng_py.choice(args.chimera_arm_lens)
                cls_name, identity, gap_len = pick_chimera_class(rng_py)

                # pick a source start from allowed (non-repeat) positions
                idx = rng_py.randrange(len(allowed_starts))
                s = int(allowed_starts[idx])
                needed = arm_len + gap_len + arm_len + args.chimera_trail_max
                if s + needed > L:
                    s = max(0, L - needed)

                lead_len = rng_py.randint(0, args.chimera_lead_max)
                lead_start = max(0, s - lead_len)
                lead = genome[lead_start:s]

                arm = genome[s:s + arm_len]
                gap_seq = genome[s + arm_len: s + arm_len + gap_len] if gap_len else ""
                arm2 = mutate_to_identity(revcomp(arm), identity, rng_py)

                # Trailing flank: sampled from an *independent* genomic locus, not
                # contiguous with the arm's own source region. If it were left
                # contiguous (genome[s+arm_len+gap_len:...]), a reference aligner
                # can "bridge over" the inverted arm as a single large insertion
                # and chain lead+trail into one colinear alignment -- silently
                # defeating any reference-based chimera detector and making this
                # simulation unrealistically easy. Real template-switch resumption
                # has no reason to resume exactly where the forward arm left off.
                trail_len = rng_py.randint(0, args.chimera_trail_max)
                trail_idx = rng_py.randrange(len(allowed_starts))
                trail_start = int(allowed_starts[trail_idx])
                if trail_start + trail_len > L:
                    trail_start = max(0, L - trail_len)
                trail = genome[trail_start:trail_start + trail_len]

                seq = lead + arm + gap_seq + arm2 + trail
                split_in_read = len(lead) + arm_len  # true junction (left-arm end)
                seq = apply_errors(seq, args.error_rate, rng_np)

                read_id = f"chimera_{i}_{cls_name}"
                # genome_start/end here describe the *fold-arm* source locus only
                # (the trail comes from an unrelated locus, so a single
                # genome_end would be misleading).
                tt.write(f"{read_id}\tchimera\t{s}\t{s + arm_len + gap_len}\t{strand}\t"
                         f"NA\t0\t{cls_name}\t{split_in_read}\t{identity:.4f}\t{gap_len}\t{len(seq)}\n")
                n_chim_written += 1

            fa.write(f">{read_id}\n")
            for k in range(0, len(seq), 80):
                fa.write(seq[k:k + 80] + "\n")
            n_written += 1

    meta_path = os.path.join(args.outdir, "sim_meta.json")
    with open(meta_path, "w") as f:
        json.dump(dict(n_reads=n_written, n_chimeras=n_chim_written, args=vars(args)), f, indent=2)

    print(f"[simulate_reads] wrote {n_written:,} reads ({n_chim_written:,} chimeras, "
          f"{100*n_chim_written/max(1,n_written):.1f}%) -> {fa_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
