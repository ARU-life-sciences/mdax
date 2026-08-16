#!/usr/bin/env python3
"""Evaluate one or more assemblies against the synth_bench truth genome.

For each assembly this reports:
  - basic contig stats: N50, L50, #contigs, total length, longest contig
  - genome fraction covered (union of aligned reference bases / genome length)
  - per real-inverted-repeat-locus outcome: spanned / broken / missing
      spanned -- a single contig alignment block fully contains the locus
                 (+ flank) => the assembler resolved straight through it
      broken  -- the locus is covered, but no single alignment block spans
                 it end-to-end => something truncated a read/contig inside
                 the repeat
      missing -- no alignment overlaps the locus at all => that genomic
                 region never made it into any contig

This is the core hypothesis test: if mdax's cross-read support-count gate
is misclassifying real repeats as artefacts and cutting them, the
mdax-corrected arm should show a much higher broken/missing rate at real
inverted-repeat loci (especially the low-coverage-tier ones) than the raw
or reference-based arms -- and that should track directly with lower N50.

Usage:
  evaluate.py --truth-genome genome.fa --repeats-bed repeats_truth.bed \
      --assembly raw=raw.p_ctg.fa --assembly mdax_balanced=mdax_balanced.p_ctg.fa \
      ... --outdir eval_out --minimap2 /path/to/minimap2
"""
import argparse
import os
import subprocess
import sys


def load_fasta_lengths(path):
    lens = {}
    name = None
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = n
                name = line[1:].strip().split()[0]
                n = 0
            else:
                n += len(line.strip())
    if name is not None:
        lens[name] = n
    return lens


def contig_stats(lens: dict):
    sizes = sorted(lens.values(), reverse=True)
    total = sum(sizes)
    n_contigs = len(sizes)
    longest = sizes[0] if sizes else 0
    n50 = l50 = 0
    running = 0
    for i, s in enumerate(sizes, start=1):
        running += s
        if running >= total / 2 and n50 == 0:
            n50 = s
            l50 = i
    return dict(n_contigs=n_contigs, total_len=total, longest=longest, n50=n50, l50=l50)


def load_bed(path):
    rows = []
    with open(path) as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rows.append(dict(
                locus_id=parts[0], kind=parts[1],
                start=int(parts[2]), end=int(parts[3]),
                arm_len=int(parts[4]), gap_len=int(parts[5]),
                identity=float(parts[6]), cov_tier=float(parts[7]),
            ))
    return rows


def run_minimap2_paf(minimap2, ref, query, out_paf, threads):
    if os.path.exists(out_paf) and os.path.getsize(out_paf) > 0:
        return
    with open(out_paf, "w") as out:
        subprocess.run(
            [minimap2, "-x", "asm20", "--secondary=no", "-t", str(threads), ref, query],
            stdout=out, check=True, stderr=subprocess.DEVNULL,
        )


def load_paf_blocks(paf_path, min_block_len=500, min_mapq=1):
    """Return list of (rstart, rend, qname) primary alignment blocks."""
    blocks = []
    with open(paf_path) as f:
        for line in f:
            f_ = line.rstrip("\n").split("\t")
            qname = f_[0]
            rstart, rend = int(f_[7]), int(f_[8])
            mapq = int(f_[11])
            if rend - rstart < min_block_len or mapq < min_mapq:
                continue
            blocks.append((rstart, rend, qname))
    return blocks


def genome_fraction_covered(blocks, genome_len):
    if not blocks:
        return 0.0
    ivs = sorted((b[0], b[1]) for b in blocks)
    covered = 0
    cur_s, cur_e = ivs[0]
    for s, e in ivs[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            covered += cur_e - cur_s
            cur_s, cur_e = s, e
    covered += cur_e - cur_s
    return covered / genome_len


def classify_loci(blocks, bed_rows, flank=300):
    """spanned / broken / missing per real inverted-repeat locus."""
    out = {}
    ir_rows = [r for r in bed_rows if r["kind"] == "inverted_repeat"]
    for r in ir_rows:
        lo, hi = r["start"] - flank, r["end"] + flank
        spanned = any(b[0] <= lo and b[1] >= hi for b in blocks)
        if spanned:
            out[r["locus_id"]] = "spanned"
            continue
        overlapping = [b for b in blocks if not (b[1] <= r["start"] or b[0] >= r["end"])]
        out[r["locus_id"]] = "broken" if overlapping else "missing"
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--truth-genome", required=True)
    ap.add_argument("--repeats-bed", required=True)
    ap.add_argument("--assembly", action="append", required=True,
                     help="NAME=PATH, repeatable")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    genome_len = sum(load_fasta_lengths(args.truth_genome).values())
    bed_rows = load_bed(args.repeats_bed)
    ir_loci = [r for r in bed_rows if r["kind"] == "inverted_repeat"]

    arms = []
    for spec in args.assembly:
        name, path = spec.split("=", 1)
        arms.append((name, path))

    summary_rows = []
    per_locus_matrix = {r["locus_id"]: {} for r in ir_loci}

    for name, path in arms:
        if not os.path.exists(path):
            print(f"[evaluate] WARNING: {name} assembly not found at {path}, skipping", file=sys.stderr)
            continue
        lens = load_fasta_lengths(path)
        stats = contig_stats(lens)

        paf_path = os.path.join(args.outdir, f"{name}.vs_truth.paf")
        run_minimap2_paf(args.minimap2, args.truth_genome, path, paf_path, args.threads)
        blocks = load_paf_blocks(paf_path)
        gfrac = genome_fraction_covered(blocks, genome_len)
        locus_class = classify_loci(blocks, bed_rows)

        n_spanned = sum(1 for v in locus_class.values() if v == "spanned")
        n_broken = sum(1 for v in locus_class.values() if v == "broken")
        n_missing = sum(1 for v in locus_class.values() if v == "missing")

        for locus_id, cls in locus_class.items():
            per_locus_matrix[locus_id][name] = cls

        summary_rows.append(dict(
            arm=name, **stats, genome_frac_covered=gfrac,
            ir_spanned=n_spanned, ir_broken=n_broken, ir_missing=n_missing,
            ir_total=len(ir_loci),
        ))

    summary_path = os.path.join(args.outdir, "summary.tsv")
    with open(summary_path, "w") as f:
        cols = ["arm", "n_contigs", "total_len", "longest", "n50", "l50",
                "genome_frac_covered", "ir_spanned", "ir_broken", "ir_missing", "ir_total"]
        f.write("\t".join(cols) + "\n")
        for row in summary_rows:
            f.write("\t".join(str(row[c]) for c in cols) + "\n")

    locus_path = os.path.join(args.outdir, "per_locus.tsv")
    arm_names = [name for name, _ in arms]
    with open(locus_path, "w") as f:
        f.write("locus_id\tarm_len\tgap_len\tidentity\tcov_tier\t" + "\t".join(arm_names) + "\n")
        by_id = {r["locus_id"]: r for r in ir_loci}
        for locus_id in sorted(per_locus_matrix, key=lambda x: by_id[x]["start"]):
            r = by_id[locus_id]
            vals = [per_locus_matrix[locus_id].get(a, "NA") for a in arm_names]
            f.write(f"{locus_id}\t{r['arm_len']}\t{r['gap_len']}\t{r['identity']}\t{r['cov_tier']}\t"
                    + "\t".join(vals) + "\n")

    # ---- human-readable summary ----
    print(f"\n{'arm':<20}{'n_ctg':>8}{'N50':>12}{'L50':>6}{'total_len':>12}"
          f"{'genome%':>10}{'IR spanned':>12}{'broken':>8}{'missing':>9}")
    for row in summary_rows:
        print(f"{row['arm']:<20}{row['n_contigs']:>8}{row['n50']:>12,}{row['l50']:>6}"
              f"{row['total_len']:>12,}{100*row['genome_frac_covered']:>9.1f}%"
              f"{row['ir_spanned']:>12}{row['ir_broken']:>8}{row['ir_missing']:>9}")
    print(f"\n[evaluate] wrote {summary_path}")
    print(f"[evaluate] wrote {locus_path}")


if __name__ == "__main__":
    main()
