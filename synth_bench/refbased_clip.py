#!/usr/bin/env python3
"""Reference-based chimera-removal proxy (the "reference based removal" arm
the user's colleagues compared mdax against).

Approach: align every read to the TRUE synthetic genome with minimap2
(split/supplementary alignment enabled). A genuine MDA foldback chimera,
built from single-copy source sequence that is NOT duplicated anywhere in
the reference, produces a primary + supplementary alignment pair whose
reference intervals *overlap* but sit on *opposite strands* -- exactly the
signature of "this read contains a local inverted duplication that isn't
in the genome". A real inverted-repeat locus baked into the reference, by
contrast, maps as a single contiguous linear alignment (no supplementary at
all), because the fold is literally part of the linear reference sequence.

This gives a defensible, cheap stand-in for a reference-based
dechimerization step without requiring a specific external tool: it will
never touch a read that truly matches the reference, so it should leave
real repeat-spanning reads untouched by construction -- which is exactly
the property we're testing whether mdax has.

Usage:
  refbased_clip.py --genome genome.fa --reads reads.fa --outdir out/ \
      --minimap2 /path/to/minimap2 --threads 8
"""
import argparse
import os
import subprocess
import sys


def parse_cigar(cigar):
    """Return (leading_soft, trailing_soft, ref_len, query_len_incl_clips)."""
    if cigar == "*":
        return 0, 0, 0, 0
    ops = []
    num = ""
    for c in cigar:
        if c.isdigit():
            num += c
        else:
            ops.append((int(num), c))
            num = ""
    ref_len = sum(n for n, op in ops if op in "MDN=X")
    leading_soft = ops[0][0] if ops and ops[0][1] == "S" else 0
    trailing_soft = ops[-1][0] if ops and ops[-1][1] == "S" else 0
    query_len = sum(n for n, op in ops if op in "MIS=X")
    return leading_soft, trailing_soft, ref_len, query_len


def parse_sa_tag(sa):
    """SA:Z: rname,pos,strand,CIGAR,mapQ,NM;... -> list of (rname, start0, end0, strand)."""
    out = []
    for entry in sa.strip(";").split(";"):
        if not entry:
            continue
        rname, pos, strand, cigar, mapq, nm = entry.split(",")
        _, _, ref_len, _ = parse_cigar(cigar)
        start0 = int(pos) - 1
        out.append((rname, start0, start0 + ref_len, strand))
    return out


def load_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-overlap", type=int, default=50,
                     help="min bp overlap between primary/supplementary ref "
                          "intervals to call it a foldback-like chimera")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sam_path = os.path.join(args.outdir, "aln.sam")

    if not os.path.exists(sam_path) or os.path.getsize(sam_path) == 0:
        print("[refbased_clip] running minimap2...", file=sys.stderr)
        with open(sam_path, "w") as out:
            subprocess.run(
                [args.minimap2, "-a", "-x", "map-hifi", "-t", str(args.threads),
                 args.genome, args.reads],
                stdout=out, check=True,
                stderr=subprocess.DEVNULL,
            )
    else:
        print("[refbased_clip] reusing existing aln.sam", file=sys.stderr)

    reads = load_fasta(args.reads)

    decisions = {}  # read_id -> (action, reason)
    primary_seen = set()

    with open(sam_path) as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            qname, flag, rname = fields[0], int(fields[1]), fields[2]
            pos, mapq, cigar = int(fields[3]), int(fields[4]), fields[5]
            seq = fields[9]

            is_unmapped = bool(flag & 0x4)
            is_secondary = bool(flag & 0x100)
            is_supplementary = bool(flag & 0x800)
            is_reverse = bool(flag & 0x10)

            if is_secondary or is_supplementary:
                continue  # we only need the primary record (it carries the SA tag)

            primary_seen.add(qname)

            if is_unmapped:
                decisions[qname] = ("unchanged", "unmapped")
                continue

            leading_s, trailing_s, ref_len, _ = parse_cigar(cigar)
            prim_start, prim_end = pos - 1, pos - 1 + ref_len
            prim_strand = "-" if is_reverse else "+"

            sa_field = None
            for tag in fields[11:]:
                if tag.startswith("SA:Z:"):
                    sa_field = tag[5:]
                    break

            if not sa_field:
                decisions[qname] = ("unchanged", "no_supplementary")
                continue

            sa_entries = parse_sa_tag(sa_field)

            # Check *all pairs* of alignment blocks for this read (primary +
            # every supplementary), not just primary-vs-supplementary. A read
            # can have >2 blocks (e.g. two fold arms plus an unrelated
            # trailing flank that happens to be the longest/primary block),
            # in which case the two conflicting arms never involve the
            # primary directly -- only each other.
            blocks = [(rname, prim_start, prim_end, prim_strand)]
            blocks.extend(sa_entries)

            is_foldback_like = False
            for bi in range(len(blocks)):
                for bj in range(bi + 1, len(blocks)):
                    rn_i, s_i, e_i, st_i = blocks[bi]
                    rn_j, s_j, e_j, st_j = blocks[bj]
                    if rn_i != rn_j or st_i == st_j:
                        continue
                    overlap = min(e_i, e_j) - max(s_i, s_j)
                    if overlap >= args.min_overlap:
                        is_foldback_like = True
                        break
                if is_foldback_like:
                    break

            if not is_foldback_like:
                decisions[qname] = ("unchanged", "supplementary_not_foldback_like")
                continue

            # Clip to the primary alignment's aligned query span (soft-clip only,
            # since SEQ for the primary record is the full, non-hard-clipped read).
            qlen = len(seq)
            qstart = leading_s
            qend = qlen - trailing_s
            if qend <= qstart:
                decisions[qname] = ("unchanged", "degenerate_clip")
                continue
            # Orientation of the clipped sequence is irrelevant for assembly.
            clipped_seq = seq[qstart:qend]
            decisions[qname] = ("clipped", f"foldback_like;kept={qend-qstart}/{qlen}")
            reads[qname] = clipped_seq

    fa_out = os.path.join(args.outdir, "reads.refbased.fa")
    tsv_out = os.path.join(args.outdir, "refbased_decisions.tsv")

    n_clipped = 0
    n_unchanged = 0
    n_unmapped = 0
    with open(fa_out, "w") as fa, open(tsv_out, "w") as tt:
        tt.write("read_id\taction\treason\n")
        for read_id, seq in reads.items():
            action, reason = decisions.get(read_id, ("unchanged", "not_in_sam"))
            tt.write(f"{read_id}\t{action}\t{reason}\n")
            if action == "clipped":
                n_clipped += 1
            elif reason == "unmapped":
                n_unmapped += 1
            else:
                n_unchanged += 1
            fa.write(f">{read_id}\n")
            for k in range(0, len(seq), 80):
                fa.write(seq[k:k + 80] + "\n")

    print(f"[refbased_clip] {len(reads):,} reads: clipped={n_clipped:,} "
          f"unchanged={n_unchanged:,} unmapped={n_unmapped:,}", file=sys.stderr)
    print(f"[refbased_clip] wrote {fa_out}", file=sys.stderr)


if __name__ == "__main__":
    main()
