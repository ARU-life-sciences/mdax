#!/usr/bin/env python3
"""
Extract a curated set of MDA foldback reads from a real FASTQ and independently
characterise each one.

Outputs:
  chris_artefacts.fasta        — the raw sequences (input to mdax)
  chris_artefacts_truth.tsv    — independent ground-truth per read

Independent characterisation (no mdax code used):
  split_pos  — position of best palindromic junction found by Hamming scan
  arm_len    — palindromic arm length from BioPython global alignment
  identity   — alignment matches / alignment length (from BioPython)
  gap_est    — alignment offset × 2 (template-switch gap estimate)

Strategy:
  1. Hamming scan: for each candidate split s (step 5), count complementary
     pairs seq[s-1-i] / seq[s+i] over min(s, n-s) bases.  Score = matches²/len
     so longer arms beat short perfect ones.  Peak = rough split.
  2. Fine scan ±200 bp around the rough peak (step 1).
  3. BioPython global alignment of left_arm vs RC(right_arm) at the refined
     split to get precise identity and gap estimate.

Usage:
  python3 extract_and_characterise.py \\
      --fastq ../chris_data/ERR12263839.fastq.gz \\
      --tsv   ../chris_data/ERR12263839.tsv \\
      --out-fasta chris_artefacts.fasta \\
      --out-truth chris_artefacts_truth.tsv
"""

import argparse
import gzip
import sys
from pathlib import Path

import subprocess
import tempfile

COMP = bytes.maketrans(b'ACGTacgtNn', b'TGCAtgcaNn')

def revcomp(seq: bytes) -> bytes:
    return seq.translate(COMP)[::-1]

# ---------------------------------------------------------------------------
# Read IDs to extract
# ---------------------------------------------------------------------------

TARGETS = [
    # -----------------------------------------------------------------------
    # Artefacts — identity ≥ 0.95  (support_n=1 in prior run)
    # -----------------------------------------------------------------------
    ("art_hi_3820bp",   "ERR12263839.4269 m64125e_231028_235749/329581/ccs"),
    ("art_hi_5653bp",   "ERR12263839.251 m64125e_231028_235749/964/ccs"),
    ("art_hi_6156bp",   "ERR12263839.69 m64125e_231028_235749/290/ccs"),
    ("art_hi_9193bp",   "ERR12263839.925 m64125e_231028_235749/66191/ccs"),
    ("art_hi_11012bp",  "ERR12263839.25 m64125e_231028_235749/110/ccs"),
    ("art_hi_13019bp",  "ERR12263839.8 m64125e_231028_235749/33/ccs"),
    # -----------------------------------------------------------------------
    # Artefacts — identity 0.80–0.95
    # -----------------------------------------------------------------------
    ("art_md_3728bp",   "ERR12263839.33522 m64125e_231028_235749/2950656/ccs"),
    ("art_md_5167bp",   "ERR12263839.1124 m64125e_231028_235749/66965/ccs"),
    ("art_md_5944bp",   "ERR12263839.905 m64125e_231028_235749/66129/ccs"),
    ("art_md_6356bp",   "ERR12263839.1898 m64125e_231028_235749/132505/ccs"),
    ("art_md_6437bp",   "ERR12263839.1176 m64125e_231028_235749/67141/ccs"),
    ("art_md_9605bp",   "ERR12263839.114 m64125e_231028_235749/449/ccs"),
    ("art_md_10782bp",  "ERR12263839.30 m64125e_231028_235749/134/ccs"),
    ("art_md_14746bp",  "ERR12263839.134 m64125e_231028_235749/540/ccs"),
    # -----------------------------------------------------------------------
    # Artefacts — identity 0.65–0.80  (includes the user's example at 0.661)
    # -----------------------------------------------------------------------
    ("art_lo_5094bp",   "ERR12263839.1426 m64125e_231028_235749/68064/ccs"),
    ("art_lo_6432bp",   "ERR12263839.1175 m64125e_231028_235749/67140/ccs"),  # user example
    ("art_lo_8105bp",   "ERR12263839.179 m64125e_231028_235749/712/ccs"),
    ("art_lo_8708bp",   "ERR12263839.130 m64125e_231028_235749/503/ccs"),
    ("art_lo_10715bp",  "ERR12263839.143 m64125e_231028_235749/569/ccs"),
    ("art_lo_13413bp",  "ERR12263839.17 m64125e_231028_235749/88/ccs"),
    # -----------------------------------------------------------------------
    # Artefacts — identity 0.50–0.65
    # -----------------------------------------------------------------------
    ("art_vlo_3906bp",  "ERR12263839.197518 m64125e_231028_235749/19071853/ccs"),
    ("art_vlo_5230bp",  "ERR12263839.757 m64125e_231028_235749/65562/ccs"),
    ("art_vlo_6497bp",  "ERR12263839.1897 m64125e_231028_235749/132502/ccs"),
    ("art_vlo_9989bp",  "ERR12263839.1174 m64125e_231028_235749/67137/ccs"),
    ("art_vlo_10708bp", "ERR12263839.16 m64125e_231028_235749/74/ccs"),
    ("art_vlo_13247bp", "ERR12263839.213 m64125e_231028_235749/807/ccs"),
    ("art_vlo_14291bp", "ERR12263839.28 m64125e_231028_235749/130/ccs"),
    # -----------------------------------------------------------------------
    # Low-identity foldbacks (identity < 0.50 → decision=low_ident)
    # These are detected but fall below min_identity; with cut_low_ident=false
    # they pass through unchanged in the output FASTA.
    # -----------------------------------------------------------------------
    ("li_3279bp",   "ERR12263839.135 m64125e_231028_235749/542/ccs"),
    ("li_5619bp",   "ERR12263839.192 m64125e_231028_235749/744/ccs"),
    ("li_6296bp",   "ERR12263839.67 m64125e_231028_235749/282/ccs"),
    ("li_7152bp",   "ERR12263839.133 m64125e_231028_235749/531/ccs"),
    ("li_9340bp",   "ERR12263839.19 m64125e_231028_235749/96/ccs"),
    ("li_11367bp",  "ERR12263839.80 m64125e_231028_235749/314/ccs"),
    ("li_12397bp",  "ERR12263839.12 m64125e_231028_235749/43/ccs"),
    ("li_14193bp",  "ERR12263839.33 m64125e_231028_235749/147/ccs"),
    # -----------------------------------------------------------------------
    # Real palindromes — high identity (support_n ≥ 10)
    # -----------------------------------------------------------------------
    ("real_hi_4085bp",  "ERR12263839.695008 m64125e_231028_235749/62458161/ccs"),
    ("real_hi_4203bp",  "ERR12263839.889341 m64125e_231028_235749/79235406/ccs"),
    ("real_hi_4469bp",  "ERR12263839.893556 m64125e_231028_235749/79627339/ccs"),
    ("real_hi_4797bp",  "ERR12263839.893902 m64125e_231028_235749/79628666/ccs"),
    ("real_hi_4942bp",  "ERR12263839.233629 m64125e_231028_235749/22219073/ccs"),
    ("real_hi_5421bp",  "ERR12263839.146530 m64125e_231028_235749/14157119/ccs"),
    ("real_hi_5879bp",  "ERR12263839.2980 m64125e_231028_235749/199424/ccs"),
    ("real_hi_8673bp",  "ERR12263839.20402 m64125e_231028_235749/1771300/ccs"),
    # -----------------------------------------------------------------------
    # Real palindromes — low/medium identity (support_n ≥ 3)
    # -----------------------------------------------------------------------
    ("real_lo_5920bp",  "ERR12263839.79675 m64125e_231028_235749/7408339/ccs"),
    ("real_lo_5954bp",  "ERR12263839.58164 m64125e_231028_235749/5310758/ccs"),
    ("real_lo_6526bp",  "ERR12263839.45243 m64125e_231028_235749/4065031/ccs"),
    ("real_lo_7250bp",  "ERR12263839.509 m64125e_231028_235749/1929/ccs"),
    ("real_lo_7284bp",  "ERR12263839.37580 m64125e_231028_235749/3343348/ccs"),
    ("real_lo_9927bp",  "ERR12263839.12885 m64125e_231028_235749/1115149/ccs"),
    ("real_lo_10434bp", "ERR12263839.1562 m64125e_231028_235749/131208/ccs"),
    ("real_lo_10612bp", "ERR12263839.3605 m64125e_231028_235749/264431/ccs"),
    ("real_lo_10650bp", "ERR12263839.1194 m64125e_231028_235749/67191/ccs"),
    ("real_lo_11491bp", "ERR12263839.1652 m64125e_231028_235749/131562/ccs"),
    ("real_lo_11727bp", "ERR12263839.8054 m64125e_231028_235749/657769/ccs"),
    ("real_lo_12897bp", "ERR12263839.16730 m64125e_231028_235749/1443507/ccs"),
    ("real_lo_13009bp", "ERR12263839.17422 m64125e_231028_235749/1508870/ccs"),
    ("real_lo_14005bp", "ERR12263839.1683 m64125e_231028_235749/131687/ccs"),
]


# ---------------------------------------------------------------------------
# BLAST-based independent characterisation
# ---------------------------------------------------------------------------
# Run each read against its own reverse complement using blastn.
# For a foldback [A][RC(A)], the best minus-strand self-hit gives:
#   qend (1-based)  = junction / split position
#   pident          = arm alignment identity
#   length          = palindromic arm length
#
# This is fully independent of mdax: it uses a completely different
# algorithm (seeded k-mer extension + Smith-Waterman) to locate the arm.

BLASTN     = "/Users/mc9148/Documents/software/bin/blastn"
BLAST_PERC = 55    # minimum % identity for BLAST to report a hit
BLAST_WSIZE = 11   # word size (default for blastn; 7 for very short arms)
MIN_ARM_LEN = 200  # discard hits shorter than this (adapter/primer noise)


def _run_blastn_self(seq: bytes) -> list[dict]:
    """
    BLAST seq against itself with -strand minus to find inverted repeats.

    For a foldback [A][RC(A)] at junction P, the left arm A = seq[P-L:P]
    aligns in minus-strand orientation to seq[P:P+L]:
      qstart = P-L+1  qend = P      (1-based, forward left arm)
      sstart = P+L    send = P+1    (1-based, minus strand: sstart > send)

    So:
      split_pos (0-based) = qend          (qend 1-based = 0-based split index
                                           because 1-based P = 0-based P as
                                           the first right-arm base is at 0-based P)
      gap_est             = send - 1 - qend  (0 for clean foldback)

    Hits are filtered by arm length and sorted longest-first (the longest
    inverted repeat is the main MDA foldback arm, not a short hairpin).
    """
    seq_str = seq.decode('ascii', errors='replace')

    with tempfile.TemporaryDirectory() as td:
        fa = f"{td}/read.fa"
        with open(fa, 'w') as f:
            f.write(f">read\n{seq_str}\n")

        cmd = [
            BLASTN,
            "-query",   fa,
            "-subject", fa,
            "-outfmt",  "6 qstart qend sstart send pident length",
            "-dust",    "no",
            "-word_size",     str(BLAST_WSIZE),
            "-perc_identity", str(BLAST_PERC),
            "-strand",  "minus",   # find inverted repeats within the read
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)

    hits = []
    for line in result.stdout.splitlines():
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        qstart, qend, sstart, send = (int(parts[i]) for i in range(4))
        pident = float(parts[4])
        length = int(parts[5])
        if length < MIN_ARM_LEN:
            continue
        hits.append(dict(qstart=qstart, qend=qend, sstart=sstart,
                         send=send, pident=pident, length=length))

    hits.sort(key=lambda h: h['length'], reverse=True)
    return hits


def characterise(label: str, seq: bytes) -> dict:
    """
    Find the palindromic junction and measure arm identity using BLAST.

    For a foldback [A at P-L..P][RC(A) at P..P+L], blastn -strand minus
    reports the inverted repeat as two symmetric hits:
      Hit type A  (left arm as query):
        qstart=P-L+1  qend=P      sstart=P+L  send=P+1
      Hit type B  (right arm as query, the mirror image):
        qstart=P+1    qend=P+L    sstart=P    send=P-L+1

    In BOTH cases the split (junction) is the midpoint of the gap between
    the end of one arm and the start of the other:
        split_pos = (qend + send - 1) // 2       [≡ P for both hit types]
        gap_est   = send - 1 - qend              [bp between arms]

    BLAST also reports the alignment identity over the matched arm length.
    We take the LONGEST hit (most of the arm aligned) for robustness.
    """
    n = len(seq)
    hits = _run_blastn_self(seq)
    if not hits:
        return {"label": label, "length": n, "split_pos": n // 2,
                "arm_len": 0, "identity": 0.0, "gap_est": 0, "blast_hits": 0}

    best = hits[0]
    qend = best['qend']
    send = best['send']

    # Midpoint formula works for both hit types and handles gaps between arms
    split_pos = (qend + send - 1) // 2
    gap_est   = max(0, send - 1 - qend)
    arm_len   = best['length']
    identity  = best['pident'] / 100.0

    return {
        "label":      label,
        "length":     n,
        "split_pos":  split_pos,
        "arm_len":    arm_len,
        "identity":   round(identity, 4),
        "gap_est":    gap_est,
        "blast_hits": len(hits),
    }


# ---------------------------------------------------------------------------
# FASTQ extraction + caching
# ---------------------------------------------------------------------------

def open_maybe_gz(path: str):
    p = Path(path)
    return gzip.open(path, 'rb') if p.suffix == '.gz' else open(path, 'rb')


def extract_reads(fastq_path: str, target_ids: set) -> dict:
    """Stream through a (possibly gzipped) FASTQ and collect target sequences."""
    found = {}
    with open_maybe_gz(fastq_path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            if not header.startswith(b'@'):
                continue
            rid = header[1:].rstrip().decode('utf-8', errors='replace')
            if rid in target_ids:
                found[rid] = seq.upper()
                if len(found) == len(target_ids):
                    break
    return found


def load_extracted_cache(cache_path: str) -> dict:
    """Read a small uncompressed FASTQ written by save_extracted_cache."""
    found = {}
    with open(cache_path, 'rb') as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq  = fh.readline().rstrip()
            fh.readline()
            fh.readline()
            if not header.startswith(b'@'):
                continue
            rid = header[1:].rstrip().decode('utf-8', errors='replace')
            found[rid] = seq.upper()
    return found


def save_extracted_cache(cache_path: str, seqs: dict) -> None:
    """Write extracted reads to a small uncompressed FASTQ for reuse."""
    with open(cache_path, 'w') as f:
        for rid, seq in seqs.items():
            seq_str = seq.decode('ascii', errors='replace')
            f.write(f"@{rid}\n{seq_str}\n+\n{'I' * len(seq_str)}\n")
    print(f"  cached {len(seqs)} reads → {cache_path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq',     required=True,
                    help='full (possibly gzipped) FASTQ — only read if --extracted cache is absent')
    ap.add_argument('--tsv',       required=True,
                    help='prior mdax TSV — used only to cross-check decision')
    ap.add_argument('--extracted', default='chris_artefacts_raw.fastq',
                    help='cache file for the extracted reads (created on first run, '
                         'reused on subsequent runs to avoid re-scanning the big FASTQ)')
    ap.add_argument('--out-fasta', default='chris_artefacts.fasta')
    ap.add_argument('--out-truth', default='chris_artefacts_truth.tsv')
    args = ap.parse_args()

    target_map = {full_id: label for label, full_id in TARGETS}
    target_ids = set(target_map.keys())

    cache = Path(args.extracted)
    if cache.exists():
        print(f"Loading extracted reads from cache {cache} ...", file=sys.stderr)
        seqs = load_extracted_cache(str(cache))
        print(f"  loaded {len(seqs)} reads", file=sys.stderr)
    else:
        print(f"Extracting {len(target_ids)} reads from {args.fastq} ...", file=sys.stderr)
        seqs = extract_reads(args.fastq, target_ids)
        print(f"  found {len(seqs)}/{len(target_ids)}", file=sys.stderr)
        save_extracted_cache(str(cache), seqs)

    missing = target_ids - set(seqs.keys())
    if missing:
        print("WARNING: reads not found:", file=sys.stderr)
        for rid in sorted(missing):
            print(f"  {rid}", file=sys.stderr)

    # Load prior mdax decisions for cross-check (not used as truth)
    mdax_decision = {}
    mdax_split    = {}
    with open(args.tsv) as f:
        hdr  = f.readline().strip().split('\t')
        ridx = hdr.index('read_id')
        didx = hdr.index('decision')
        sidx = hdr.index('refined_split')
        for line in f:
            row = line.strip().split('\t')
            if len(row) <= didx:
                continue
            mdax_decision[row[ridx]] = row[didx]
            mdax_split[row[ridx]]    = row[sidx]

    truth_rows  = []
    fasta_lines = []

    for label, full_id in TARGETS:
        if full_id not in seqs:
            print(f"  SKIP {label} (not found in FASTQ)", file=sys.stderr)
            continue

        seq = seqs[full_id]
        print(f"  characterising {label} ({len(seq)} bp) ...",
              file=sys.stderr, end='', flush=True)
        ch = characterise(label, seq)

        mdax_dec = mdax_decision.get(full_id, 'unknown')
        mdax_sp  = mdax_split.get(full_id, '')

        print(
            f"  split={ch['split_pos']:5d}  arm={ch['arm_len']:4d}  "
            f"ident={ch['identity']:.3f}  gap={ch['gap_est']:4d}"
            f"  [prior mdax: {mdax_dec} split={mdax_sp}]",
            file=sys.stderr
        )

        truth_rows.append(ch | {"mdax_prior_decision": mdax_dec, "read_id": full_id})
        fasta_lines.append(f">{label}  {full_id}")
        fasta_lines.append(seq.decode('ascii'))

    # Write FASTA
    out_fasta = Path(args.out_fasta)
    out_fasta.write_text('\n'.join(fasta_lines) + '\n')
    print(f"\nWrote {out_fasta}", file=sys.stderr)

    # Write truth TSV
    out_truth = Path(args.out_truth)
    cols = ['label', 'read_id', 'length', 'split_pos', 'arm_len',
            'identity', 'gap_est', 'blast_hits', 'mdax_prior_decision']
    with open(out_truth, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for row in truth_rows:
            f.write('\t'.join(str(row[c]) for c in cols) + '\n')
    print(f"Wrote {out_truth}", file=sys.stderr)


if __name__ == '__main__':
    main()
