#!/usr/bin/env python3
"""
probe_identity.py — per-position arm match rate diagnostic.

For each detected foldback artefact, compares base-by-base:
  seq[split - 1 - i]  vs  complement(seq[split + jump + i])
for i = 0 .. min(MAX_DIST, available)

Outputs:
  - per-position mean match rate vs distance from junction
  - estimated palindromic arm length distribution
  - short-arm (first 200bp) identity vs reported identity_est
"""
import csv, statistics

REPORT   = "out/nrCaeEleg92.report.tsv"
FASTA    = "out/nrCaeEleg92.subset.fasta"
MAX_DIST = 1200
N_READS  = 2000

RC = str.maketrans("ACGTacgt", "TGCAtgca")

def complement(b):
    return b.translate(RC)

# ── 1. load report ────────────────────────────────────────────────────────────
splits = {}
with open(REPORT) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("decision") != "artefact":
            continue
        try:
            rid      = row["read_id"].split()[0]
            sp       = int(row["refined_split"])
            jump     = int(row["gap_est"]) // 2   # gap_est = 2δ; δ is the symmetric half-shift
            ident    = float(row["identity_est"])
            arm_est  = int(row["arm_len_est"]) if row.get("arm_len_est") else None
        except (ValueError, KeyError):
            continue
        splits[rid] = (sp, jump, ident, arm_est)

print(f"Artefact reads in report: {len(splits)}", flush=True)

# ── 2. load sequences ─────────────────────────────────────────────────────────
seqs = {}
cur_id = None
cur_seq = []
with open(FASTA) as fh:
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if cur_id is not None and cur_id in splits:
                seqs[cur_id] = "".join(cur_seq)
            cur_id  = line[1:].split()[0]
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_id is not None and cur_id in splits:
        seqs[cur_id] = "".join(cur_seq)

print(f"Sequences loaded: {len(seqs)}", flush=True)

# ── 3. per-position match profile ─────────────────────────────────────────────
match_sum = [0] * MAX_DIST
match_cnt = [0] * MAX_DIST

short_idents    = []   # identity over first 200 bp (jump-corrected)
reported_idents = []
arm_lengths     = []   # estimated palindromic arm length per read

n = 0
reported_arm_lens = []  # arm_len_est from TSV (for comparison with probe's estimate)
for rid, (sp, jump, ident, arm_est) in splits.items():
    if rid not in seqs:
        continue
    if n >= N_READS:
        break
    n += 1

    seq  = seqs[rid]
    # Symmetric model: left arm ends at (sp - jump), right arm starts at (sp + jump).
    # jump = δ = gap_est/2 is the half-gap shift applied symmetrically around split_pos.
    left_end   = sp - jump          # true end of forward arm (approx.)
    right_start = sp + jump         # true start of RC arm (approx.)
    dist = min(MAX_DIST, left_end, len(seq) - right_start)

    if dist < 50:
        continue

    reported_idents.append(ident)
    if arm_est is not None:
        reported_arm_lens.append(arm_est)

    # Walk outward from the corrected edges:
    # position i is i bp outward from the arm edges (away from the gap).
    match_run = []
    for i in range(dist):
        lb = seq[left_end - 1 - i]             # i bp left of corrected junction
        rb = complement(seq[right_start + i])  # i bp right of corrected junction
        m  = int(lb == rb)
        match_sum[i] += m
        match_cnt[i] += 1
        match_run.append(m)

    short_d = min(200, dist)
    short_idents.append(sum(match_run[:short_d]) / short_d)

    # Palindromic arm estimate: first position where 50-bp rolling mean < 0.70.
    W = 50
    arm_len = dist
    for i in range(W, dist):
        win = sum(match_run[i - W : i]) / W
        if win < 0.70:
            arm_len = i - W // 2
            break
    arm_lengths.append(arm_len)

# ── 4. print results ──────────────────────────────────────────────────────────
print(f"\nReads profiled: {n}  (used in stats: {len(short_idents)})")
print(f"Reported identity_est:       "
      f"mean={statistics.mean(reported_idents):.3f}  "
      f"median={statistics.median(reported_idents):.3f}")
print(f"Short-arm identity (≤200 bp): "
      f"mean={statistics.mean(short_idents):.3f}  "
      f"median={statistics.median(short_idents):.3f}")

arm_sorted = sorted(arm_lengths)
p = lambda q: arm_sorted[int(len(arm_sorted) * q)]
print(f"Palindromic arm length (probe rolling-window): "
      f"mean={statistics.mean(arm_lengths):.0f}  "
      f"p25={p(0.25):.0f}  median={statistics.median(arm_lengths):.0f}  p75={p(0.75):.0f}  max={max(arm_lengths)}")
if reported_arm_lens:
    ra = sorted(reported_arm_lens)
    rp = lambda q: ra[int(len(ra) * q)]
    print(f"Palindromic arm length (TSV arm_len_est):      "
          f"mean={statistics.mean(reported_arm_lens):.0f}  "
          f"p25={rp(0.25):.0f}  median={statistics.median(reported_arm_lens):.0f}  p75={rp(0.75):.0f}  max={max(reported_arm_lens)}")

# ── 5. position profile ───────────────────────────────────────────────────────
print(f"\nPer-position match rate (distance from junction → match %):")
print(f"{'dist':>5}  {'rate':>5}  bar")
checkpoints = (list(range(0, 101, 5)) +
               list(range(100, 501, 25)) +
               list(range(500, MAX_DIST + 1, 100)))
for i in sorted(set(checkpoints)):
    if i >= MAX_DIST or match_cnt[i] == 0:
        continue
    rate = match_sum[i] / match_cnt[i]
    bar  = "█" * int(rate * 50)
    print(f"{i:>5}  {rate:>4.1%}  {bar}")
