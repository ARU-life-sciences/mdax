#!/usr/bin/env python3
"""What mdax leaves behind, and how well it places its cuts, measured against
a reference.

Inputs (all for the same random read sample):
  TRUTH_SAM   raw sampled reads aligned to the reference (cols 1-6, primary +
              supplementary) -- the same file refbased_recall_eval.py scores
  CORR_SAM    the *corrected* versions of those reads, aligned the same way
  REPORT      mdax's report TSV rows (any superset of the sample is fine)
  RAW_FA      raw sampled reads
  CORR_FA     corrected sampled reads (ids without the _chopped_at_ suffix)

Sections printed:
  1. how many true foldback junctions a read has
  2. reads that still carry a foldback after correction, by cause
  3. breakpoint precision (refined split, and where the cut actually fell)
  4. the identity estimate vs an independent arm identity at the true junction
  5. unique sequence discarded with the dropped side
"""
import collections
import csv
import re
import sys

truth_sam, corr_sam, report, raw_fa, corr_fa = sys.argv[1:6]
MIN_SEG = 200
FOLD_REF_GAP = 5000
cig = re.compile(r"(\d+)([MIDNSHP=X])")


def load_sam(path):
    segs = collections.defaultdict(list)
    rlen = {}
    seen = set()
    for line in open(path):
        if line[0] == "@":
            continue
        f = line.rstrip("\n").split("\t", 6)
        flag = int(f[1])
        q = f[0]
        seen.add(q)
        if flag & 4 or flag & 256:
            continue
        ops = [(int(n), o) for n, o in cig.findall(f[5])]
        lead = ops[0][0] if ops[0][1] in "SH" else 0
        trail = ops[-1][0] if ops[-1][1] in "SH" else 0
        ql = sum(n for n, o in ops if o in "MI=X")
        rl = sum(n for n, o in ops if o in "MDN=X")
        rev = bool(flag & 16)
        qs = trail if rev else lead
        rs = int(f[3]) - 1
        segs[q].append((qs, qs + ql, f[2], rs, rs + rl, rev))
        rlen[q] = lead + ql + trail
    out = {}
    for q in seen:
        ss = sorted(s for s in segs.get(q, []) if s[1] - s[0] >= MIN_SEG)
        out[q] = ss
    return out, rlen


def fold_junctions(ss):
    """Adjacent segment pairs that form a foldback: same chromosome, opposite
    strand, reference intervals overlapping or within FOLD_REF_GAP."""
    js = []
    for a, b in zip(ss, ss[1:]):
        if a[2] == b[2] and a[5] != b[5]:
            gap = max(a[3], b[3]) - min(a[4], b[4])
            if gap < FOLD_REF_GAP:
                js.append((a, b))
    return js


def junction_geometry(a, b):
    """Read-coordinate junction and the fold's displacement.

    Returns (left_end, right_start, d, spacer_side): the left arm ends at
    left_end, the right arm starts at right_start (read coordinates), and d is
    how far apart on the reference the two arms' junction-proximal ends are.
    spacer_side says which arm carries those d un-mirrored bases next to the
    junction ('L', 'R', or '-' when d is ~0).
    """
    xa = a[3] if a[5] else a[4]  # ref coordinate where the left arm stops
    xb = b[4] if b[5] else b[3]  # ref coordinate where the right arm starts
    d = abs(xa - xb)
    if d < 5:
        side = "-"
    elif not a[5]:  # left arm ascending: right arm starting short of xa => spacer on the left
        side = "L" if xb < xa else "R"
    else:  # left arm descending
        side = "L" if xb > xa else "R"
    return a[1], b[0], d, side


def load_fa(path):
    d = {}
    name = None
    buf = []
    for line in open(path):
        if line[0] == ">":
            if name is not None:
                d[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name is not None:
        d[name] = "".join(buf)
    return d


COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def banded_ed(a, b, band):
    """Edit distance between a and b with |i-j| <= band (capped at band+1)."""
    n, m = len(a), len(b)
    if abs(n - m) > band:
        return band + 1
    INF = 10 ** 9
    prev = {j: j for j in range(0, min(m, band) + 1)}
    for i in range(1, n + 1):
        cur = {}
        lo, hi = max(0, i - band), min(m, i + band)
        ai = a[i - 1]
        for j in range(lo, hi + 1):
            if j == 0:
                cur[j] = i
                continue
            best = prev.get(j - 1, INF) + (ai != b[j - 1])
            v = prev.get(j, INF) + 1
            if v < best:
                best = v
            v = cur.get(j - 1, INF) + 1
            if v < best:
                best = v
            cur[j] = best
        prev = cur
    return min(prev.get(m, INF), band + 1)


def pct(xs, ps=(5, 25, 50, 75, 95)):
    if not xs:
        return "n=0"
    xs = sorted(xs)
    return "  ".join(f"p{p}={xs[min(len(xs) - 1, int(p / 100 * len(xs)))]}" for p in ps)


def share(n, d):
    return f"{100 * n / d:5.1f}%" if d else "  n/a"


# ---------------------------------------------------------------- load
tsegs, trlen = load_sam(truth_sam)
csegs, _ = load_sam(corr_sam)
rows = {}
for r in csv.DictReader(open(report), delimiter="\t"):
    q = r["read_id"].split()[0]
    if q in tsegs:
        rows[q] = r
raw = load_fa(raw_fa)
corr = load_fa(corr_fa)

tj = {q: fold_junctions(ss) for q, ss in tsegs.items()}
fold_reads = [q for q, js in tj.items() if js]
N = len(tsegs)
print(f"sample reads={N}  true-foldback reads={len(fold_reads)} ({share(len(fold_reads), N)})")

# ---------------------------------------------------------------- 1
print("\n[1] true foldback junctions per read")
h = collections.Counter(min(len(tj[q]), 4) for q in fold_reads)
for k in sorted(h):
    print(f"  {k}{'+' if k == 4 else ' '} junction(s): {h[k]:7d} {share(h[k], len(fold_reads))} of true-foldback reads")
multi = [q for q in fold_reads if len(tj[q]) > 1]
print(f"  >1 junction: {len(multi)} = {share(len(multi), len(fold_reads))} of true-foldback reads, {share(len(multi), N)} of all reads")
seg3 = [q for q, ss in tsegs.items() if len(ss) >= 3]
print(f"  reads with >=3 aligned segments: {len(seg3)} ({share(len(seg3), N)} of all); of those with >=1 foldback junction: "
      f"{sum(1 for q in seg3 if tj[q])}, with >=2: {sum(1 for q in seg3 if len(tj[q]) > 1)}")

# ---------------------------------------------------------------- kept ranges
keep = {}
unlocated = 0
for q in tsegs:
    r, c = raw.get(q), corr.get(q)
    if r is None or c is None:
        continue
    if len(c) == len(r):
        keep[q] = (0, len(r))
        continue
    i = r.find(c)
    if i < 0:
        unlocated += 1
        continue
    keep[q] = (i, i + len(c))
cut_reads = [q for q, (s, e) in keep.items() if (s, e) != (0, len(raw[q]))]
print(f"\nreads with raw+corrected sequence: {len(keep)}; cut by mdax: {len(cut_reads)}; corrected read not locatable in raw: {unlocated}")

# ---------------------------------------------------------------- 2
print("\n[2] foldbacks surviving correction (corrected read re-aligned, same foldback test)")
still = [q for q in csegs if fold_junctions(csegs[q])]
print(f"  corrected reads still carrying a foldback: {len(still)} = {share(len(still), N)} of all reads")
orig_fold_still = [q for q in still if tj.get(q)]
print(f"  of the {len(fold_reads)} true-foldback reads, still carrying one: {len(orig_fold_still)} = {share(len(orig_fold_still), len(fold_reads))}")
print(f"  (foldback only visible after correction, i.e. raw read was not classed foldback: {len(still) - len(orig_fold_still)})")
cause = collections.Counter()
resid_arm = collections.defaultdict(list)
for q in orig_fold_still:
    r = rows.get(q)
    if r is None:
        c = "never detected"
    elif r["decision"] == "real":
        c = "kept as real"
    elif r["decision"] in ("low_ident", "unknown"):
        c = "detected, not cut (" + r["decision"] + ")"
    elif q in keep and keep[q] == (0, len(raw[q])):
        c = "called artefact but not actually cut"
    elif len(tj[q]) > 1:
        c = "cut, read had >1 junction (another junction left)"
    else:
        c = "cut, single junction (residual of the same fold)"
    cause[c] += 1
    js = fold_junctions(csegs[q])
    resid_arm[c].append(max(min(a[1] - a[0], b[1] - b[0]) for a, b in js))
for c, n in cause.most_common():
    print(f"    {c:52s} {n:6d} {share(n, len(orig_fold_still))} of surviving, {share(n, len(fold_reads))} of true-foldback reads"
          f"   residual shorter-arm bp: {pct(resid_arm[c], (25, 50, 75))}")
# among multi-junction reads specifically
m_cut = [q for q in multi if q in rows and rows[q]["decision"] == "artefact"]
m_still = [q for q in m_cut if q in csegs and fold_junctions(csegs[q])]
print(f"  multi-junction reads mdax cut: {len(m_cut)}; still carrying a foldback afterwards: {len(m_still)} ({share(len(m_still), len(m_cut))})")
s_cut = [q for q in fold_reads if len(tj[q]) == 1 and q in rows and rows[q]["decision"] == "artefact"]
s_still = [q for q in s_cut if q in csegs and fold_junctions(csegs[q])]
print(f"  single-junction reads mdax cut: {len(s_cut)}; still carrying a foldback afterwards: {len(s_still)} ({share(len(s_still), len(s_cut))})")

# ---------------------------------------------------------------- 3
print("\n[3] breakpoint precision, true foldbacks mdax cut (nearest true junction to mdax's split)")
split_err, cut_err, cut_err_by = [], [], collections.defaultdict(list)
geom = {}
for q in fold_reads:
    r = rows.get(q)
    if r is None or r["decision"] != "artefact" or q not in keep:
        continue
    sp = int(r["refined_split"])
    best = None
    for a, b in tj[q]:
        le, rs_, d, side = junction_geometry(a, b)
        # the symmetry axis mdax aims at: junction midpoint, shifted by half the spacer
        axis = (le + rs_) / 2 + (-d / 2 if side == "L" else d / 2 if side == "R" else 0)
        if best is None or abs(sp - axis) < abs(sp - best[0]):
            best = (axis, le, rs_, d, side, a, b)
    axis, le, rs_, d, side, a, b = best
    geom[q] = best
    split_err.append(int(round(sp - axis)))
    ks, ke = keep[q]
    L = len(raw[q])
    if (ks, ke) == (0, L):
        continue
    # Chimeric bases retained (+) or good bases lost (-) at the cut, relative
    # to the true arm boundary on the side that was kept. Only meaningful for
    # the first cut at this junction; reads cut on both sides are skipped.
    if ks == 0 and ke < L:
        e = ke - le
        kind = "kept-left"
    elif ke == L and ks > 0:
        e = rs_ - ks
        kind = "kept-right"
    else:
        continue
    if abs(e) > 5000:
        continue  # the cut that produced this boundary was at a different junction
    cut_err.append(e)
    cut_err_by[(kind, "spacer-" + side)].append(e)


def within(xs, t):
    return share(sum(1 for x in xs if abs(x) <= t), len(xs))


print(f"  refined_split - true symmetry axis (bp), n={len(split_err)}: {pct(split_err)}")
print(f"     |err|<=10: {within(split_err, 10)}  <=50: {within(split_err, 50)}  <=200: {within(split_err, 200)}  <=1000: {within(split_err, 1000)}")
print(f"  actual cut vs true arm boundary on the kept side (+ = chimeric bases retained, - = good bases lost), n={len(cut_err)}: {pct(cut_err)}")
print(f"     |err|<=10: {within(cut_err, 10)}  <=50: {within(cut_err, 50)}  <=200: {within(cut_err, 200)}")
print(f"     retains >50 bp of chimeric sequence: {share(sum(1 for x in cut_err if x > 50), len(cut_err))}   "
      f">200 bp: {share(sum(1 for x in cut_err if x > 200), len(cut_err))}   loses >200 bp of good sequence: {share(sum(1 for x in cut_err if x < -200), len(cut_err))}")
for k in sorted(cut_err_by):
    xs = cut_err_by[k]
    print(f"     {k[0]:10s} {k[1]:9s} n={len(xs):6d} {pct(xs)}   retains>50bp: {share(sum(1 for x in xs if x > 50), len(xs))}")
kept_bases = sum(keep[q][1] - keep[q][0] for q in geom if q in keep)
print(f"     chimeric bases retained next to the cut (sum of + errors): {sum(x for x in cut_err if x > 0)} = "
      f"{share(sum(x for x in cut_err if x > 0), kept_bases)} of the bases kept from cut true-foldback reads")
ds = [g[3] for g in geom.values()]
print(f"  fold displacement d (bp) at the true junction: {pct(ds)}   d<5: {share(sum(1 for x in ds if x < 5), len(ds))}  d>200: {share(sum(1 for x in ds if x > 200), len(ds))}")
gap_diff = [int(rows[q]['gap_est']) - geom[q][3] for q in geom]
print(f"  gap_est - true displacement (bp): {pct(gap_diff)}   |diff|<=50: {within(gap_diff, 50)}")

# ---------------------------------------------------------------- 4
print("\n[4] identity_est for true foldbacks mdax detected")
ids = collections.Counter()


def ib(x):
    return "<0.60" if x < .6 else "0.60-0.90" if x < .9 else "0.90-0.97" if x < .97 else "0.97-0.99" if x < .99 else "0.99-0.995" if x < .995 else ">=0.995"


order = ["<0.60", "0.60-0.90", "0.90-0.97", "0.97-0.99", "0.99-0.995", ">=0.995"]
det = [q for q in fold_reads if q in rows]
for q in det:
    ids[ib(float(rows[q]["identity_est"]))] += 1
print("  " + "  ".join(f"{k}:{share(ids[k], len(det))}" for k in order) + f"   (n={len(det)})")

# independent identity at the reference-derived junction
W = 300
MAXN = 1200
recs = []
for q in sorted(geom)[:MAXN]:
    axis, le, rs_, d, side, a, b = geom[q]
    r = raw.get(q)
    if r is None:
        continue
    l_end = le - (d if side == "L" else 0)
    r_start = rs_ + (d if side == "R" else 0)
    w = min(W, l_end, len(r) - r_start, a[1] - a[0] - (d if side == "L" else 0), b[1] - b[0] - (d if side == "R" else 0))
    if w < 100:
        continue
    left = r[l_end - w:l_end].upper()
    # The alignment-derived arm boundaries are only good to a few bp, and a
    # 2 bp misregistration alone costs ~1% identity over a few hundred bp, so
    # take the best of a small range of right-arm offsets (cheap narrow band),
    # and only fall back to one wide-band comparison for divergent arms.
    ed = 10 ** 9
    for off in range(-6, 7):
        rs2 = r_start + off
        if rs2 < 0 or rs2 + w > len(r):
            continue
        e2 = banded_ed(left, revcomp(r[rs2:rs2 + w]).upper(), 24)
        if e2 < ed:
            ed = e2
    if ed > 24:
        ed = banded_ed(left, revcomp(r[r_start:r_start + w]).upper(), max(24, int(0.3 * w)))
    true_id = 1 - ed / w
    row = rows[q]
    recs.append((q, float(row["identity_est"]), true_id, int(row["refined_split"]) - axis, int(row["gap_est"]), d, side, w, len(r), int(row["refined_split"])))
if recs:
    tid = collections.Counter(ib(x[2]) for x in recs)
    print(f"  independent arm identity at the true junction ({W} bp windows, n={len(recs)}):")
    print("  " + "  ".join(f"{k}:{share(tid[k], len(recs))}" for k in order))
    low = [x for x in recs if x[1] < 0.995]
    print(f"  reads with identity_est < 0.995: {len(low)} ({share(len(low), len(recs))})")
    hi_true = [x for x in low if x[2] >= 0.995]
    print(f"    of those, independent identity >= 0.995 (estimator under-reads a pristine fold): {len(hi_true)} ({share(len(hi_true), len(low))})")
    cat = collections.Counter()
    for x in hi_true:
        _, est, tru, serr, gap, d, side, w, L, sp = x
        near_end = min(sp, L - sp) < 1200 + 200
        if abs(gap - d) > 100:
            cat["gap_est off by >100 bp from true displacement"] += 1
        elif abs(serr) > 100:
            cat["split off by >100 bp from true axis"] += 1
        elif near_end:
            cat["split within refine arm+window of a read end (clamped arm)"] += 1
        else:
            cat["split and gap both close; window content differs"] += 1
    for k, v in cat.most_common():
        print(f"      {k:62s} {v:6d} {share(v, len(hi_true))}")
    # The clear-cut failures: estimate says "poor" where the arms are in fact near-identical.
    bad = [x for x in recs if x[1] < 0.90 and x[2] >= 0.97]
    print(f"    identity_est < 0.90 while independent identity >= 0.97: {len(bad)} ({share(len(bad), len(recs))} of scored reads); "
          f"independent identity < 0.90: {sum(1 for x in recs if x[2] < 0.90)} ({share(sum(1 for x in recs if x[2] < 0.90), len(recs))})")
    cat = collections.Counter()
    for x in bad:
        _, est, tru, serr, gap, d, side, w, L, sp = x
        if abs(gap - d) > 100:
            cat["gap_est off by >100 bp from true displacement"] += 1
        elif abs(serr) > 100:
            cat["split off by >100 bp from true axis"] += 1
        elif min(sp, L - sp) < 1400:
            cat["split within refine arm+window of a read end (clamped arm)"] += 1
        else:
            cat["split and gap both close; window content differs"] += 1
    for k, v in cat.most_common():
        print(f"      {k:62s} {v:6d} {share(v, len(bad))}")
    under = [x[4] - x[5] for x in bad]
    print(f"      gap_est - true displacement in those reads (bp): {pct(under)}; gap_est at/near the 1000 bp scan ceiling or 0: "
          f"{share(sum(1 for x in bad if x[4] == 0 or 950 <= x[4] <= 1000), len(bad))}")
    lo_true = [x for x in low if x[2] < 0.995]
    print(f"    of those, independent identity also < 0.995 (arms genuinely differ at the junction, or junction mis-derived): {len(lo_true)} ({share(len(lo_true), len(low))})")
    print("  examples (identity_est < 0.995, sorted by estimate):")
    print("    read\tidentity_est\tindep_identity\tsplit_minus_axis\tgap_est\ttrue_d\tspacer_side\tread_len\trefined_split")
    for x in sorted(low, key=lambda x: x[1])[:6] + sorted(low, key=lambda x: x[1])[len(low) // 2:len(low) // 2 + 6]:
        print(f"    {x[0]}\t{x[1]:.3f}\t{x[2]:.3f}\t{x[3]:+.0f}\t{x[4]}\t{x[5]}\t{x[6]}\t{x[8]}\t{x[9]}")
    # does the estimate track the independent value at all?
    both = sum(1 for x in recs if x[1] >= 0.995 and x[2] >= 0.995)
    print(f"  agreement: est>=0.995 & indep>=0.995: {both}; est>=0.995 & indep<0.995: {sum(1 for x in recs if x[1] >= 0.995 and x[2] < 0.995)}; "
          f"est<0.995 & indep>=0.995: {len(hi_true)}; both <0.995: {len(lo_true)}")

# ---------------------------------------------------------------- 5
print("\n[5] sequence discarded with the dropped side")
tot_bases = sum(len(raw[q]) for q in keep)
removed = uniq = unaligned_removed = 0
uniq_reads = 0
per_read_uniq = []
for q in cut_reads:
    ks, ke = keep[q]
    L = len(raw[q])
    ss = tsegs[q]
    kept_iv = collections.defaultdict(list)
    for (qs, qe, ch, rs, re_, rev) in ss:
        s, e = max(qs, ks), min(qe, ke)
        if e <= s:
            continue
        f0, f1 = (s - qs) / (qe - qs), (e - qs) / (qe - qs)
        if rev:
            f0, f1 = 1 - f1, 1 - f0
        kept_iv[ch].append((rs + f0 * (re_ - rs), rs + f1 * (re_ - rs)))
    rem_total = L - (ke - ks)
    removed += rem_total
    rem_aln = 0
    u = 0.0
    for (qs, qe, ch, rs, re_, rev) in ss:
        for (xs, xe) in ((0, ks), (ke, L)):
            s, e = max(qs, xs), min(qe, xe)
            if e <= s:
                continue
            rem_aln += e - s
            f0, f1 = (s - qs) / (qe - qs), (e - qs) / (qe - qs)
            if rev:
                f0, f1 = 1 - f1, 1 - f0
            a0, a1 = rs + f0 * (re_ - rs), rs + f1 * (re_ - rs)
            cov = 0.0
            for (b0, b1) in sorted(kept_iv.get(ch, [])):
                o = min(a1, b1) - max(a0, b0)
                if o > 0:
                    cov += o
            cov = min(cov, a1 - a0)
            scale = (e - s) / max(1.0, a1 - a0)
            u += (a1 - a0 - cov) * scale
    unaligned_removed += rem_total - rem_aln
    uniq += u
    per_read_uniq.append(int(u))
    if u >= 500:
        uniq_reads += 1
print(f"  bases in sample: {tot_bases}; removed by mdax: {removed} ({share(removed, tot_bases)})")
print(f"  removed bases whose reference position is NOT covered by the kept side: {int(uniq)} = {share(uniq, removed)} of removed, {share(uniq, tot_bases)} of all bases")
print(f"  removed bases with no alignment at all: {unaligned_removed} ({share(unaligned_removed, removed)} of removed)")
print(f"  cut reads discarding >=500 bp of such unique sequence: {uniq_reads} ({share(uniq_reads, len(cut_reads))} of cut reads); per-read unique bp: {pct(per_read_uniq, (50, 75, 90, 95, 99))}")
