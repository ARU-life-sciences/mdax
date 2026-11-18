#!/usr/bin/env python3
"""Reference-based truth for a RANDOM read sample (not just mdax-detected reads),
joined to mdax's report: detection recall, classification loss, non-foldback
chimera share, and a duplicate-aware support re-count."""
import sys, re, csv, collections
sam, report = sys.argv[1], sys.argv[2]
cig = re.compile(r'(\d+)([MIDNSHP=X])')
segs = collections.defaultdict(list); rlen = {}; unmapped = set(); clip = {}
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t', 11); flag = int(f[1]); q = f[0]
    if flag & 4: unmapped.add(q); continue
    if flag & 256: continue
    ops = [(int(n), o) for n, o in cig.findall(f[5])]
    lead = ops[0][0] if ops[0][1] in 'SH' else 0
    trail = ops[-1][0] if ops[-1][1] in 'SH' else 0
    ql = sum(n for n, o in ops if o in 'MI=X'); rl = sum(n for n, o in ops if o in 'MDN=X')
    tot = lead + ql + trail; rev = bool(flag & 16)
    qs = trail if rev else lead           # read-forward coordinates
    segs[q].append((qs, qs + ql, f[2], int(f[3]) - 1, int(f[3]) - 1 + rl, rev))
    rlen[q] = tot
    if not flag & 2048: clip[q] = (lead if not rev else trail, trail if not rev else lead)
truth = {}; info = {}
for q, ss in segs.items():
    ss = [s for s in ss if s[1] - s[0] >= 200]
    ss.sort()
    if len(ss) <= 1: truth[q] = 'linear'; continue
    cls = 'linear_split'; best = None
    for a, b in zip(ss, ss[1:]):
        same_chr = a[2] == b[2]
        gap = max(a[3], b[3]) - min(a[4], b[4])   # <0 => ref intervals overlap
        if same_chr and a[5] != b[5] and gap < 5000:
            arm = min(a[1] - a[0], b[1] - b[0]); junc = (a[1] + b[0]) // 2
            if best is None or arm > best[0]:
                best = (arm, junc, min(junc, rlen[q] - junc), gap < 0)
            cls = 'foldback'
        elif cls != 'foldback' and not (same_chr and a[5] == b[5] and abs(gap) < 20000):
            cls = 'other_chimera'
    truth[q] = cls
    if best: info[q] = best
for q in unmapped: truth.setdefault(q, 'unmapped')

rows = {}; clusters = collections.defaultdict(list)
for r in csv.DictReader(open(report), delimiter='\t'):
    q = r['read_id'].split()[0]; L = int(r['len']); s = int(r['refined_split'])
    if r['sig_band'] != 'NA':
        clusters[r['sig_band']].append((L // 30, min(s, L - s) // 30))
    if q in truth: rows[q] = r
distinct = {k: len(set(v)) for k, v in clusters.items()}

open(sam + ".truth.tsv","w").writelines(f"{q}\t{t}\n" for q,t in truth.items())
N = len(truth); c = collections.Counter(truth.values())
print(f"sample reads={N}")
for k, v in c.most_common(): print(f"  truth {k:14s} {v:7d} {100*v/N:5.1f}%")
print("\nmdax outcome by truth class (row %):")
for k in ('foldback', 'other_chimera', 'linear_split', 'linear', 'unmapped'):
    o = collections.Counter(rows[q]['decision'] if q in rows else 'not_detected' for q, t in truth.items() if t == k)
    n = sum(o.values()) or 1
    print(f"  {k:14s} n={n:6d} " + "  ".join(f"{d}={100*x/n:.1f}%" for d, x in o.most_common()))
print("\ntrue foldbacks: detection rate by shorter-arm length (aligned bp):")
bins = [(0, 500), (500, 1000), (1000, 2000), (2000, 4000), (4000, 10**9)]
for lo, hi in bins:
    qs = [q for q in info if lo <= info[q][0] < hi]
    det = sum(q in rows for q in qs); cut = sum(q in rows and rows[q]['decision'] == 'artefact' for q in qs)
    if qs: print(f"  arm {lo:>5}-{hi if hi<10**9 else 'inf':<6} n={len(qs):6d} share={100*len(qs)/len(info):5.1f}% detected={100*det/len(qs):5.1f}% cut={100*cut/len(qs):5.1f}%")
print("\ntrue foldbacks: detection by junction distance to nearest read end:")
for lo, hi in [(0, 500), (500, 1000), (1000, 2000), (2000, 10**9)]:
    qs = [q for q in info if lo <= info[q][2] < hi]
    det = sum(q in rows for q in qs)
    if qs: print(f"  dist {lo:>5}-{hi if hi<10**9 else 'inf':<6} n={len(qs):6d} share={100*len(qs)/len(info):5.1f}% detected={100*det/len(qs):5.1f}%")
print("\nduplicate-aware support (distinct (len,split) per cluster) for reads mdax called 'real':")
for t in ('foldback', 'linear'):
    qs = [q for q, tt in truth.items() if tt == t and q in rows and rows[q]['decision'] == 'real']
    h = collections.Counter()
    for q in qs:
        n = int(rows[q]['support_n']); d = distinct.get(rows[q]['sig_band'], 0)
        h[('n<=5' if n <= 5 else 'n6-50' if n <= 50 else 'n>50', 'distinct<=1' if d <= 1 else 'distinct2' if d == 2 else 'distinct>=3')] += 1
    print(f"  truth={t} called real: {len(qs)}")
    for k, v in sorted(h.items()): print(f"     {k[0]:6s} {k[1]:12s} {v:6d} {100*v/max(1,len(qs)):5.1f}%")

print("\ntruth=linear reads: is mdax's cut explained by an unaligned (soft-clipped) tail the aligner could not place?")
for dec in ('artefact', 'real', 'not_detected'):
    qs = [q for q, t in truth.items() if t == 'linear' and (rows[q]['decision'] if q in rows else 'not_detected') == dec]
    big = [q for q in qs if max(clip.get(q, (0, 0))) >= 200]
    near = 0
    for q in big:
        if q in rows:
            s = int(rows[q]['refined_split']); L = rlen[q]; a, b = clip[q]
            # split lies inside or within 300bp of the clipped tail
            if (a >= 200 and s <= a + 300) or (b >= 200 and s >= L - b - 300): near += 1
    print(f"  {dec:13s} n={len(qs):6d} softclip>=200bp: {100*len(big)/max(1,len(qs)):5.1f}%  split at/inside clip: {100*near/max(1,len(qs)):5.1f}%")

print("\nlinear reads mdax cut, with NO unaligned tail (the residual false-positive candidates):")
qs = [q for q, t in truth.items() if t == 'linear' and q in rows and rows[q]['decision'] == 'artefact' and max(clip.get(q, (0, 0))) < 200]
h = collections.Counter(min(int(rows[q]['support_n']), 5) for q in qs)
print(f"  n={len(qs)}  support_n: " + "  ".join(f"{k}{'+' if k==5 else ''}:{v}" for k, v in sorted(h.items())))
def ib(x): return '<0.6' if x < .6 else '0.6-0.9' if x < .9 else '0.9-0.99' if x < .99 else '>=0.99'
h = collections.Counter(ib(float(rows[q]['identity_est'])) for q in qs)
print("  identity_est: " + "  ".join(f"{k}:{v}" for k, v in sorted(h.items())))
h = collections.Counter('arm<500' if int(rows[q]['span_p1']) < 1000 else 'arm<1000' if int(rows[q]['span_p1']) < 2000 else 'arm>=1000' for q in qs)
print("  footprint: " + "  ".join(f"{k}:{v}" for k, v in sorted(h.items())))
print("true foldbacks mdax kept as real: support_n " + str(sorted(collections.Counter(min(int(rows[q]['support_n']), 20) for q, t in truth.items() if t == 'foldback' and q in rows and rows[q]['decision'] == 'real').items())))
print("linear reads mdax kept as real:   support_n " + str(sorted(collections.Counter(min(int(rows[q]['support_n']), 20) for q, t in truth.items() if t == 'linear' and q in rows and rows[q]['decision'] == 'real').items())))
