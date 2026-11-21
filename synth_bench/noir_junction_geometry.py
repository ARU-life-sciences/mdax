#!/usr/bin/env python3
"""Geometry of true foldback junctions from a reference alignment, and why some have no inverted copy.

usage: geom.py TRUTH.sam [mdax_report.tsv] [--dump OUT.tsv]

For each read, consecutive aligned segments (>=200 bp) on the same sequence, opposite strands, <5 kb apart
are a foldback junction. Arm A is the segment before the junction in read order, arm B after it.
  p1 = reference position where A ends (in read order), p2 = where B starts.
  d  = displacement of p2 from p1, measured along A's direction of travel (>0: landed ahead of A's end;
       <0: landed behind it, inside or before A).
  overlap = bp of reference covered by both arms = length of the inverted copy present in the read.
An inverted copy exists only if the arm that has to "come back" is longer than |d|:
  d >= 0: B must be longer than d;   d < 0: A must be longer than |d|.
"""
import sys, re, csv, collections
args = [a for a in sys.argv[1:] if not a.startswith('--')]
dump = sys.argv[sys.argv.index('--dump') + 1] if '--dump' in sys.argv else None
if dump: args.remove(dump)
sam = args[0]; report = args[1] if len(args) > 1 else None
cig = re.compile(r'(\d+)([MIDNSHP=X])')
segs = collections.defaultdict(list); rlen = {}
for line in open(sam):
    if line[0] == '@': continue
    f = line.split('\t', 11); flag = int(f[1])
    if flag & 4 or flag & 256: continue
    ops = [(int(n), o) for n, o in cig.findall(f[5])]
    lead = ops[0][0] if ops[0][1] in 'SH' else 0; trail = ops[-1][0] if ops[-1][1] in 'SH' else 0
    ql = sum(n for n, o in ops if o in 'MI=X'); rl = sum(n for n, o in ops if o in 'MDN=X')
    rev = bool(flag & 16); qs = trail if rev else lead
    segs[f[0]].append((qs, qs + ql, f[2], int(f[3]) - 1, int(f[3]) - 1 + rl, rev)); rlen[f[0]] = lead + ql + trail

dec = {}
if report:
    for r in csv.DictReader(open(report), delimiter='\t'):
        q = r['read_id'].split()[0]
        if dec.get(q) != 'artefact': dec[q] = r['decision']

J = []   # one record per junction
per_read = collections.defaultdict(list)
for q, ss in segs.items():
    ss = sorted(s for s in ss if s[1] - s[0] >= 200)
    for i, (a, b) in enumerate(zip(ss, ss[1:])):
        if a[2] != b[2] or a[5] == b[5]: continue
        gap = max(a[3], b[3]) - min(a[4], b[4])          # <0 when the reference intervals overlap
        if gap >= 5000: continue
        ov = max(0, min(a[4], b[4]) - max(a[3], b[3]))
        if not a[5]: p1, p2 = a[4], b[4]; d = p2 - p1     # A forward, B reverse: B starts at its high end
        else:        p1, p2 = a[3], b[3]; d = p1 - p2     # A reverse, B forward
        la, lb = a[4] - a[3], b[4] - b[3]
        a_first = i == 0; b_last = i + 1 == len(ss) - 1
        a_to_start = a[0]; b_to_end = rlen[q] - b[1]
        rec = dict(q=q, chr=a[2], p1=p1, p2=p2, d=d, ov=ov, la=la, lb=lb, a_first=a_first, b_last=b_last,
                   a_to_start=a_to_start, b_to_end=b_to_end, rl=rlen[q], qj=a[1], qgap=b[0] - a[1])
        J.append(rec); per_read[q].append(rec)

def pct(x, n): return f"{100 * x / max(1, n):5.1f}%"
def q3(v):
    v = sorted(v)
    return "n=0" if not v else f"p10={v[len(v)//10]} p25={v[len(v)//4]} p50={v[len(v)//2]} p75={v[3*len(v)//4]} p90={v[9*len(v)//10]}"

N = len(per_read)
print(f"true-foldback reads: {N}   junctions: {len(J)}")
# read-level: best junction = the one with the longest inverted copy
best = {q: max(v, key=lambda r: r['ov']) for q, v in per_read.items()}
def cls(r):
    if r['ov'] >= 100: return 'A inverted copy >=100 bp'
    if r['ov'] > 0:    return 'B inverted copy 1-99 bp'
    return 'C no inverted copy'
print("\n[1] inverted copy present in the read? (per read, best junction)")
by = collections.defaultdict(list)
for q, r in best.items(): by[cls(r)].append(r)
for k in sorted(by):
    v = by[k]; line = f"  {k:26s} {len(v):6d} {pct(len(v), N)}"
    if report:
        c = collections.Counter(dec.get(r['q'], 'not_detected') for r in v)
        line += "   mdax: " + "  ".join(f"{x}={pct(c[x], len(v))}" for x in ('artefact', 'real', 'not_detected'))
    print(line)

C = by['C no inverted copy']
print(f"\n[2] reads with no inverted copy (n={len(C)}): why the arms do not overlap")
def why(r):
    # the arm that would have to come back over the other
    if r['d'] >= 0: arm, to_end, outer = r['lb'], r['b_to_end'], r['b_last']
    else:           arm, to_end, outer = r['la'], r['a_to_start'], r['a_first']
    side = 'landed ahead of the first arm (d>=0); second arm' if r['d'] >= 0 else 'landed behind the first arm (d<0); first arm'
    lim = 'runs to the read end' if (outer and to_end <= 100) else ('ends at another junction' if not outer else 'ends in unaligned sequence')
    return f"{side} {lim}"
c = collections.Counter(why(r) for r in C)
for k, v in c.most_common(): print(f"  {v:6d} {pct(v, len(C))}  {k}")
print("  |d| (bp):                    ", q3([abs(r['d']) for r in C]))
print("  arm that would have to return:", q3([r['lb'] if r['d'] >= 0 else r['la'] for r in C]))
print("  shortfall |d| - arm (bp):     ", q3([abs(r['d']) - (r['lb'] if r['d'] >= 0 else r['la']) for r in C]))
print("  unaligned read bp between the arms:", q3([r['qgap'] for r in C]))

A = by['A inverted copy >=100 bp']
print(f"\n[3] for comparison, reads WITH an inverted copy (n={len(A)})")
print("  |d| (bp):          ", q3([abs(r['d']) for r in A]))
print("  inverted copy (bp):", q3([r['ov'] for r in A]))
print("  d sign: ahead (d>=0) " + pct(sum(r['d'] >= 0 for r in A), len(A)) + "   behind (d<0) " + pct(sum(r['d'] < 0 for r in A), len(A)))
print("  in no-copy reads:  ahead " + pct(sum(r['d'] >= 0 for r in C), len(C)) + "   behind " + pct(sum(r['d'] < 0 for r in C), len(C)))

print("\n[4] share of junctions with no inverted copy, by |d|")
for lo, hi in [(0, 50), (50, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 5000)]:
    v = [r for r in best.values() if lo <= abs(r['d']) < hi]
    if v: print(f"  |d| {lo:>4}-{hi:<5} n={len(v):6d} ({pct(len(v), N)} of reads)  no copy: {pct(sum(r['ov'] == 0 for r in v), len(v))}")

print("\n[5] do the same junctions recur across reads? (same sequence, p1 and p2 each within 50 bp, same sign of d)")
def recur(v):
    idx = collections.defaultdict(list)
    for r in v: idx[(r['chr'], r['p1'] // 50, r['d'] >= 0)].append(r)
    out = []
    for r in v:
        n = 0
        for k in (-1, 0, 1):
            for o in idx.get((r['chr'], r['p1'] // 50 + k, r['d'] >= 0), ()):
                if o['q'] != r['q'] and abs(o['p1'] - r['p1']) <= 50 and abs(o['p2'] - r['p2']) <= 50: n += 1
        out.append(n)
    return out
allb = list(best.values()); rc = dict(zip([r['q'] for r in allb], recur(allb)))
for k in sorted(by):
    v = by[k]; h = collections.Counter(min(rc[r['q']], 3) for r in v)
    print(f"  {k:26s} other reads at the same junction: 0: {pct(h[0], len(v))}  1: {pct(h[1], len(v))}  2: {pct(h[2], len(v))}  3+: {pct(h[3], len(v))}")

if dump:
    with open(dump, 'w') as o:
        o.write("read\tclass\tchr\tp1\tp2\td\tov\tla\tlb\tqj\tqgap\trl\tdecision\trecur\n")
        for q, r in best.items():
            o.write(f"{q}\t{cls(r)[0]}\t{r['chr']}\t{r['p1']}\t{r['p2']}\t{r['d']}\t{r['ov']}\t{r['la']}\t{r['lb']}\t{r['qj']}\t{r['qgap']}\t{r['rl']}\t{dec.get(q, 'not_detected')}\t{rc[q]}\n")
