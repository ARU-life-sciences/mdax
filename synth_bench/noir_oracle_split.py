#!/usr/bin/env python3
"""Oracle cut: split corrected reads at reference-detected foldback junctions.
usage: split.py PAF READS.fa OUTPREFIX   -> OUTPREFIX.{nocopy,all,allpieces}.fa + stats on stderr
  nocopy     cut only at junctions whose arms share no reference sequence; keep the longest piece
  all        cut at every foldback junction; keep the longest piece
  allpieces  cut at every foldback junction; keep every piece >= 1 kb
"""
import sys, collections
paf, reads, pre = sys.argv[1:4]
segs = collections.defaultdict(list)
for l in open(paf):
    f = l.split('\t', 12)
    if int(f[3]) - int(f[2]) >= 200:
        segs[f[0]].append((int(f[2]), int(f[3]), f[5], int(f[7]), int(f[8]), f[4] == '-'))
cuts = {}
for q, ss in segs.items():
    ss.sort(); c = []
    for a, b in zip(ss, ss[1:]):
        if a[2] != b[2] or a[5] == b[5]: continue
        if max(a[3], b[3]) - min(a[4], b[4]) >= 5000: continue
        ov = max(0, min(a[4], b[4]) - max(a[3], b[3]))
        c.append(((a[1] + b[0]) // 2, ov == 0))
    if c: cuts[q] = c
out = {m: open(f"{pre}.{m}.fa", 'w') for m in ('nocopy', 'all', 'allpieces')}
st = collections.Counter()
def pieces(s, pts):
    pts = [0] + sorted(p for p in pts if 0 < p < len(s)) + [len(s)]
    return [s[a:b] for a, b in zip(pts, pts[1:])]
def emit(name, seq):
    q = name.split()[0]; c = cuts.get(q)
    st['reads'] += 1; st['bases'] += len(seq)
    if not c:
        for o in out.values(): o.write(f">{q}\n{seq}\n")
        for m in out: st[m + '_bases'] += len(seq)
        return
    st['reads_with_junction'] += 1
    nc = [p for p, n in c if n]; al = [p for p, _ in c]
    if nc: st['reads_with_nocopy_junction'] += 1
    p1 = max(pieces(seq, nc), key=len); out['nocopy'].write(f">{q}\n{p1}\n"); st['nocopy_bases'] += len(p1)
    p2 = max(pieces(seq, al), key=len); out['all'].write(f">{q}\n{p2}\n"); st['all_bases'] += len(p2)
    for i, p in enumerate(x for x in pieces(seq, al) if len(x) >= 1000):
        out['allpieces'].write(f">{q}/p{i}\n{p}\n"); st['allpieces_bases'] += len(p)
name = None; buf = []
for l in open(reads):
    if l[0] == '>':
        if name: emit(name, ''.join(buf))
        name = l[1:].strip(); buf = []
    else: buf.append(l.strip())
if name: emit(name, ''.join(buf))
for k, v in sorted(st.items()): print(k, v, file=sys.stderr)
