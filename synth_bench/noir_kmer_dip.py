#!/usr/bin/env python3
"""Library k-mer count profile around true junctions vs elsewhere.
usage: kdip.py GEOM.tsv TRUTH.tsv READS.fa.gz COUNTS(stdin: 'kmer count' per line from jellyfish query -s, same read order as SEL.fa)
Step 1 (select):  kdip.py select GEOM.tsv TRUTH.tsv READS.fa.gz SEL.fa   -> writes SEL.fa + SEL.fa.meta
Step 2 (score):   jellyfish query lib.jf -s SEL.fa | kdip.py score SEL.fa.meta
A 'dip' is a run of consecutive 25-mers with library count <= 1 (jellyfish built with -L 2, so they read 0),
with the median count of the 50 k-mers on each side >= 5."""
import sys, gzip, random, collections, statistics
K = 25
if sys.argv[1] == 'select':
    geom, truth, reads, sel = sys.argv[2:6]; random.seed(7)
    g = {}
    for l in list(open(geom))[1:]:
        f = l.rstrip('\n').split('\t'); g[f[0]] = (f[1], int(f[9]), f[12])   # class, qj, decision
    lin = [l.split()[0] for l in open(truth) if l.split()[1] == 'linear']
    pick = {}
    for cl in 'AC':
        ids = [q for q, v in g.items() if v[0] == cl]; random.shuffle(ids)
        for q in ids[:2500]: pick[q] = (cl, g[q][1], g[q][2])
    random.shuffle(lin)
    for q in lin[:2500]: pick[q] = ('L', -1, 'NA')
    with open(sel, 'w') as o, open(sel + '.meta', 'w') as m:
        name = None
        def flush(name, seq):
            q = name.split()[0]
            if q in pick and 'N' not in seq and len(seq) > 600:
                o.write(f">{q}\n{seq}\n"); m.write(f"{q}\t{len(seq)}\t" + "\t".join(map(str, pick[q])) + "\n")
        buf = []
        for l in gzip.open(reads, 'rt'):
            if l[0] == '>':
                if name: flush(name, ''.join(buf).upper())
                name = l[1:].strip(); buf = []
            else: buf.append(l.strip())
        if name: flush(name, ''.join(buf).upper())
    sys.exit()

meta = [l.rstrip('\n').split('\t') for l in open(sys.argv[2])]
def runs(c):
    out = []; i = 0; n = len(c)
    while i < n:
        if c[i] <= 1:
            j = i
            while j < n and c[j] <= 1: j += 1
            L = c[max(0, i - 50):i]; R = c[j:j + 50]
            if len(L) >= 20 and len(R) >= 20 and statistics.median(L) >= 5 and statistics.median(R) >= 5:
                out.append((i, j - i, statistics.median(L), statistics.median(R)))
            i = j
        else: i += 1
    return out
S = collections.defaultdict(list); junc_len = []; other_len = collections.Counter(); cov = collections.defaultdict(list)
perkb = collections.defaultdict(list); found = collections.Counter(); tot = collections.Counter()
jflank = []
for q, L, cl, qj, dec in meta:
    L = int(L); qj = int(qj); n = L - K + 1
    c = [int(sys.stdin.readline().split()[1]) for _ in range(n)]
    cov[cl].append(statistics.median(c))
    rr = runs(c); tot[cl] += 1
    hit = None
    if qj >= 0:
        # k-mers spanning read position qj start in [qj-K+1, qj-1]; allow +-15 bp for alignment slack
        for r in rr:
            if r[0] <= qj + 15 and r[0] + r[1] >= qj - K - 15: hit = r; break
        if hit: found[cl] += 1; junc_len.append((cl, hit[1])); jflank.append((cl, hit[2], hit[3]))
        # is the junction itself a low region at all (ignoring flank requirement)?
        w = c[max(0, qj - K + 1):max(1, qj)]
        S[cl + '_minspan'].append(min(w) if w else -1)
    others = [r for r in rr if r is not hit]
    perkb[cl].append(1000 * len(others) / n)
    for r in others: other_len[min(r[1], 60)] += 1
print("class: A = true foldback with inverted copy, C = true foldback with no inverted copy, L = linear")
for cl in 'ACL':
    print(f"  {cl}: reads={tot[cl]}  median library count of the read's 25-mers: p25/p50/p75 = " + "/".join(str(int(x)) for x in statistics.quantiles(cov[cl], n=4)))
print("\n[1] is there a unique-k-mer dip (with solid flanks) at the true junction?")
for cl in 'AC':
    ms = S[cl + '_minspan']
    print(f"  {cl}: dip at junction {100*found[cl]/tot[cl]:.1f}%   (junction-spanning 25-mers all >1 in library: {100*sum(x>1 for x in ms)/len(ms):.1f}%)")
print("\n[2] dips NOT at a known junction, per kb of read (sequencing errors, other chimeras, low coverage)")
for cl in 'ACL':
    v = perkb[cl]; print(f"  {cl}: mean {statistics.mean(v):.2f}/kb   reads with none: {100*sum(x==0 for x in v)/len(v):.1f}%   median per 10 kb read: {10*statistics.median(v):.1f}")
print("\n[3] run length of the dip (number of consecutive unique 25-mers)")
def hist(v):
    h = collections.Counter('<15' if x < 15 else '15-22' if x < 23 else '23-25' if x <= 25 else '26-40' if x <= 40 else '>40' for x in v)
    return "  ".join(f"{k}: {100*h[k]/max(1,len(v)):.1f}%" for k in ('<15', '15-22', '23-25', '26-40', '>40'))
for cl in 'AC': print(f"  junction dips {cl}: " + hist([x for c2, x in junc_len if c2 == cl]))
allo = [k for k, n in other_len.items() for _ in range(n)]
print("  other dips:      " + hist(allo))
