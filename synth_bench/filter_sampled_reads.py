#!/usr/bin/env python3
"""Stream mdax's corrected FASTA (from stdin, typically a FIFO so the full
output never touches disk) and keep only the reads whose ids are listed in
IDS. mdax appends `_chopped_at_<pos>` to the header of cut reads; that suffix
is stripped so ids match the raw reads again.

usage: filter_sampled_reads.py IDS < corrected.fa > sampled_corrected.fa
"""
import re
import sys

ids = {l.strip() for l in open(sys.argv[1]) if l.strip()}
suffix = re.compile(r"_chopped_at_\d+$")
out = sys.stdout
keep = False
for line in sys.stdin:
    if line.startswith(">"):
        tok = suffix.sub("", line[1:].rstrip("\n")).split()
        name = tok[0] if tok else ""
        keep = name in ids
        if keep:
            out.write(">" + name + "\n")
    elif keep:
        out.write(line)
