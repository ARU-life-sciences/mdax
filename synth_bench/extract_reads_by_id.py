#!/usr/bin/env python3
"""Stream a (possibly gzipped) FASTA once, writing only records whose ID
(first whitespace-delimited token of the header) is in the given ID list."""
import gzip
import sys


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    ids_path, in_fa, out_fa = sys.argv[1], sys.argv[2], sys.argv[3]
    with open(ids_path) as f:
        wanted = set(line.strip() for line in f if line.strip())
    found = 0
    with opener(in_fa) as fin, open(out_fa, "w") as fout:
        writing = False
        for line in fin:
            if line.startswith(">"):
                rid = line[1:].strip().split()[0]
                writing = rid in wanted
                if writing:
                    found += 1
                    fout.write(line)
            elif writing:
                fout.write(line)
    print(f"wrote {found} / {len(wanted)} requested reads", file=sys.stderr)


if __name__ == "__main__":
    main()
