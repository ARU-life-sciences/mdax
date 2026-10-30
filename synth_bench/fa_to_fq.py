#!/usr/bin/env python3
"""Convert a FASTA to FASTQ with uniform dummy quality (Q30), for tools
(CADECT) that require FASTQ input. No real quality information exists for
simulated reads, so this is a format shim, not a simulation of error."""
import sys

QCHAR = "?"  # Phred 30


def main():
    if len(sys.argv) != 3:
        sys.stderr.write(f"usage: {sys.argv[0]} in.fa out.fq\n")
        sys.exit(1)
    in_fa, out_fq = sys.argv[1], sys.argv[2]
    with open(in_fa) as fin, open(out_fq, "w") as fout:
        name = None
        seq_chunks = []

        def flush():
            if name is None:
                return
            seq = "".join(seq_chunks)
            fout.write(f"@{name}\n{seq}\n+\n{QCHAR * len(seq)}\n")

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)
        flush()


if __name__ == "__main__":
    main()
