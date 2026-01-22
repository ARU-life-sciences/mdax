#!/usr/bin/env python3
import random

def rand_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))

def revcomp(s):
    comp = {"A":"T","C":"G","G":"C","T":"A"}
    return "".join(comp[c] for c in reversed(s))

def write_fasta(f, name, seq, width=80):
    f.write(f">{name}\n")
    for i in range(0, len(seq), width):
        f.write(seq[i:i+width] + "\n")

random.seed(42)

with open("test_foldback.fasta", "w") as f:
    # --- Read 1: single artefact foldback (unique junction)
    A = rand_dna(3000)
    seq1 = A + revcomp(A)
    write_fasta(f, "artefact_foldback", seq1)

    # --- Read 2a + 2b: true palindrome (same junction twice)
    B = rand_dna(3000)
    pal = B + revcomp(B)
    write_fasta(f, "true_palindrome_1", pal)
    write_fasta(f, "true_palindrome_2", pal)

    # --- Read 3: normal read
    C = rand_dna(6000)
    write_fasta(f, "normal_read", C)

    # --- Read 4: some differing length read (test difference to above)
    D = rand_dna(2473)
    rev = D + revcomp(D)
    write_fasta(f, "8000bp", rev)

