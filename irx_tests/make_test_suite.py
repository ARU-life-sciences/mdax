#!/usr/bin/env python3
"""
make_irx_truthset.py

Generate a synthetic FASTA truth set for validating `irx`, along with a TSV
describing the embedded ground-truth inverted repeats.

The generator creates a wide range of outcomes, including:

1. Clean high-identity IRs
2. Diverged IRs
3. Asymmetric arms
4. overlap / immediate / spaced / wide classes
5. Window-edge cases
6. Negative controls
7. Micro-IR / low-complexity cases

Each contig contains either:
- exactly one planted inverted repeat structure, or
- no planted IR (negative control)

Output:
- synthetic FASTA
- truth TSV
- manifest TSV (all contigs, including negatives)

Example:
    python make_irx_truthset.py \
        --out-fasta irx_truth.fa \
        --out-truth irx_truth.tsv \
        --out-manifest irx_manifest.tsv \
        --seed 42

Designed for downstream validation of:
- detection
- coordinate recovery
- spacer recovery
- `tir_ident` calibration
- class assignment
- false positives on negatives
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, List, Optional


DNA_ALPHABET = "ACGT"


# -----------------------------
# sequence helpers
# -----------------------------

def random_dna(length: int, rng: random.Random, gc: float = 0.5) -> str:
    """Generate random DNA with approximately the given GC fraction."""
    gc = max(0.0, min(1.0, gc))
    at = (1.0 - gc) / 2.0
    gc_half = gc / 2.0
    alphabet = ["A", "C", "G", "T"]
    probs = [at, gc_half, gc_half, at]

    seq = []
    for _ in range(length):
        r = rng.random()
        cum = 0.0
        for base, p in zip(alphabet, probs):
            cum += p
            if r <= cum:
                seq.append(base)
                break
    return "".join(seq)


def low_complexity_repeat(unit: str, length: int) -> str:
    """Repeat a short unit to the requested length."""
    if not unit:
        raise ValueError("unit must be non-empty")
    reps = math.ceil(length / len(unit))
    return (unit * reps)[:length]


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]


def mutate_sequence(
    seq: str,
    target_identity: float,
    rng: random.Random,
    preserve_len: bool = True,
) -> str:
    """
    Mutate a sequence to roughly the requested identity.

    This introduces substitutions only by default, preserving length.
    That makes the truth cleaner and easier to interpret for `tir_ident`.
    """
    if not 0.0 <= target_identity <= 1.0:
        raise ValueError("target_identity must be in [0,1]")

    if len(seq) == 0:
        return seq

    n_mut = round((1.0 - target_identity) * len(seq))
    n_mut = max(0, min(len(seq), n_mut))

    seq_list = list(seq)
    idxs = list(range(len(seq)))
    rng.shuffle(idxs)

    for i in idxs[:n_mut]:
        old = seq_list[i]
        choices = [b for b in DNA_ALPHABET if b != old]
        seq_list[i] = rng.choice(choices)

    return "".join(seq_list)


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


# -----------------------------
# truth model
# -----------------------------

@dataclass
class TruthRow:
    contig: str
    contig_len: int
    scenario_group: str
    scenario_name: str
    planted: bool

    # planted IR coordinates
    la0: Optional[int] = None
    la1: Optional[int] = None
    ra0: Optional[int] = None
    ra1: Optional[int] = None

    arm_len_left: Optional[int] = None
    arm_len_right: Optional[int] = None
    arm_len_min: Optional[int] = None
    spacer: Optional[int] = None
    true_identity: Optional[float] = None
    ir_class: Optional[str] = None

    # metadata
    gc: Optional[float] = None
    low_complexity: Optional[bool] = None
    edge_case: Optional[str] = None
    notes: Optional[str] = None


@dataclass
class ContigRecord:
    name: str
    seq: str
    truth: TruthRow


# -----------------------------
# classification helpers
# -----------------------------

def classify_ir(spacer: int, immediate_bp: int = 2000, wide_bp: int = 100000) -> str:
    if spacer <= 0:
        return "overlap"
    if spacer <= immediate_bp:
        return "immediate"
    if spacer >= wide_bp:
        return "wide"
    return "spaced"


# -----------------------------
# contig builder
# -----------------------------

def build_contig_with_ir(
    *,
    name: str,
    contig_len: int,
    left_arm_len: int,
    right_arm_len: int,
    spacer_len: int,
    target_identity: float,
    rng: random.Random,
    gc: float = 0.5,
    left_arm_seq: Optional[str] = None,
    low_complexity: bool = False,
    edge_mode: str = "center",
    scenario_group: str,
    scenario_name: str,
    notes: str = "",
) -> ContigRecord:
    """
    Build a contig containing one planted spaced inverted repeat.

    Layout:
        prefix + left_arm + spacer + right_arm + suffix

    where right_arm is reverse-complement related to left_arm (possibly diverged).
    """
    if left_arm_len <= 0 or right_arm_len <= 0:
        raise ValueError("arm lengths must be positive")

    if contig_len <= 0:
        raise ValueError("contig_len must be positive")

    min_payload = left_arm_len + max(0, spacer_len) + right_arm_len
    if min_payload >= contig_len:
        raise ValueError(
            f"contig_len={contig_len} too small for payload={min_payload}"
        )

    # Build the left arm.
    if left_arm_seq is not None:
        if len(left_arm_seq) != left_arm_len:
            raise ValueError("left_arm_seq length does not match left_arm_len")
        left_arm = left_arm_seq
    else:
        if low_complexity:
            unit = rng.choice(["AT", "TA", "GC", "CG", "AAT", "TTA", "GCG", "CGC"])
            left_arm = low_complexity_repeat(unit, left_arm_len)
        else:
            left_arm = random_dna(left_arm_len, rng, gc=gc)

    # Derive the sequence that will become the right arm before RC.
    # Start from left arm, mutate to requested identity, then RC.
    # If right arm length differs, trim or extend cleanly.
    base_for_right = left_arm
    if right_arm_len != left_arm_len:
        if right_arm_len < left_arm_len:
            base_for_right = left_arm[:right_arm_len]
        else:
            extra = random_dna(right_arm_len - left_arm_len, rng, gc=gc)
            base_for_right = left_arm + extra

    mutated = mutate_sequence(base_for_right, target_identity, rng)
    right_arm = revcomp(mutated)

    # Spacer sequence
    if spacer_len < 0:
        raise ValueError("negative spacer not supported in this generator")
    if low_complexity:
        spacer = low_complexity_repeat(rng.choice(["AT", "TA", "A", "T"]), spacer_len)
    else:
        spacer = random_dna(spacer_len, rng, gc=gc)

    payload = left_arm + spacer + right_arm
    slack = contig_len - len(payload)

    if slack < 0:
        raise ValueError("payload exceeds contig length")

    # Edge placement
    # We intentionally allow placements near ends / near likely window boundaries.
    if edge_mode == "center":
        prefix_len = slack // 2
    elif edge_mode == "left":
        prefix_len = min(500, slack)  # near left edge
    elif edge_mode == "right":
        prefix_len = max(0, slack - 500)  # near right edge
    elif edge_mode == "quarter":
        prefix_len = slack // 4
    elif edge_mode == "three_quarter":
        prefix_len = (3 * slack) // 4
    elif edge_mode == "random":
        prefix_len = rng.randint(0, slack)
    else:
        raise ValueError(f"unknown edge_mode: {edge_mode}")

    suffix_len = slack - prefix_len
    prefix = random_dna(prefix_len, rng, gc=gc)
    suffix = random_dna(suffix_len, rng, gc=gc)

    seq = prefix + payload + suffix

    la0 = prefix_len
    la1 = la0 + left_arm_len
    ra0 = la1 + spacer_len
    ra1 = ra0 + right_arm_len

    assert len(seq) == contig_len

    truth = TruthRow(
        contig=name,
        contig_len=contig_len,
        scenario_group=scenario_group,
        scenario_name=scenario_name,
        planted=True,
        la0=la0,
        la1=la1,
        ra0=ra0,
        ra1=ra1,
        arm_len_left=left_arm_len,
        arm_len_right=right_arm_len,
        arm_len_min=min(left_arm_len, right_arm_len),
        spacer=spacer_len,
        true_identity=target_identity,
        ir_class=classify_ir(spacer_len),
        gc=gc,
        low_complexity=low_complexity,
        edge_case=edge_mode,
        notes=notes,
    )
    return ContigRecord(name=name, seq=seq, truth=truth)


def build_negative_contig(
    *,
    name: str,
    contig_len: int,
    rng: random.Random,
    gc: float,
    mode: str,
) -> ContigRecord:
    """
    Build negative controls:
    - random
    - tandem_repeat
    - low_complexity
    - direct_repeat
    """
    if mode == "random":
        seq = random_dna(contig_len, rng, gc=gc)
        notes = "random DNA only"
    elif mode == "low_complexity":
        seq = low_complexity_repeat(rng.choice(["AT", "A", "T", "TAA", "AAT"]), contig_len)
        notes = "low-complexity negative"
    elif mode == "tandem_repeat":
        unit = rng.choice(["ATGC", "AATG", "GATA", "TATC", "AC"])
        seq = low_complexity_repeat(unit, contig_len)
        notes = "tandem repeat negative"
    elif mode == "direct_repeat":
        # Plant direct repeats, not inverted repeats.
        flank = contig_len // 4
        arm_len = min(1000, contig_len // 10)
        spacer_len = contig_len - (2 * flank + 2 * arm_len)
        spacer_len = max(100, spacer_len)
        left = random_dna(arm_len, rng, gc=gc)
        seq = (
            random_dna(flank, rng, gc=gc)
            + left
            + random_dna(spacer_len, rng, gc=gc)
            + left  # direct repeat, not RC
            + random_dna(contig_len - flank - arm_len - spacer_len - arm_len, rng, gc=gc)
        )
        seq = seq[:contig_len]
        notes = "direct repeat negative"
    else:
        raise ValueError(f"unknown negative mode: {mode}")

    truth = TruthRow(
        contig=name,
        contig_len=contig_len,
        scenario_group="negative",
        scenario_name=mode,
        planted=False,
        gc=gc,
        low_complexity=(mode in {"low_complexity", "tandem_repeat"}),
        notes=notes,
    )
    return ContigRecord(name=name, seq=seq, truth=truth)


# -----------------------------
# scenario generator
# -----------------------------

def generate_truthset(rng: random.Random) -> List[ContigRecord]:
    """
    Generate a broad panel of positive and negative synthetic contigs.

    The truth set is intentionally diverse rather than minimal.
    """
    records: List[ContigRecord] = []
    c = 1

    def next_name(prefix: str) -> str:
        nonlocal c
        name = f"{prefix}_{c:03d}"
        c += 1
        return name

    # ------------------------------------------------------------------
    # 1. Clean high-identity IRs
    # ------------------------------------------------------------------
    for arm_len, spacer_len, contig_len in [
        (100, 100, 10000),
        (500, 1000, 20000),
        (2000, 10000, 50000),
        (10000, 50000, 150000),
    ]:
        records.append(build_contig_with_ir(
            name=next_name("clean"),
            contig_len=contig_len,
            left_arm_len=arm_len,
            right_arm_len=arm_len,
            spacer_len=spacer_len,
            target_identity=0.99,
            rng=rng,
            gc=0.5,
            edge_mode="center",
            scenario_group="clean_high_identity",
            scenario_name=f"clean_arm{arm_len}_spacer{spacer_len}",
            notes="clean high-identity IR",
        ))

    # ------------------------------------------------------------------
    # 2. Diverged IRs (legacy simple series)
    # ------------------------------------------------------------------
    for identity in [0.95, 0.85, 0.70, 0.50]:
        records.append(build_contig_with_ir(
            name=next_name("diverged"),
            contig_len=30000,
            left_arm_len=1000,
            right_arm_len=1000,
            spacer_len=3000,
            target_identity=identity,
            rng=rng,
            gc=0.5,
            edge_mode="center",
            scenario_group="diverged",
            scenario_name=f"diverged_ident{identity:.2f}",
            notes="diverged IR series",
        ))

    # ------------------------------------------------------------------
    # 3. Length × identity × spacer matrix
    #
    # This is the key extension for benchmarking sensitivity to longer,
    # more diverged IRs. It asks:
    # - for a given arm length, how much divergence can irx tolerate?
    # - does spacer size affect recovery?
    # ------------------------------------------------------------------
    matrix_arm_lengths = [1000, 2000, 5000, 10000, 50000, 100000]
    matrix_identities = [0.95, 0.85, 0.70, 0.55, 0.40]
    matrix_spacers = [1000, 10000, 100000]

    for arm_len in matrix_arm_lengths:
        for identity in matrix_identities:
            for spacer_len in matrix_spacers:
                # Give enough slack so the structure is well embedded.
                # Keep contigs moderate in size while accommodating wide spacers.
                contig_len = max(
                    50000,
                    2 * arm_len + spacer_len + 20000
                )

                records.append(build_contig_with_ir(
                    name=next_name("matrix"),
                    contig_len=contig_len,
                    left_arm_len=arm_len,
                    right_arm_len=arm_len,
                    spacer_len=spacer_len,
                    target_identity=identity,
                    rng=rng,
                    gc=0.5,
                    edge_mode="center",
                    scenario_group="length_identity_matrix",
                    scenario_name=(
                        f"matrix_arm{arm_len}"
                        f"_ident{identity:.2f}"
                        f"_spacer{spacer_len}"
                    ),
                    notes="joint arm-length × identity × spacer benchmark",
                ))

    # ------------------------------------------------------------------
    # 4. Asymmetric arms
    # ------------------------------------------------------------------
    for left_len, right_len in [(1000, 900), (1000, 700), (3000, 2500)]:
        records.append(build_contig_with_ir(
            name=next_name("asym"),
            contig_len=40000,
            left_arm_len=left_len,
            right_arm_len=right_len,
            spacer_len=5000,
            target_identity=0.95,
            rng=rng,
            gc=0.5,
            edge_mode="center",
            scenario_group="asymmetric",
            scenario_name=f"asym_L{left_len}_R{right_len}",
            notes="asymmetric arm lengths",
        ))

    # ------------------------------------------------------------------
    # 5. overlap / immediate / spaced / wide
    # overlap is awkward in a simple planted model unless we literally overlap.
    # Here we cover immediate / spaced / wide cleanly.
    # ------------------------------------------------------------------
    for spacer_len, label in [
        (50, "immediate_small"),
        (500, "immediate"),
        (5000, "spaced"),
        (120000, "wide"),
    ]:
        contig_len = max(30000, 2 * 2000 + spacer_len + 10000)
        records.append(build_contig_with_ir(
            name=next_name("class"),
            contig_len=contig_len,
            left_arm_len=2000,
            right_arm_len=2000,
            spacer_len=spacer_len,
            target_identity=0.98,
            rng=rng,
            gc=0.45,
            edge_mode="center",
            scenario_group="class_spacing",
            scenario_name=label,
            notes="spacing-class test",
        ))

    # ------------------------------------------------------------------
    # 6. Window-edge / boundary cases
    # Default irx window is 250 kb; these tests place IRs near edges.
    # ------------------------------------------------------------------
    for edge_mode in ["left", "right", "quarter", "three_quarter"]:
        records.append(build_contig_with_ir(
            name=next_name("edge"),
            contig_len=250000,
            left_arm_len=3000,
            right_arm_len=3000,
            spacer_len=10000,
            target_identity=0.97,
            rng=rng,
            gc=0.5,
            edge_mode=edge_mode,
            scenario_group="edge_cases",
            scenario_name=f"edge_{edge_mode}",
            notes="IR placed near contig/window edge",
        ))

    # ------------------------------------------------------------------
    # 7. Negative controls
    # ------------------------------------------------------------------
    for mode, contig_len, gc in [
        ("random", 30000, 0.5),
        ("random", 80000, 0.35),
        ("low_complexity", 20000, 0.5),
        ("tandem_repeat", 25000, 0.5),
        ("direct_repeat", 40000, 0.5),
    ]:
        records.append(build_negative_contig(
            name=next_name("neg"),
            contig_len=contig_len,
            rng=rng,
            gc=gc,
            mode=mode,
        ))

    # ------------------------------------------------------------------
    # 8. Micro-IR / low-complexity cases
    # These are especially useful because `irx` may call them strongly.
    # ------------------------------------------------------------------
    for arm_len, spacer_len, low_complexity, identity in [
        (40, 3000, False, 0.99),
        (60, 10000, False, 0.99),
        (80, 200000, False, 0.98),
        (50, 50000, True, 0.99),
        (80, 100000, True, 0.95),
    ]:
        contig_len = max(25000, 2 * arm_len + spacer_len + 10000)
        records.append(build_contig_with_ir(
            name=next_name("micro"),
            contig_len=contig_len,
            left_arm_len=arm_len,
            right_arm_len=arm_len,
            spacer_len=spacer_len,
            target_identity=identity,
            rng=rng,
            gc=0.4,
            low_complexity=low_complexity,
            edge_mode="center",
            scenario_group="micro_ir",
            scenario_name=f"micro_arm{arm_len}_spacer{spacer_len}",
            notes="micro-IR / low-complexity test",
        ))

    # ------------------------------------------------------------------
    # 9. GC-biased cases
    # ------------------------------------------------------------------
    for gc in [0.2, 0.5, 0.8]:
        records.append(build_contig_with_ir(
            name=next_name("gc"),
            contig_len=50000,
            left_arm_len=1500,
            right_arm_len=1500,
            spacer_len=7000,
            target_identity=0.96,
            rng=rng,
            gc=gc,
            edge_mode="center",
            scenario_group="gc_bias",
            scenario_name=f"gc_{gc:.2f}",
            notes="GC composition stress test",
        ))

    # ------------------------------------------------------------------
    # 10. Reused-arm family cases
    # Useful for seeing whether repeated families create extra calls.
    # ------------------------------------------------------------------
    family_arm = random_dna(1200, rng, gc=0.45)
    for spacer_len in [500, 5000, 20000]:
        records.append(build_contig_with_ir(
            name=next_name("family"),
            contig_len=50000,
            left_arm_len=1200,
            right_arm_len=1200,
            spacer_len=spacer_len,
            target_identity=0.97,
            rng=rng,
            gc=0.45,
            left_arm_seq=family_arm,
            edge_mode="center",
            scenario_group="reused_family",
            scenario_name=f"family_spacer{spacer_len}",
            notes="shared-arm family test",
        ))

    return records

# -----------------------------
# writing outputs
# -----------------------------

def write_fasta(records: Iterable[ContigRecord], path: Path, wrap: int = 80) -> None:
    with path.open("w") as fh:
        for rec in records:
            fh.write(f">{rec.name}\n")
            fh.write(wrap_fasta(rec.seq, width=wrap))
            fh.write("\n")


def write_truth_tsv(records: Iterable[ContigRecord], path: Path) -> None:
    rows = [r.truth for r in records if r.truth.planted]
    if not rows:
        raise ValueError("no planted truth rows to write")

    fieldnames = list(asdict(rows[0]).keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))


def write_manifest_tsv(records: Iterable[ContigRecord], path: Path) -> None:
    rows = [r.truth for r in records]
    fieldnames = list(asdict(rows[0]).keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(asdict(row))


# -----------------------------
# CLI
# -----------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a synthetic FASTA truth set for validating irx."
    )
    p.add_argument("--out-fasta", required=True, help="Output FASTA path")
    p.add_argument("--out-truth", required=True, help="Output TSV of planted IR truth rows")
    p.add_argument("--out-manifest", required=True, help="Output TSV of all contigs")
    p.add_argument("--wrap", type=int, default=80, help="FASTA line width [default: 80]")
    p.add_argument("--seed", type=int, default=42, help="Random seed [default: 42]")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    rng = random.Random(args.seed)

    records = generate_truthset(rng)

    out_fasta = Path(args.out_fasta)
    out_truth = Path(args.out_truth)
    out_manifest = Path(args.out_manifest)

    write_fasta(records, out_fasta, wrap=args.wrap)
    write_truth_tsv(records, out_truth)
    write_manifest_tsv(records, out_manifest)

    n_total = len(records)
    n_pos = sum(1 for r in records if r.truth.planted)
    n_neg = n_total - n_pos

    print(f"[ok] wrote FASTA:    {out_fasta}")
    print(f"[ok] wrote truth:    {out_truth}")
    print(f"[ok] wrote manifest: {out_manifest}")
    print(f"[ok] contigs: total={n_total} planted={n_pos} negative={n_neg}")


if __name__ == "__main__":
    main()
