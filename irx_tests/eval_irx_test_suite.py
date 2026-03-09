#!/usr/bin/env python3
"""
eval_irx_truthset.py

Compare `irx` calls against a synthetic truth set produced by
`make_irx_truthset.py`.

This script performs:

1. Truth-to-call matching on the same contig
2. Overlap-based detection scoring
3. Per-scenario-group recovery summaries
4. Per-record TP / FN table
5. Optional call-centric FP table
6. Error summaries for recovered calls:
   - arm length error
   - spacer error
   - tir_ident vs true_identity

Definitions
-----------
A truth record is counted as "recovered" if there is at least one `irx` call on
the same contig whose emitted interval overlaps the truth interval by at least:

- `--min-recip-overlap` on both truth and call interval lengths, OR
- `--min-overlap-bp` absolute bp overlap

By default:
- reciprocal overlap >= 0.25 OR overlap >= 100 bp

This is intentionally permissive for first-pass validation, because `irx`
intervals may be conservative or expanded relative to truth.

Outputs
-------
- summary TSV by scenario group
- truth-vs-call TSV (one row per planted truth record)
- false-positive TSV (optional; one row per unmatched call)
- optional stdout summary

Example
-------
python eval_irx_truthset.py \
  --truth irx_truth.tsv \
  --calls irx_calls.tsv \
  --out-summary irx_eval.summary.tsv \
  --out-truth irx_eval.truth.tsv \
  --out-fp irx_eval.fp.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# -----------------------------
# models
# -----------------------------

@dataclass
class TruthRecord:
    contig: str
    contig_len: int
    scenario_group: str
    scenario_name: str
    planted: bool
    la0: Optional[int]
    la1: Optional[int]
    ra0: Optional[int]
    ra1: Optional[int]
    arm_len_left: Optional[int]
    arm_len_right: Optional[int]
    arm_len_min: Optional[int]
    spacer: Optional[int]
    true_identity: Optional[float]
    ir_class: Optional[str]
    gc: Optional[float]
    low_complexity: Optional[bool]
    edge_case: Optional[str]
    notes: Optional[str]

    @property
    def start(self) -> Optional[int]:
        if self.la0 is None or self.ra0 is None:
            return None
        return min(self.la0, self.ra0)

    @property
    def end(self) -> Optional[int]:
        if self.la1 is None or self.ra1 is None:
            return None
        return max(self.la1, self.ra1)

    @property
    def interval_len(self) -> Optional[int]:
        if self.start is None or self.end is None:
            return None
        return self.end - self.start


@dataclass
class CallRecord:
    contig: str
    start: int
    end: int
    name: str
    score: int
    strand: str
    break_pos: int
    identity_est: float
    tir_ident: Optional[float]
    matches: int
    span: int
    la0: int
    la1: int
    ra0: int
    ra1: int
    contig_len: int
    win_start: int
    win_end: int
    kept_pts: int
    bin: int
    arm_len: Optional[int] = None
    spacer: Optional[int] = None
    ir_class: Optional[str] = None

    @property
    def interval_len(self) -> int:
        return self.end - self.start


@dataclass
class TruthMatchRow:
    contig: str
    scenario_group: str
    scenario_name: str
    planted: bool

    truth_start: Optional[int]
    truth_end: Optional[int]
    truth_interval_len: Optional[int]
    truth_arm_len_min: Optional[int]
    truth_spacer: Optional[int]
    truth_identity: Optional[float]
    truth_ir_class: Optional[str]

    detected: bool
    matched_call_index: Optional[int]

    overlap_bp: int = 0
    recip_truth: float = 0.0
    recip_call: float = 0.0

    call_start: Optional[int] = None
    call_end: Optional[int] = None
    call_interval_len: Optional[int] = None
    call_arm_len: Optional[int] = None
    call_spacer: Optional[int] = None
    call_tir_ident: Optional[float] = None
    call_identity_est: Optional[float] = None
    call_ir_class: Optional[str] = None
    call_kept_pts: Optional[int] = None
    call_matches: Optional[int] = None

    arm_len_error: Optional[int] = None
    spacer_error: Optional[int] = None
    tir_ident_error: Optional[float] = None


@dataclass
class SummaryRow:
    scenario_group: str
    n_truth: int
    n_detected: int
    recovery_rate: float

    median_overlap_bp: Optional[float]
    median_recip_truth: Optional[float]
    median_recip_call: Optional[float]

    median_arm_len_error: Optional[float]
    median_spacer_error: Optional[float]
    median_tir_ident_error: Optional[float]


@dataclass
class FPRow:
    call_index: int
    contig: str
    start: int
    end: int
    interval_len: int
    arm_len: Optional[int]
    spacer: Optional[int]
    tir_ident: Optional[float]
    identity_est: float
    kept_pts: int
    matches: int
    ir_class: Optional[str]


# -----------------------------
# parsing
# -----------------------------

def parse_int(value: str) -> Optional[int]:
    if value is None or value == "" or value == "None":
        return None
    return int(value)


def parse_float(value: str) -> Optional[float]:
    if value is None or value == "" or value == "None":
        return None
    return float(value)


def parse_bool(value: str) -> Optional[bool]:
    if value is None or value == "":
        return None
    if value.lower() in {"true", "t", "1"}:
        return True
    if value.lower() in {"false", "f", "0"}:
        return False
    raise ValueError(f"cannot parse bool: {value}")


def read_truth(path: Path) -> List[TruthRecord]:
    out: List[TruthRecord] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(
                TruthRecord(
                    contig=row["contig"],
                    contig_len=int(row["contig_len"]),
                    scenario_group=row["scenario_group"],
                    scenario_name=row["scenario_name"],
                    planted=parse_bool(row["planted"]) or False,
                    la0=parse_int(row["la0"]),
                    la1=parse_int(row["la1"]),
                    ra0=parse_int(row["ra0"]),
                    ra1=parse_int(row["ra1"]),
                    arm_len_left=parse_int(row["arm_len_left"]),
                    arm_len_right=parse_int(row["arm_len_right"]),
                    arm_len_min=parse_int(row["arm_len_min"]),
                    spacer=parse_int(row["spacer"]),
                    true_identity=parse_float(row["true_identity"]),
                    ir_class=row.get("ir_class"),
                    gc=parse_float(row["gc"]),
                    low_complexity=parse_bool(row["low_complexity"]),
                    edge_case=row.get("edge_case"),
                    notes=row.get("notes"),
                )
            )
    return out


def read_calls(path: Path) -> List[CallRecord]:
    out: List[CallRecord] = []
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(
                CallRecord(
                    contig=row["#contig"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    name=row["name"],
                    score=int(row["score"]),
                    strand=row["strand"],
                    break_pos=int(row["break_pos"]),
                    identity_est=float(row["identity_est"]),
                    tir_ident=parse_float(row["tir_ident"]),
                    matches=int(row["matches"]),
                    span=int(row["span"]),
                    la0=int(row["la0"]),
                    la1=int(row["la1"]),
                    ra0=int(row["ra0"]),
                    ra1=int(row["ra1"]),
                    contig_len=int(row["contig_len"]),
                    win_start=int(row["win_start"]),
                    win_end=int(row["win_end"]),
                    kept_pts=int(row["kept_pts"]),
                    bin=int(row["bin"]),
                    arm_len=parse_int(row.get("arm_len", "")),
                    spacer=parse_int(row.get("spacer", "")),
                    ir_class=row.get("ir_class"),
                )
            )
    return out


# -----------------------------
# overlap / matching
# -----------------------------

def interval_overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def reciprocal_overlap(
    truth_start: int,
    truth_end: int,
    call_start: int,
    call_end: int,
) -> Tuple[int, float, float]:
    ov = interval_overlap(truth_start, truth_end, call_start, call_end)
    truth_len = max(1, truth_end - truth_start)
    call_len = max(1, call_end - call_start)
    return ov, ov / truth_len, ov / call_len


def match_passes(
    ov: int,
    recip_truth: float,
    recip_call: float,
    *,
    min_overlap_bp: int,
    min_recip_overlap: float,
) -> bool:
    return (
        ov >= min_overlap_bp
        or (recip_truth >= min_recip_overlap and recip_call >= min_recip_overlap)
    )


def choose_best_call(
    truth: TruthRecord,
    calls_on_contig: List[Tuple[int, CallRecord]],
    *,
    min_overlap_bp: int,
    min_recip_overlap: float,
) -> Tuple[Optional[int], Optional[CallRecord], int, float, float]:
    if truth.start is None or truth.end is None:
        return None, None, 0, 0.0, 0.0

    best = None
    best_key = None

    for idx, call in calls_on_contig:
        ov, rt, rc = reciprocal_overlap(truth.start, truth.end, call.start, call.end)
        if not match_passes(
            ov, rt, rc,
            min_overlap_bp=min_overlap_bp,
            min_recip_overlap=min_recip_overlap,
        ):
            continue

        # rank by overlap, then reciprocal truth overlap, then tir_ident, then shorter interval
        key = (
            ov,
            rt,
            call.tir_ident if call.tir_ident is not None and math.isfinite(call.tir_ident) else -1.0,
            -call.interval_len,
        )

        if best is None or key > best_key:
            best = (idx, call, ov, rt, rc)
            best_key = key

    if best is None:
        return None, None, 0, 0.0, 0.0
    return best


# -----------------------------
# summaries
# -----------------------------

def median_or_none(values: Iterable[Optional[float]]) -> Optional[float]:
    vals = [v for v in values if v is not None and math.isfinite(v)]
    if not vals:
        return None
    vals.sort()
    n = len(vals)
    mid = n // 2
    if n % 2 == 1:
        return vals[mid]
    return 0.5 * (vals[mid - 1] + vals[mid])


# -----------------------------
# main evaluation
# -----------------------------

def evaluate(
    truth_records: List[TruthRecord],
    call_records: List[CallRecord],
    *,
    min_overlap_bp: int,
    min_recip_overlap: float,
) -> Tuple[List[TruthMatchRow], List[SummaryRow], List[FPRow]]:
    calls_by_contig: Dict[str, List[Tuple[int, CallRecord]]] = defaultdict(list)
    for i, c in enumerate(call_records):
        calls_by_contig[c.contig].append((i, c))

    matched_call_indices = set()
    truth_rows: List[TruthMatchRow] = []

    planted_truth = [t for t in truth_records if t.planted]

    for truth in planted_truth:
        idx, call, ov, rt, rc = choose_best_call(
            truth,
            calls_by_contig.get(truth.contig, []),
            min_overlap_bp=min_overlap_bp,
            min_recip_overlap=min_recip_overlap,
        )

        detected = call is not None
        if idx is not None:
            matched_call_indices.add(idx)

        row = TruthMatchRow(
            contig=truth.contig,
            scenario_group=truth.scenario_group,
            scenario_name=truth.scenario_name,
            planted=truth.planted,
            truth_start=truth.start,
            truth_end=truth.end,
            truth_interval_len=truth.interval_len,
            truth_arm_len_min=truth.arm_len_min,
            truth_spacer=truth.spacer,
            truth_identity=truth.true_identity,
            truth_ir_class=truth.ir_class,
            detected=detected,
            matched_call_index=idx,
            overlap_bp=ov,
            recip_truth=rt,
            recip_call=rc,
        )

        if call is not None:
            row.call_start = call.start
            row.call_end = call.end
            row.call_interval_len = call.interval_len
            row.call_arm_len = call.arm_len
            row.call_spacer = call.spacer
            row.call_tir_ident = call.tir_ident
            row.call_identity_est = call.identity_est
            row.call_ir_class = call.ir_class
            row.call_kept_pts = call.kept_pts
            row.call_matches = call.matches

            if truth.arm_len_min is not None and call.arm_len is not None:
                row.arm_len_error = call.arm_len - truth.arm_len_min
            if truth.spacer is not None and call.spacer is not None:
                row.spacer_error = call.spacer - truth.spacer
            if truth.true_identity is not None and call.tir_ident is not None:
                row.tir_ident_error = call.tir_ident - truth.true_identity

        truth_rows.append(row)

    # scenario summaries
    grouped: Dict[str, List[TruthMatchRow]] = defaultdict(list)
    for row in truth_rows:
        grouped[row.scenario_group].append(row)

    summary_rows: List[SummaryRow] = []
    for scenario_group, rows in sorted(grouped.items()):
        n_truth = len(rows)
        n_detected = sum(r.detected for r in rows)
        summary_rows.append(
            SummaryRow(
                scenario_group=scenario_group,
                n_truth=n_truth,
                n_detected=n_detected,
                recovery_rate=(n_detected / n_truth if n_truth else float("nan")),
                median_overlap_bp=median_or_none(r.overlap_bp for r in rows if r.detected),
                median_recip_truth=median_or_none(r.recip_truth for r in rows if r.detected),
                median_recip_call=median_or_none(r.recip_call for r in rows if r.detected),
                median_arm_len_error=median_or_none(r.arm_len_error for r in rows if r.detected),
                median_spacer_error=median_or_none(r.spacer_error for r in rows if r.detected),
                median_tir_ident_error=median_or_none(r.tir_ident_error for r in rows if r.detected),
            )
        )

    # false positives = calls not matched to any planted truth
    fp_rows: List[FPRow] = []
    for i, call in enumerate(call_records):
        if i in matched_call_indices:
            continue
        fp_rows.append(
            FPRow(
                call_index=i,
                contig=call.contig,
                start=call.start,
                end=call.end,
                interval_len=call.interval_len,
                arm_len=call.arm_len,
                spacer=call.spacer,
                tir_ident=call.tir_ident,
                identity_est=call.identity_est,
                kept_pts=call.kept_pts,
                matches=call.matches,
                ir_class=call.ir_class,
            )
        )

    return truth_rows, summary_rows, fp_rows


# -----------------------------
# writing
# -----------------------------

def write_tsv(path: Path, rows: List[object]) -> None:
    if not rows:
        # still write an empty file with no rows? easiest: just create empty file
        path.write_text("")
        return

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
    p = argparse.ArgumentParser(description="Evaluate irx calls against a synthetic truth set.")
    p.add_argument("--truth", required=True, help="Truth TSV from make_irx_truthset.py")
    p.add_argument("--calls", required=True, help="irx call TSV")
    p.add_argument("--out_summary", required=True, help="Output per-scenario summary TSV")
    p.add_argument("--out_truth", required=True, help="Output truth-vs-call TSV")
    p.add_argument("--out_fp", required=True, help="Output false-positive TSV")
    p.add_argument(
        "--min_overlap_bp",
        type=int,
        default=100,
        help="Minimum absolute overlap to count as detected [default: 100]",
    )
    p.add_argument(
        "--min_recip_overlap",
        type=float,
        default=0.25,
        help="Minimum reciprocal overlap on both truth and call [default: 0.25]",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    truth = read_truth(Path(args.truth))
    calls = read_calls(Path(args.calls))

    truth_rows, summary_rows, fp_rows = evaluate(
        truth,
        calls,
        min_overlap_bp=args.min_overlap_bp,
        min_recip_overlap=args.min_recip_overlap,
    )

    write_tsv(Path(args.out_summary), summary_rows)
    write_tsv(Path(args.out_truth), truth_rows)
    write_tsv(Path(args.out_fp), fp_rows)

    n_truth = len(truth_rows)
    n_detected = sum(r.detected for r in truth_rows)
    n_calls = len(calls)
    n_fp = len(fp_rows)

    print(f"[ok] truth rows:      {n_truth}")
    print(f"[ok] detected truth:  {n_detected} ({(100*n_detected/n_truth if n_truth else 0):.1f}%)")
    print(f"[ok] calls:           {n_calls}")
    print(f"[ok] unmatched calls: {n_fp}")
    print(f"[ok] wrote summary:   {args.out_summary}")
    print(f"[ok] wrote truth:     {args.out_truth}")
    print(f"[ok] wrote fp:        {args.out_fp}")

    print("\nPer-scenario recovery:")
    for row in summary_rows:
        print(
            f"  {row.scenario_group:16s} "
            f"{row.n_detected:>3d}/{row.n_truth:<3d} "
            f"({100*row.recovery_rate:5.1f}%)"
        )


if __name__ == "__main__":
    main()
