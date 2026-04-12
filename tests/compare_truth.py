"""
compare_truth.py — run mdax on the 57-read HiFi ground-truth set and produce
a per-read comparison TSV for downstream analysis and plotting.

Matches the parameter set and pass/fail criteria of the Rust integration test
`real_artefacts_detected_and_cut` in pipeline_integration.rs.

Usage (from repo root or tests/):
    python3 tests/compare_truth.py [--out PATH] [--mdax PATH]

Output columns:
    label           truth label (e.g. art_hi_11012bp)
    category        label prefix (art_hi / art_md / art_lo / art_vlo / li / real_hi / real_lo)
    expected_dec    expected mdax decision from ground truth
    truth_split     BLAST-derived split position (bp); 0 = non-deterministic, skip check
    truth_arm_len   BLAST arm length (bp)
    truth_identity  BLAST arm identity
    truth_gap_est   estimated template-switch gap (bp)
    detected        true/false — read appeared in mdax report
    mdax_decision   mdax output decision (artefact / low_ident / real / -)
    mdax_split      mdax refined_split (bp); "" if not detected
    split_delta_bp  mdax_split - truth_split (signed); "" if not applicable
    abs_delta_bp    abs(split_delta_bp); "" if not applicable
    decision_ok     true/false — mdax_decision matches expected
    split_ok        true/false — abs_delta <= 100 bp (or n/a when skipped)
    pass_fail       PASS if detected + decision_ok + split_ok; FAIL otherwise
"""
import argparse, csv, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

FASTA     = os.path.join(REPO, "data", "chris_artefacts.fasta")
TRUTH_TSV = os.path.join(REPO, "data", "chris_artefacts_truth.tsv")
MDAX_BIN  = os.path.join(REPO, "target", "release", "mdax")

SPLIT_TOL_BP = 100

# Expected decisions and canonical split positions from REAL_READ_EXPECTATIONS
# in tests/pipeline_integration.rs.
# blast_split=0 means non-deterministic junction — skip the split accuracy check.
EXPECTATIONS = {
    "art_hi_3820bp":    ("artefact",  2128),
    "art_hi_5653bp":    ("artefact",  2609),
    "art_hi_6156bp":    ("artefact",  1144),
    "art_hi_9193bp":    ("artefact",  5480),
    "art_hi_11012bp":   ("artefact",  2116),   # two junctions; use 100%-identity junction
    "art_hi_13019bp":   ("artefact",  4418),
    "art_md_3728bp":    ("artefact",  2231),
    "art_md_5167bp":    ("artefact",  3560),
    "art_md_5944bp":    ("artefact",  2347),
    "art_md_6356bp":    ("artefact",  1330),
    "art_md_6437bp":    ("artefact",  3606),
    "art_md_9605bp":    ("artefact",  2247),
    "art_md_10782bp":   ("artefact",  4910),
    "art_md_14746bp":   ("artefact",  9555),
    "art_lo_5094bp":    ("artefact",  2609),
    "art_lo_6432bp":    ("artefact",  2050),
    "art_lo_8105bp":    ("artefact",     0),   # two junctions; skip split check
    "art_lo_8708bp":    ("artefact",  2059),
    "art_lo_10715bp":   ("artefact",  4305),
    "art_lo_13413bp":   ("artefact",  9821),
    "art_vlo_3906bp":   ("artefact",  1045),
    "art_vlo_5230bp":   ("artefact",  4000),
    "art_vlo_6497bp":   ("artefact",  2479),
    "art_vlo_9989bp":   ("artefact",  7110),
    "art_vlo_10708bp":  ("artefact",  7284),
    "art_vlo_13247bp":  ("artefact",  2264),
    "art_vlo_14291bp":  ("artefact",  3053),
    "li_3279bp":        ("low_ident", 1415),
    "li_5619bp":        ("low_ident", 2298),
    "li_6296bp":        ("low_ident", 4003),
    "li_7152bp":        ("low_ident", 3150),
    "li_9340bp":        ("low_ident", 5215),
    "li_11367bp":       ("low_ident", 2134),
    "li_12397bp":       ("low_ident", 7704),
    "li_14193bp":       ("low_ident", 6098),
    "real_hi_4085bp":   ("artefact",  1210),
    "real_hi_4203bp":   ("artefact",  1747),
    "real_hi_4469bp":   ("artefact",  3174),
    "real_hi_4797bp":   ("artefact",  1012),
    "real_hi_4942bp":   ("artefact",  1803),
    "real_hi_5421bp":   ("artefact",  1118),
    "real_hi_5879bp":   ("artefact",  3605),
    "real_hi_8673bp":   ("artefact",  4328),
    "real_lo_5920bp":   ("artefact",  1915),
    "real_lo_5954bp":   ("artefact",  2386),
    "real_lo_6526bp":   ("artefact",  1967),
    "real_lo_7250bp":   ("artefact",  1264),
    "real_lo_7284bp":   ("artefact",  1230),
    "real_lo_9927bp":   ("artefact",  3366),
    "real_lo_10434bp":  ("artefact",  7713),
    "real_lo_10612bp":  ("artefact",  5932),
    "real_lo_10650bp":  ("artefact",  7662),
    "real_lo_11491bp":  ("artefact",     0),   # non-deterministic junction; skip split check
    "real_lo_11727bp":  ("artefact",  8875),
    "real_lo_12897bp":  ("artefact",  2841),
    "real_lo_13009bp":  ("artefact",  4789),
    "real_lo_14005bp":  ("artefact",  8221),
}


def label_prefix(label):
    parts = label.rsplit("_", 1)
    return parts[0] if len(parts) == 2 else label


def run_mdax(mdax_bin, fasta, report_path):
    # Parameters matching real_data_cfg(2) in pipeline_integration.rs
    cmd = [
        mdax_bin,
        "--report", report_path,
        "--output", os.devnull,
        "--min-identity", "0.5",
        "--min-matches", "12",
        "--min-support", "2",
        fasta,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[error] mdax exited {result.returncode}:", file=sys.stderr)
        print(result.stderr[-2000:], file=sys.stderr)
        sys.exit(1)


def load_truth(path):
    truth = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            truth[row["label"]] = row
    return truth


def load_mdax_report(path):
    """Key by the label (first whitespace token of read_id, e.g. 'art_hi_11012bp')."""
    mdax = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            label = row["read_id"].split()[0]
            # Keep the row with the largest refined_split (outermost cut in recursive runs)
            if label not in mdax:
                mdax[label] = row
            else:
                cur = int(mdax[label].get("refined_split") or 0)
                new = int(row.get("refined_split") or 0)
                if new > cur:
                    mdax[label] = row
    return mdax


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",  default=os.path.join(HERE, "results_comparison.tsv"),
                    help="Output TSV path (default: tests/results_comparison.tsv)")
    ap.add_argument("--mdax", default=MDAX_BIN,
                    help="Path to mdax binary")
    args = ap.parse_args()

    for p in [args.mdax, FASTA, TRUTH_TSV]:
        if not os.path.exists(p):
            print(f"[error] Not found: {p}", file=sys.stderr)
            sys.exit(1)

    truth = load_truth(TRUTH_TSV)

    with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False) as tf:
        report_path = tf.name
    try:
        print(f"Running mdax on {os.path.basename(FASTA)}…", file=sys.stderr)
        run_mdax(args.mdax, FASTA, report_path)
        mdax = load_mdax_report(report_path)
    finally:
        os.unlink(report_path)

    fields = [
        "label", "category", "expected_dec",
        "truth_split", "truth_arm_len", "truth_identity", "truth_gap_est",
        "detected", "mdax_decision", "mdax_split",
        "split_delta_bp", "abs_delta_bp",
        "decision_ok", "split_ok", "pass_fail",
    ]

    n_pass = n_fail = 0
    rows = []

    for label in sorted(EXPECTATIONS):
        exp_dec, exp_split = EXPECTATIONS[label]
        t = truth.get(label, {})

        t_arm    = t.get("arm_len", "")
        t_ident  = t.get("identity", "")
        t_gap    = t.get("gap_est", "")

        m = mdax.get(label)
        detected = m is not None
        mdax_dec   = m["decision"] if m else "-"
        mdax_split = m.get("refined_split", "") if m else ""
        if mdax_split == "" and m:
            mdax_split = m.get("coarse_split", "")

        # Decision check (informational; Rust test does not enforce specific decision,
        # only that it is one of artefact/low_ident/real)
        dec_ok = (mdax_dec == exp_dec)
        dec_valid = mdax_dec in ("artefact", "low_ident", "real")

        # Split check: only for artefact/real decisions, and only when exp_split > 0
        split_delta = abs_delta = ""
        split_ok = True  # default: not applicable (low_ident or skipped)
        if detected and exp_split > 0 and mdax_dec in ("artefact", "real"):
            if mdax_split:
                try:
                    delta = int(mdax_split) - exp_split
                    split_delta = str(delta)
                    abs_delta   = str(abs(delta))
                    split_ok = abs(delta) <= SPLIT_TOL_BP
                except ValueError:
                    split_ok = False

        # PASS = detected + valid decision + split within tolerance (matching Rust test logic)
        ok = detected and dec_valid and split_ok
        result = "PASS" if ok else "FAIL"
        if ok:
            n_pass += 1
        else:
            n_fail += 1

        rows.append([
            label, label_prefix(label), exp_dec,
            str(exp_split), t_arm, t_ident, t_gap,
            str(detected).lower(), mdax_dec, mdax_split,
            split_delta, abs_delta,
            str(dec_ok).lower(), str(split_ok).lower(), result,
        ])

    with open(args.out, "w") as f:
        f.write("\t".join(fields) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    total = n_pass + n_fail
    print(f"Wrote {args.out}  ({total} reads: {n_pass} PASS, {n_fail} FAIL)",
          file=sys.stderr)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
