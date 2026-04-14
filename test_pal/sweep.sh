#!/usr/bin/env bash
# sweep.sh — parameter sensitivity analysis for mdax on the SRR24201687 ONT subsample.
#
# For each parameter combination, runs mdax, normalises calls, and collects
# precision/recall/F1 (strict and generous) into a formatted table.
#
# Usage:
#   ./sweep.sh                   # run all sweep groups
#   ./sweep.sh modes             # only the --mode preset sweep
#   ./sweep.sh identity          # only the --min-identity sweep
#   ./sweep.sh refine            # only the --refine-mode sweep
#   ./sweep.sh matches           # only the --min-matches sweep
#   ./sweep.sh support           # only the --min-support sweep
#   ./sweep.sh combined          # curated combos (ont + low identity etc.)
#
# Output: a TSV table printed to stdout; per-run mdax stderr is suppressed.
# Cached runs: results are written to sweep_cache/<tag>.tsv and reused if present.
#              Delete sweep_cache/ to force a full re-run.
#
# Requires: ../target/release/mdax  (cargo build --release)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MDAX_BIN="${MDAX_BIN:-$SCRIPT_DIR/../target/release/mdax}"
INPUT_FQ="$SCRIPT_DIR/out/subsample/SRR24201687.first20k.fastq.gz"
BASELINE_DIR="$SCRIPT_DIR/out/baseline/SRR24201687_first20k"
CACHE_DIR="$SCRIPT_DIR/sweep_cache"
TOTAL_READS=20000
THREADS="${THREADS:-8}"

mkdir -p "$CACHE_DIR"

[[ -x "$MDAX_BIN" ]]   || { echo "mdax binary not found: $MDAX_BIN"; exit 1; }
[[ -f "$INPUT_FQ"  ]]  || { echo "Input not found: $INPUT_FQ"; exit 1; }

# ── run one mdax configuration ─────────────────────────────────────────────────
# Args: TAG MDAX_EXTRA_FLAGS...
# Writes normalised calls to sweep_cache/<TAG>.tsv (cached).
# Prints a single TSV metrics row to stdout.
run_one() {
    local tag="$1"; shift
    # Sanitize tag for use as a filename (replace any non-alphanumeric chars with _)
    local safe_tag="${tag//[^a-zA-Z0-9=.,_-]/_}"
    local cached="$CACHE_DIR/${safe_tag}.calls.tsv"

    if [[ ! -f "$cached" ]]; then
        local report="$CACHE_DIR/${tag}.report.tsv"
        "$MDAX_BIN" --threads "$THREADS" "$@" \
            --report  "$report" \
            --output  /dev/null \
            "$INPUT_FQ" 2>/dev/null

        # Normalise: resolve column indices from the header row so this works
        # regardless of whether the binary was built with support_rank_frac or not.
        awk -F'\t' '
          NR==1{
            for(i=1;i<=NF;i++){
              if ($i=="refined_split") rs_col=i;
              if ($i=="decision")      dec_col=i;
            }
            next
          }
          {
            split($1, a, " "); id=a[1];
            sub(/-dup[0-9]+(\.[0-9]+)?$/, "", id);
            rs  = rs_col  ? $rs_col  : "";
            dec = dec_col ? $dec_col : "";
            called = (dec=="artefact") ? 1 : 0;
            print id "\t" called "\t" rs "\t" dec
          }
        ' "$report" | LC_ALL=C sort -u > "$cached"
    fi

    # Run test_pal.py against the cached calls
    (
        cd "$SCRIPT_DIR"
        # Point test_pal.py at the cached calls file by temporarily symlinking
        ln -sf "$cached" "_sweep_calls_tmp.tsv"
        python3 - "$tag" "$TOTAL_READS" <<'PYEOF'
import sys, csv, statistics
from collections import Counter

tag = sys.argv[1]
total_reads = int(sys.argv[2])

mdax = {}
dec_counts = Counter()
with open("_sweep_calls_tmp.tsv") as f:
    for r in csv.reader(f, delimiter="\t"):
        mdax[r[0]] = {"called": int(r[1]), "split": r[2], "decision": r[3] if len(r)>3 else ""}
        dec_counts[mdax[r[0]]["decision"]] += 1

baseline_called = set()
with open("baseline.flagged.canon.ids") as f:
    for l in f:
        rid = l.strip()
        if rid: baseline_called.add(rid)

def compute_cm(called_decs):
    tp=fp=fn=tn=0
    called_ids = {rid for rid,v in mdax.items() if v["decision"] in called_decs}
    known = set(mdax)|baseline_called
    for rid in known:
        m = 1 if rid in called_ids else 0
        b = 1 if rid in baseline_called else 0
        if m and b: tp+=1
        elif m:     fp+=1
        elif b:     fn+=1
        else:       tn+=1
    tn += total_reads - len(known)
    return tp,fp,fn,tn

def metrics(tp,fp,fn,tn):
    def d(a,b): return a/b if b else float("nan")
    return d(tp,tp+fp), d(tp,tp+fn), d(2*tp,2*tp+fp+fn), d(tp,tp+fp+fn)

tp_s,fp_s,fn_s,tn_s = compute_cm({"artefact"})
tp_g,fp_g,fn_g,tn_g = compute_cm({"artefact","low_ident"})

ps,rs,f1s,js = metrics(tp_s,fp_s,fn_s,tn_s)
pg,rg,f1g,jg = metrics(tp_g,fp_g,fn_g,tn_g)

n_det = len(mdax)
n_art = dec_counts["artefact"]
n_low = dec_counts["low_ident"]
n_real= dec_counts["real"]

# baseline reads not detected at all
detected_baseline = sum(1 for rid in mdax if rid in baseline_called)
missed = len(baseline_called) - detected_baseline

print("\t".join([
    tag,
    f"{ps:.3f}", f"{rs:.3f}", f"{f1s:.3f}", f"{js:.3f}",
    f"{pg:.3f}", f"{rg:.3f}", f"{f1g:.3f}", f"{jg:.3f}",
    str(n_det), str(n_art), str(n_low), str(n_real),
    str(missed),
]))
PYEOF
        rm -f "_sweep_calls_tmp.tsv"
    )
}

# ── print helpers ─────────────────────────────────────────────────────────────

GROUP="${1:-all}"

SUMMARY_TSV="$SCRIPT_DIR/sweep_summary.tsv"

# Always write/overwrite header so sweep_summary.tsv is always valid
printf "tag\tstrict_P\tstrict_R\tstrict_F1\tstrict_J\tgenerous_P\tgenerous_R\tgenerous_F1\tgenerous_J\tdetected\tartefact\tlow_ident\treal\tmissed\n" \
    > "$SUMMARY_TSV"

_print_header() {
    printf "%-40s  %5s %5s %5s %5s  |  %5s %5s %5s %5s  | %6s %6s %6s %6s  %6s\n" \
        "tag" "sP" "sR" "sF1" "sJ" "gP" "gR" "gF1" "gJ" \
        "det" "art" "low" "real" "missed"
    printf "%s\n" "$(printf '─%.0s' {1..105})"
}

# begin_group NAME — print group header when the group matches; return 0 if
# this group should run, 1 if it should be skipped.
begin_group() {
    local name="$1"
    [[ "$GROUP" == "all" || "$GROUP" == "$name" ]] || return 1
    echo
    echo "══════════════════════════════════════════════════════════════════════════════"
    echo "  SWEEP: $name"
    echo "══════════════════════════════════════════════════════════════════════════════"
    echo "  s=strict call (artefact only) | g=generous call (artefact+low_ident)"
    echo "  P=Precision  R=Recall  F1  J=Jaccard  det=detected  art=artefact  low=low_ident  missed=baseline reads not detected"
    echo
    _print_header
    return 0
}

# emit TAG MDAX_FLAGS... — run one config, print its formatted row, and append to summary TSV.
emit() {
    local row
    row="$(run_one "$@")"
    IFS=$'\t' read -r tag sp sr sf1 sj gp gr gf1 gj det art low real missed <<< "$row"
    printf "%-40s  %5s %5s %5s %5s  |  %5s %5s %5s %5s  | %6s %6s %6s %6s  %6s\n" \
        "$tag" "$sp" "$sr" "$sf1" "$sj" "$gp" "$gr" "$gf1" "$gj" \
        "$det" "$art" "$low" "$real" "$missed"
    # Append to machine-readable summary (append so partial sweeps accumulate)
    printf "%s\n" "$row" >> "$SUMMARY_TSV"
}

# ── sweep groups ───────────────────────────────────────────────────────────────

# group: mode presets
sweep_modes() {
    begin_group "modes" || return 0
    for mode in strict balanced permissive; do
        emit "mode=$mode" --mode "$mode"
    done
}

# group: refine-mode (hifi vs ont) at each preset
sweep_refine() {
    begin_group "refine" || return 0
    for mode in strict balanced permissive; do
        for rmode in hifi ont; do
            emit "mode=$mode,refine=$rmode" --mode "$mode" --refine-mode "$rmode"
        done
    done
}

# group: min-identity sweep (balanced base, both refine modes)
sweep_identity() {
    begin_group "identity" || return 0
    for rmode in hifi ont; do
        for ident in 0.30 0.33 0.36 0.40 0.45 0.50 0.55 0.60 0.65 0.70; do
            emit "balanced,refine=$rmode,ident=$ident" \
                --mode balanced --refine-mode "$rmode" --min-identity "$ident"
        done
    done
}

# group: min-matches sweep (balanced, hifi)
sweep_matches() {
    begin_group "matches" || return 0
    for mm in 8 10 12 15 18 20 25 30; do
        emit "balanced,hifi,min_matches=$mm" \
            --mode balanced --refine-mode hifi --min-matches "$mm"
    done
}

# group: min-support sweep (balanced, hifi, identity=0.45)
sweep_support() {
    begin_group "support" || return 0
    for sup in 2 3 4 5 6; do
        emit "balanced,hifi,ident=0.45,support=$sup" \
            --mode balanced --refine-mode hifi --min-identity 0.45 --min-support "$sup"
    done
}

# group: refine-arm sweep (ont mode, two identity thresholds)
# Larger arm gives the Levenshtein aligner more context on noisy ONT reads,
# which raises identity estimates and converts low_ident → artefact calls.
sweep_arm() {
    begin_group "arm" || return 0
    for ident in 0.35 0.40; do
        for arm in 200 400 600 800 1000; do
            emit "ont,ident=$ident,arm=$arm" \
                --mode balanced --refine-mode ont \
                --min-identity "$ident" --refine-arm "$arm"
        done
    done
}

# group: curated combos — best guesses for ONT + low-ident data
sweep_combined() {
    begin_group "combined" || return 0
    # Current defaults (baseline reference)
    emit "balanced,hifi [DEFAULT]" \
        --mode balanced --refine-mode hifi

    # ONT mode only
    emit "balanced,ont" \
        --mode balanced --refine-mode ont

    # Permissive
    emit "permissive,hifi" \
        --mode permissive --refine-mode hifi

    emit "permissive,ont" \
        --mode permissive --refine-mode ont

    # Balanced + ont + lower identity thresholds
    for ident in 0.35 0.38 0.40 0.42 0.45; do
        emit "balanced,ont,ident=$ident" \
            --mode balanced --refine-mode ont --min-identity "$ident"
    done

    # Permissive + ont + lower identity
    for ident in 0.30 0.33 0.35; do
        emit "permissive,ont,ident=$ident" \
            --mode permissive --refine-mode ont --min-identity "$ident"
    done

    # More sensitive matching (lower min-matches) + ont
    emit "balanced,ont,ident=0.40,mm=15" \
        --mode balanced --refine-mode ont --min-identity 0.40 --min-matches 15
    emit "balanced,ont,ident=0.40,mm=12" \
        --mode balanced --refine-mode ont --min-identity 0.40 --min-matches 12

    # Bigger refine arm — more context for ONT Levenshtein alignment
    for arm in 400 600 1000; do
        emit "balanced,ont,ident=0.40,arm=$arm" \
            --mode balanced --refine-mode ont --min-identity 0.40 --refine-arm "$arm"
    done
    for arm in 400 600 1000; do
        emit "balanced,ont,ident=0.35,arm=$arm" \
            --mode balanced --refine-mode ont --min-identity 0.35 --refine-arm "$arm"
    done

    # Best candidate: permissive + ont + low ident + large arm
    for arm in 600 1000; do
        emit "permissive,ont,ident=0.35,arm=$arm" \
            --mode permissive --refine-mode ont --min-identity 0.35 --refine-arm "$arm"
    done
}

# group: --cut-low-ident effect on ONT data
# Compares default low_ident behaviour vs cutting low-identity reads for ONT mode.
# Shows impact of Priority 2 fix across a range of arm sizes and identity thresholds.
sweep_cut_low() {
    begin_group "cut_low" || return 0
    # Baseline: balanced+hifi (default, no cut_low_ident)
    emit "balanced,hifi,default" \
        --mode balanced --refine-mode hifi

    # ONT without cut_low_ident (show current state after ONT-fallback fix only)
    for ident in 0.35 0.40; do
        emit "ont,ident=$ident,no_cut_low" \
            --mode balanced --refine-mode ont --min-identity "$ident"
    done

    # ONT with cut_low_ident (cut reads below identity threshold too)
    for ident in 0.35 0.40; do
        emit "ont,ident=$ident,cut_low" \
            --mode balanced --refine-mode ont --min-identity "$ident" --cut-low-ident
    done

    # Same with larger arm for better ONT identity estimates
    for ident in 0.35 0.40; do
        emit "ont,ident=$ident,arm=600,cut_low" \
            --mode balanced --refine-mode ont --min-identity "$ident" \
            --refine-arm 600 --cut-low-ident
    done

    # Permissive + cut_low_ident
    emit "permissive,ont,ident=0.35,cut_low" \
        --mode permissive --refine-mode ont --min-identity 0.35 --cut-low-ident
    emit "permissive,ont,ident=0.35,arm=600,cut_low" \
        --mode permissive --refine-mode ont --min-identity 0.35 \
        --refine-arm 600 --cut-low-ident
}

# group: k-mer size sweep for ONT coarse detection sensitivity
# The 26% missed reads fail detect_foldback entirely: k=17 requires exact
# 17-mer matches, but with ~10% ONT error each arm copy has only ~2% survival
# probability per k-mer, producing ~2 matches on a 2000 bp arm vs min_matches=12.
# Smaller k raises this probability and produces more matches.
# Best results will likely use --cut-low-ident since smaller k lowers identity estimates.
sweep_k() {
    begin_group "k" || return 0
    # Baseline with default k=17 w=21 (permissive + cut_low to show ceiling)
    emit "k=17,w=21,permissive,cut_low" \
        -k 17 -w 21 --mode permissive --refine-mode ont --cut-low-ident

    # Sweep k with w=11 (odd, more minimizers than w=21)
    for k in 13 11 9; do
        emit "k=$k,w=11,permissive,cut_low" \
            -k "$k" -w 11 --mode permissive --refine-mode ont --cut-low-ident
    done

    # Also try w=7 for most sensitive detection
    for k in 11 9; do
        emit "k=$k,w=7,permissive,cut_low" \
            -k "$k" -w 7 --mode permissive --refine-mode ont --cut-low-ident
    done

    # k=9 balanced (sanity check that precision holds)
    emit "k=9,w=7,balanced,cut_low" \
        -k 9 -w 7 --mode balanced --refine-mode ont --cut-low-ident

    # k=9 without cut_low_ident to see identity picture
    emit "k=9,w=7,permissive,no_cut_low" \
        -k 9 -w 7 --mode permissive --refine-mode ont --min-identity 0.35
}

# group: final recommended configurations using new binary defaults
# (auto k=11/w=11 for ONT, weighted split estimator, --cut-low-ident)
sweep_final() {
    begin_group "final" || return 0
    # Old default (HiFi): reference baseline
    emit "hifi,balanced,default" \
        --mode balanced --refine-mode hifi

    # New recommended ONT flags: auto k=11 w=11, cut_low_ident
    emit "ont,permissive,cut_low [RECOMMENDED]" \
        --mode permissive --refine-mode ont --cut-low-ident

    emit "ont,balanced,cut_low" \
        --mode balanced --refine-mode ont --cut-low-ident

    # Permissive without cut_low to show low/artefact split
    emit "ont,permissive,no_cut_low" \
        --mode permissive --refine-mode ont --min-identity 0.35
}

sweep_modes
sweep_refine
sweep_identity
sweep_matches
sweep_support
sweep_arm
sweep_combined
sweep_cut_low
sweep_k
sweep_final

echo
echo "Cache: $CACHE_DIR"
echo "Delete cache to force re-runs: rm -rf $CACHE_DIR"
echo
