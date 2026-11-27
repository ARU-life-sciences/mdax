#!/usr/bin/env python3
"""Collect the assembly and read-level results the evaluation report reads (results/*.tsv).

Run on the farm head node (reads only small files plus `seqkit stats` on assemblies):
    python3 reports/mdax_eval/collect.py
"""
import os, re, subprocess, glob
U = "/lustre/scratch122/tol/teams/blaxter/users/mb39"
B = f"{U}/mdax_celegans_rerun"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)
busco_re = re.compile(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\]")


def stats(fa):
    if not os.path.exists(fa) or os.path.getsize(fa) == 0:
        return None
    o = subprocess.run(["seqkit", "stats", "-a", "-T", fa], capture_output=True, text=True).stdout.splitlines()[1].split("\t")
    return int(o[3].replace(",", "")), float(o[4]) / 1e6, float(o[12]) / 1e6


def busco(path):
    for p in glob.glob(path):
        m = busco_re.search(open(p).read())
        if m:
            return m.group(1), m.group(3)
    return None


def row(st, bu):
    if not st:
        return ["", "", "", "", ""]
    return [str(st[0]), f"{st[1]:.1f}", f"{st[2]:.2f}", bu[0] if bu else "", bu[1] if bu else ""]


# ---- HiFi, hifiasm -----------------------------------------------------------------
# arm directory name -> label. "now" = newest arm run with the current algorithm for that sample.
hifi = [
    ("nrCaeEleg92", "MDA", "main_base"), ("nrCaeEleg96", "MDA", "cut_resolve"), ("nrCaeEleg98", "MDA", "cut_resolve"),
    ("nrCaeEleg97", "MDA", None), ("nrCaeEleg99", "MDA", None), ("nrCaeEleg100", "MDA", None),
    ("nrCaeEleg34", "standard", None), ("nrCaeEleg_ToL", "standard", None), ("nrCaeEleg33", "standard", None),
]
def old(s, a):
    return row(stats(f"{B}/asm/{s}/{a}/{s}_{a}.p_ctg.fa"), busco(f"{B}/busco/{s}_{a}/short_summary*.txt"))
def new(s, a):
    return row(stats(f"{B}/retest92/{s}.{a}/p_ctg.fa"), busco(f"{B}/retest92/{s}.{a}/busco_summary.txt"))

with open(f"{OUT}/hifi_hifiasm.tsv", "w") as o:
    o.write("sample\tlibrary\tarm\tcontigs\tMb\tN50_Mb\tBUSCO_C\tBUSCO_D\n")
    for s, lib, now in hifi:
        arms = [("raw reads", old(s, "raw")), ("old mdax (Oct)", old(s, "balanced")), ("v5 fold-axis gate", new(s, "v5_balanced"))]
        if lib == "standard":
            arms.append(("calibration fix (v6)", new(s, "v6_balanced")))
        if now:
            arms.append(("main (rescue + cut)", new(s, now)))
        for label, r in arms:
            if r[0]:
                o.write("\t".join([s.replace("nrCaeEleg", ""), lib, label] + r) + "\n")

# ---- HiFi, deep (~720x) libraries: old full-depth arms vs a ~150x subsample ------------
# deep_sub/: mdax v0.2.0 run on the full library, then the same random read subsample
# (seqkit sample -s 11, ~15 Gb) assembled raw and corrected.
def sub(s, a):
    return row(stats(f"{B}/deep_sub/{s}/{a}/p_ctg.fa"), busco(f"{B}/deep_sub/{s}/{a}/busco_summary.txt"))

with open(f"{OUT}/hifi_deep.tsv", "w") as o:
    o.write("sample\tdepth\tarm\tcontigs\tMb\tN50_Mb\tBUSCO_C\tBUSCO_D\n")
    for s in ("nrCaeEleg93", "nrCaeEleg94", "nrCaeEleg95"):
        for depth, label, r in [("full (~720×)", "raw reads", old(s, "raw")), ("full (~720×)", "old mdax (Oct)", old(s, "balanced")),
                                ("150× subsample", "raw reads", sub(s, "raw")), ("150× subsample", "main (v0.2.0)", sub(s, "mdax"))]:
            if r[0]:
                o.write("\t".join([s.replace("nrCaeEleg", ""), depth, label] + r) + "\n")

# ---- HiFi oracle / rejected arms -----------------------------------------------------
with open(f"{OUT}/hifi_experiments.tsv", "w") as o:
    o.write("sample\tarm\tcontigs\tMb\tN50_Mb\tBUSCO_C\tBUSCO_D\n")
    for s in ("nrCaeEleg92", "nrCaeEleg96"):
        for label, a in [("main", "o_base"), ("+ second-junction recursion (rejected)", "multij"),
                         ("+ oracle: no-copy foldbacks cut", "o_nocopy"), ("+ oracle: every residual foldback cut", "o_all"),
                         ("+ oracle: every one cut, keep all pieces", "o_allpieces")]:
            r = new(s, a)
            if r[0]:
                o.write("\t".join([s.replace("nrCaeEleg", ""), label] + r) + "\n")

# ---- ONT / CLR, Flye ------------------------------------------------------------------
SN = f"{U}/mdax_single_nematode_benchmark"; GA = f"{U}/mdax_gastrotricha_benchmark"; RV = f"{U}/mdax_retest_v4"
def sn(root, sp, arm):
    st = stats(f"{root}/asm/{sp}/{arm}/assembly.fasta")
    bu = None
    for p in glob.glob(f"{root}/busco/{sp}_{arm}*/**/short_summary*.txt", recursive=True) + glob.glob(f"{root}/busco/{sp}_{arm}*/short_summary*.txt"):
        bu = busco(p)
        if bu: break
    return row(st, bu)
def rv(d):
    return row(stats(f"{RV}/{d}/assembly.fasta"), busco(f"{RV}/{d}/busco_summary.txt"))

ont = [  # species, label, data, root, current-arm dirs in mdax_retest_v4 (label, dir)
    ("C. elegans (celegans_wga)", "ONT", SN, "celegans_wga", [("v4 fold-axis gate", "celegans_wga"), ("main (rescue + cut)", "celegans_wga_cut")]),
    ("A. besseyi (aphelenchoides)", "ONT", SN, "aphelenchoides", [("v4 fold-axis gate", "aphelenchoides"), ("main (rescue + cut)", "aphelenchoides_v9")]),
    ("L. squamata SRR27483171", "PacBio CLR", GA, "SRR27483171", [("v4 fold-axis gate", "SRR27483171")]),
    ("L. squamata SRR27483172", "PacBio CLR", GA, "SRR27483172", [("v4 fold-axis gate", "SRR27483172")]),
    ("Sabatieria", "ONT", SN, "sabatieria", [("v4 fold-axis gate", "sabatieria")]),
    ("Microlaimidae", "ONT", SN, "microlaimidae", [("v4 fold-axis gate", "microlaimidae")]),
]
others = ["enoplolaimus", "mesodorylaimus", "rhynchonema", "theristus", "ptycholaimellus", "trileptium", "linhomoeus",
          "paralinhomoeus", "trissonchulus_lat", "trissonchulus_wlg"]
with open(f"{OUT}/ont_flye.tsv", "w") as o:
    o.write("species\tdata\tarm\tcontigs\tMb\tN50_Mb\tBUSCO_C\tBUSCO_D\n")
    for sp, data, root, d, cur in ont:
        for label, arm in [("raw reads", "raw"), ("Lee et al. remove_palindrome", "theirs"), ("old mdax (Oct)", "mdax")]:
            r = sn(root, d, arm)
            if r[0]: o.write("\t".join([sp, data, label] + r) + "\n")
        for label, dd in cur:
            r = rv(dd)
            if r[0]: o.write("\t".join([sp, data, label] + r) + "\n")
with open(f"{OUT}/ont_flye_suite.tsv", "w") as o:
    o.write("species\traw_contigs\traw_N50\ttheirs_contigs\ttheirs_N50\tmdax_contigs\tmdax_N50\traw_C\ttheirs_C\tmdax_C\n")
    for sp in others:
        r, t, m = sn(SN, sp, "raw"), sn(SN, sp, "theirs"), sn(SN, sp, "mdax")
        o.write("\t".join([sp, r[0], r[2], t[0], t[2], m[0], m[2], r[3], t[3], m[3]]) + "\n")

# ---- read level ------------------------------------------------------------------------
D = f"{B}/recall_check"
line_re = lambda cls: re.compile(rf"^\s+{cls}\s+n=\s*(\d+)\s+(.*)$", re.M)
def outcome(path, cls):
    if not os.path.exists(path): return None
    m = line_re(cls).search(open(path).read())
    if not m: return None
    d = dict(x.split("=") for x in m.group(2).split())
    return m.group(1), d.get("artefact", "0%").rstrip("%"), d.get("not_detected", "0%").rstrip("%"), d.get("real", "0%").rstrip("%")
sets = [("nrCaeEleg92", "HiFi", [("old mdax (Oct)", "nrCaeEleg92"), ("v5 fold-axis gate", "v5_default"), ("main", "regress/v9.nrCaeEleg92")]),
        ("nrCaeEleg94", "HiFi", [("old mdax (Oct)", "nrCaeEleg94"), ("v5 fold-axis gate", "v6_94"), ("main", "regress/v9.nrCaeEleg94")]),
        ("celegans_wga", "ONT", [("old mdax (2 Nov)", "ont_old"), ("v4 fold-axis gate", "ont_v4"), ("main", "regress/v9.celegans_wga")]),
        ("aphelenchoides", "ONT", [("old mdax (2 Nov)", "aph_old"), ("v5 fold-axis gate", "aph_v5"), ("main", "regress/v9.aphelenchoides")])]
with open(f"{OUT}/read_level.tsv", "w") as o:
    o.write("truth set\tdata\tbuild\ttrue foldbacks\tcut %\tmissed %\tkept as real %\tlinear reads\tlinear cut %\n")
    for s, data, builds in sets:
        for label, f in builds:
            fb, li = outcome(f"{D}/{f}.result.txt", "foldback"), outcome(f"{D}/{f}.result.txt", "linear")
            if fb and li:
                o.write("\t".join([s, data, label, fb[0], fb[1], fb[2], fb[3], li[0], li[1]]) + "\n")
print("wrote", sorted(os.listdir(OUT)))
