# Changelog

## v0.2.0 (2026-11-22)

This release changes what mdax cuts by default. On MDA HiFi it now cuts 74–79% of reference-confirmed foldback
chimeras (v0.1.x: 16–38%) while cutting about 2% of normal reads, and hifiasm assemblies of six amplified
C. elegans libraries have 40–62% fewer contigs and two to four times the N50 of assemblies from raw reads. See
`reports/mdax_eval/mdax_eval.pdf` for the full evaluation.

### Changed
- **Fold-axis support gate.** A foldback is kept as a genuine palindrome only if enough other reads fold at the
  same axis (same minimizer pair and offset within a tolerance), not merely in the same locus. In deep MDA
  libraries, chance co-location of artefacts at one locus made the old locus signatures keep 35–52% of artefacts
  as "real". The threshold is calibrated per library from the support histogram; on libraries with few
  singleton foldbacks (not amplified) calibration is skipped.
- **Balanced defaults widened:** `min_matches` 12, `min_span` 1000, `end_guard` 500.
- **Cut placement.** Displaced foldbacks are cut at the strand switch, not past the un-mirrored spacer between
  the arms (HiFi: spacer side resolved from library k-mer counts; ONT: spacer always dropped). Leftover foldbacks
  in cut reads fell from 62–66% to 33%.
- **Low-identity foldbacks are cut by default** (subject to the support gate), keeping the longer side. Use
  `--keep-low-ident` for the old behaviour.

### Added
- **Rescue detection** of short-arm and read-end foldbacks, confirmed by banded alignment
  (`--no-rescue`, `--rescue-*` to tune). ONT foldbacks cut: 65% to 77.5%.
- `[SUPPORT CALIBRATION]` and `[QC SUMMARY]` lines on stderr.
- Reference-based regression test (`tests/truth_regression.rs`, fixture in `data/ci_truth/`, public SRA reads),
  a farm-side scorer over four truth sets (`synth_bench/truth_regress.sh`) and a GitHub Actions workflow.
- `.cargo/config.toml`, so a plain `cargo build --release` works (AVX2 and AES-NI are required).

### Deprecated
- `--sig-flank-bp`, `--sig-take`, `--sig-value-shift` and `--cut-low-ident` are accepted but do nothing;
  each prints a note.
- The `support_rank_frac` report column is always 0 and is kept only so column positions do not change.

### Removed
- The locus band-signature support code (about 1,300 lines). Output was verified byte-identical before and after.

### Notes
- About 20% of true foldbacks have no inverted copy in the read (the strand switch lands further away than the
  read extends), so self-comparison cannot find them. Cutting them at reference-detected junctions gave no
  assembly gain, so this is not a target.
- On libraries that were not amplified, raw reads assemble marginally better; mdax is safe to run but not useful.
