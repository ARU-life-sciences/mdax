# data/

Test and development data for `mdax`.

---

## Real HiFi dataset (primary test suite)

| File | Description |
|------|-------------|
| `chris_artefacts.fasta` | 57 HiFi CCS reads curated from `ERR12263839` (PacBio HiFi WGA). Contains a mix of confirmed MDA foldback artefacts (art_hi_*, art_lo_*, art_vlo_*, li_*) and genuine reads (real_*). This is the input used by the integration test suite. |
| `chris_artefacts_truth.tsv` | Ground-truth characterisation of the 57 reads, produced independently of `mdax` by `extract_and_characterise.py`. Columns: `label`, `read_id`, `length`, `split_pos`, `arm_len`, `identity`, `gap_est`, `blast_hits`, `mdax_prior_decision`. |
| `chris_artefacts_raw.fastq` | The same 57 reads in FASTQ format (with quality scores), as extracted from the original run. Source for `extract_and_characterise.py`. |

Read name conventions:
- `art_hi_*` — artefact, high identity (≥0.99), zero gap
- `art_lo_*` / `art_vlo_*` — artefact, lower arm identity or very short arm
- `li_*` — artefact, large template-switch gap (1.6–4 kb), medium identity (~0.78)
- `real_*` — genuine biological read (not an artefact)

`chris_artefacts.fasta` and `chris_artefacts_truth.tsv` are the canonical names used by the integration tests in `tests/pipeline_integration.rs`.

---

## Synthetic datasets

| File | Generator | Description |
|------|-----------|-------------|
| `test_foldback.fasta` | `synthesise.py` | 5 synthetic reads: a single artefact foldback, two identical true palindromes, a normal read, and an 8000 bp foldback. Used by the unit tests. |
| `test_reads.fa` | `synthesise_concatemers.py` | 4 synthetic concatemer reads: `forward_concatemer`, `inverted_concatemer`, `template_switch`, `multi_switch`. |
| `truth.tsv` | — | Ground truth for the 4 concatemer reads in `test_reads.fa`. |

---

## ONT dataset

| File | Description |
|------|-------------|
| `cadect_test_file.fasta` | 8 ONT reads from a CADECT WGA experiment (flow cell FAU80014, barcode15). Useful for manual ONT-mode testing. |
| `cadect_test_file.fastq` | Same reads with quality scores. |

---

## Scripts

| File | Description |
|------|-------------|
| `extract_and_characterise.py` | Extracts and independently characterises MDA foldback reads from a raw FASTQ + TSV pair. Produces `chris_artefacts.fasta` and `chris_artefacts_truth.tsv`. Run with `--fastq <ERR12263839.fastq.gz> --tsv <ERR12263839.tsv> --out-fasta chris_artefacts.fasta --out-truth chris_artefacts_truth.tsv`. |
| `synthesise.py` | Generates `test_foldback.fasta` (synthetic foldback reads, seed 42). |
| `synthesise_concatemers.py` | Generates `test_reads.fa` (synthetic concatemer reads). |

---

## Notes

- `__pycache__/` is a Python cache directory (gitignored) and can be safely deleted if present.
