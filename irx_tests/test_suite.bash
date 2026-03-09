# run the test suite
python3 make_test_suite.py --out-fasta irx_truth.fa --out-truth irx_truth.tsv --out-manifest irx_manifest.tsv --seed 42

# run irx
../target/release/irx --version
../target/release/irx -b irx_calls.tsv irx_truth.fa

# evaluate
python3 eval_irx_test_suite.py \
  --truth irx_truth.tsv \
  --calls irx_calls.tsv \
  --out_summary irx_eval.summary.tsv \
  --out_truth irx_eval.truth.tsv \
  --out_fp irx_eval.fp.tsv
