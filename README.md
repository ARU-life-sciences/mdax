# **Fo**ldback **s**plitter

`fos` is a way to detect simple inverted duplication chimeras from 'Multiple Displacement Amplification' based whole genome amplification. Which is a common problem I hear.

Sensitive mode:
/target/release/mdax --output test ./test_pal/out/subsample/SRR24201687.first20k_2.fasta --min-span 800 --min-matches 10 --end-guard 200 --forward-only false > test_sensitive.tsv
