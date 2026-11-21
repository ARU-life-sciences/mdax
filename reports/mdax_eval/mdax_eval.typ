#import "@local/report:0.1.0": *
#show: report.with(
  title: "Is mdax worth using on MDA data?",
  subtitle: "Evaluation of mdax on C. elegans, A. besseyi, Lepidodermella and twelve other nematodes, 18–21 November 2026",
  abstract: [
    mdax removes foldback chimeras, the inverted copies that multiple displacement amplification (MDA) leaves in long reads.
    We scored the current version (main at `812d743`) against a reference-based truth on four read sets from two
    species, and assembled corrected reads from 11 C. elegans HiFi libraries and 16 ONT or PacBio CLR samples from 14 species. On MDA HiFi it cuts
    74–79% of true foldbacks while cutting about 2% of normal reads. The 95–99% of foldbacks whose inverted copy is
    still in the read are nearly all cut; the rest have no copy to find. Assembly gains are large on MDA HiFi: against
    raw reads, contigs fall by 40–62%, N50 rises 2–4-fold and BUSCO duplication drops from 6–10% to 1–4%. On MDA
    ONT and PacBio CLR the gains are real but smaller, and mdax matches or beats the palindrome removal of Lee et al. on
    every sample with a decent assembly. On libraries that were not amplified it gives nothing and costs a little. The
    answer is yes for MDA data; no for anything else.
  ],
  toc: false,
)

= Verdict

#solid *Use mdax on MDA HiFi.* On all six amplified C. elegans HiFi libraries assembled with the new code, the hifiasm
assembly from mdax-corrected reads has 40–62% fewer contigs than raw reads, two to four times the N50, higher BUSCO
completeness, and BUSCO duplication falls from 6–10% to 0.8–4.2%. The old (October) mdax cut contigs by only 17–28%.

#solid *Use it on MDA ONT and PacBio CLR from a single animal, with smaller expectations.* On C. elegans ONT, A. besseyi
and the first Lepidodermella sample, mdax beats raw reads clearly and beats or matches Lee et al.'s
`remove_palindrome.sh`. The gain over raw is mostly in contig count and N50; BUSCO completeness hardly moves.

#caveat *It does not rescue a poor assembly.* On the ten non-model nematodes whose assemblies are in tens of thousands of
contigs with BUSCO 20–50%, mdax and Lee et al.'s method tie, and both beat raw by a similar margin. Something other than
foldbacks limits those assemblies (coverage, heterozygosity, genome size).

#null *Do not use it on libraries that were not amplified.* On three standard C. elegans HiFi libraries the raw reads
assemble marginally better (N50 6–9% higher). mdax is safe to run there (it now recognises such libraries and keeps
its threshold low) but there is nothing for it to fix.

= What was tested

*Read-level truth.* For a random sample of about 40,000 reads from each library we aligned the raw reads to an
independent reference with minimap2 and called a read a true foldback when it had two aligned segments on opposite
strands of the same sequence, within 5 kb. That gave four truth sets: nrCaeEleg92 (MDA HiFi, ~140×) and nrCaeEleg94
(MDA HiFi, ~1500×) against WBcel235; celegans_wga (Lee et al.'s MDA ONT, SRR24201716) against WBcel235; and A. besseyi
strain APTV (MDA ONT, SRR24201688) against an independent assembly of strain AORJ (GCA_024699835.1) made from
unamplified reads. mdax runs on the full library, since its decision to keep a foldback as genuine depends on how many
other reads share it, and the sample reads are then scored. Every number below for "linear reads cut" is the share of
reads with a clean single alignment that mdax cut anyway.

*Assembly arms.* HiFi reads went through mdax (balanced mode), hifiasm 0.25 and BUSCO (nematoda_odb10). ONT and CLR reads
went through mdax (permissive mode, ONT refinement), Flye 2.9 (`--meta`) and BUSCO. Each arm is compared with the same
reads uncorrected and, for the ONT and CLR samples, with the palindrome removal from Lee et al. (two rounds of
`remove_palindrome.sh`, from the paper's Zenodo pipeline). All arms used identical downstream settings.

= Read level: is mdax finding the real chimeras?

#figure(tsv-table(read("results/read_level.tsv")),
  caption: [Outcome of mdax on reference-labelled reads. "cut" means mdax split the read; "kept as real" means it
  detected the foldback but judged it a genuine palindrome; "missed" means it found nothing.]) <tab-reads>

The old October build cut only 16–38% of true foldbacks on HiFi, because it judged 35–52% of them to be genuine
palindromes: deep MDA libraries contain many artefacts at the same locus by chance, and the old support count took that
recurrence as evidence of real structure. The fold-axis gate fixed this (@tab-reads), and the rescue pass added since
then lifts ONT by 11–12 points. The cost is a small rise in normal reads cut, to 1.8–2.3% in C. elegans and 4.6% in
A. besseyi. In A. besseyi, 97% of those "false" cuts fall at a read end the aligner could not place, which is what a
short-armed foldback looks like to the truth; only 33 of 26,428 normal reads (0.12%) were cut with no such unplaced
end.

#note(title: "Why 20–25% of foldbacks are never found")[
  Nearly all of the missed foldbacks are not misses. When the polymerase switches strand it lands some distance _d_
  from where it left off (median 2.6 kb in these cases), and the inverted copy only begins once the new strand has
  travelled back over that gap. In 77–94% of the no-copy reads, the read ends before that happens, so the read
  holds sequence X followed by the reverse complement of a different, nearby sequence Y, and no part of it is the
  reverse complement of another part. Among foldbacks that do carry an inverted copy of at least 100 bp, mdax cuts
  98.8% (nrCaeEleg92), 95.2% (celegans_wga) and 98.3% (A. besseyi). The no-copy class is 19–21% of true foldbacks in
  all three sets, independent of platform and species. Their junctions are unique to one read (99.5–99.9%), so they
  are artefacts, not real inversions.
]

= Assembly: MDA HiFi (hifiasm)

#figure(tsv-table(read("results/hifi_hifiasm.tsv")),
  caption: [hifiasm primary contigs. "v5" is the fold-axis gate alone (18 Nov); "main" adds alignment-verified rescue
  and the corrected cut placement. Samples 97, 99 and 100 were not re-run with main; on 92, 96 and 98, main gained
  0.10–0.20 Mb of N50 over v5.]) <tab-hifi>

The gains are large and consistent (@tab-hifi). On the six MDA libraries, mdax roughly halves the
contig count against raw reads (92: 1248 → 711; 96: 832 → 319; 98: 523 → 199) and doubles to quadruples N50. BUSCO
duplication, which in raw MDA assemblies is inflated by chimeric reads building duplicated contigs, falls from 6.5–9.9%
to 0.8–1.0% on 96 and 98. nrCaeEleg92 is the one exception: its duplication is higher with main (4.2%) than with v5
(2.9%), which we have not explained, though it is still well under raw (7.2%). The total assembly length falls from
about 111 Mb to 103–110 Mb, against a 100 Mb reference, which is the expected direction.

On the three libraries that were not amplified (34, ToL, 33), mdax does nothing useful: N50 is within 0.04 Mb of raw and
BUSCO is unchanged. Version 5 briefly made these worse (34: 638 contigs against 592 raw) because its calibration pushed
the threshold to 47–60; the fix in `643e0bf` restored them to the old build's numbers.

The three deepest C. elegans libraries (93, 94, 95, ~1000–1500×) were not assembled with the new code: the hugemem
re-runs never completed. Their old assemblies are poor in every arm (9,000–12,000 contigs, BUSCO duplication about
50%), so the question for them is coverage, not correction.

= Assembly: MDA ONT and PacBio CLR (Flye)

#figure(tsv-table(read("results/ont_flye.tsv")),
  caption: [Flye contigs for the samples re-run with the new code. "Lee et al." is two rounds of their
  `remove_palindrome.sh`. Lepidodermella (gastrotrich) BUSCO is metazoa_odb10 (n = 954).]) <tab-ont>

On C. elegans ONT, mdax main gives 563 contigs and N50 0.65 Mb, against 748 / 0.51 raw and 622 / 0.55 for Lee et al.
(@tab-ont). On A. besseyi it gives 128 contigs and N50 2.00 Mb, against 215 / 1.31 raw and 141 / 1.25 for Lee et al.
On the first Lepidodermella sample (SRR27483171) the gain over raw is dramatic (3960 → 190 contigs, BUSCO duplication
10.1% → 1.4%) and mdax edges Lee et al. (203 contigs, N50 1.15 Mb). On the second (SRR27483172) mdax has fewer contigs
than Lee et al. (2788 against 3809) but higher BUSCO duplication (13.7% against 11.1%), so neither method clearly wins.
On Sabatieria and Microlaimidae no method helps much. The old October mdax had blown up the Microlaimidae assembly to
1084 Mb; the fold-axis gate fixed that.

#figure(tsv-table(read("results/ont_flye_suite.tsv")),
  caption: [The remaining Lee et al. nematodes, run with the old October mdax only (N50 in Mb, BUSCO C in %).
  Missing cells are arms that did not complete.]) <tab-suite>

On the ten non-model nematodes run only with the old build (@tab-suite), mdax and Lee et al.'s method are level: mdax has
fewer contigs on six of the nine species with both arms, and the two methods' N50 agree to within 0.03 Mb on all of them. Both reduce contig count against raw by 2–46%. BUSCO completeness is within 1.5 points of raw throughout. These assemblies are so fragmented that the foldback rate is not what limits them.

= Changes that were tried and rejected

Several changes improved the read-level score but made assemblies worse, and are kept off main. A coverage-normalised
gate (fraction of reads covering a position that fold there) gave no gain on deep HiFi, was worse on ONT and added 30–200% to the runtime. Recursing into second foldbacks in reads with several cut the leftover foldbacks on ONT from 66% to 31%
of multi-junction reads, but cost hifiasm and Flye contiguity (92: 711 → 741 contigs; celegans_wga: N50 0.65 → 0.51).

The clearest test was an oracle: we took mdax's output and cut it again at every foldback junction the reference could
still see, which is the most any further detection could achieve (@tab-oracle). It did not help: N50 was unchanged on
both samples, and contigs and duplication rose on nrCaeEleg92. Keeping both pieces of each cut made it worse still.
The leftover junctions are unique to one read and carry no repeat, so hifiasm already discards them; cutting them only
shortens reads and throws away unique sequence.

#figure(tsv-table(read("results/hifi_experiments.tsv")),
  caption: [hifiasm on mdax main output with further cuts: the rejected recursion change, and three oracle arms that cut at
  reference-detected junctions.]) <tab-oracle>

= Limits of this evaluation

#caveat Each configuration was assembled once, except one arm on nrCaeEleg92 that was repeated and reproduced to
within one contig; differences of a few percent in contig count are within run-to-run noise. Only one assembler per
platform was tested. The read-level truth comes from two species, both nematodes with compact genomes; the A. besseyi
reference is a different strain, which may count a few real inversions as foldbacks. No polyploid or repeat-rich
genome has been tested. The oracle used hifiasm only, so the conclusion that leftover foldbacks don't matter is not
proven for Flye. Ten of the ONT nematodes and three C. elegans HiFi libraries have only been assembled with the old
build.

= Practical guidance

Run mdax on MDA reads before assembly: `mdax reads.fa -o corrected.fa -r report.tsv` (balanced mode is the default for
HiFi; use `--mode permissive --refine-mode ont` for ONT). A full library takes 3–20 minutes on 16 threads and under 4 GB
of memory. Check the `[QC SUMMARY]` line: on an amplified library the support threshold is calibrated upward (7–14 on the C. elegans HiFi libraries); on a library that was not amplified calibration is skipped and the threshold stays at 3,
which is the sign that there was little to correct. Don't expect it to fix low coverage, heterozygosity or an assembly
that is fragmented for other reasons.
