use packed_seq::AsciiSeq;
use simd_minimizers::seq_hash::NtHasher;

use crate::cfg::MinimizerCfg;

// Return (positions, minimizer_values) using simd_minimizers Builder API.
//
// - forward_only=true  -> NtHasher::<false>
// - forward_only=false -> canonical minimizers (strand-agnostic) using NtHasher::<true>

/// Fill (pos_out, val_out) with minimizers for `seq`, reusing allocations.
///
/// IMPORTANT:
/// - Clears both output vecs.
/// - Preserves capacity across calls.
///
pub fn sampled_minimizers_into(
    seq: &[u8],
    mm: &MinimizerCfg,
    pos_out: &mut Vec<u32>,
    val_out: &mut Vec<u64>,
) {
    pos_out.clear();
    val_out.clear();

    let k = mm.k;
    let w = mm.w;

    if mm.forward_only {
        let hasher = NtHasher::<false>::new(k);
        let runner = simd_minimizers::minimizers(k, w)
            .hasher(&hasher)
            .run(AsciiSeq(seq), pos_out);
        let iter = runner.values_u64();
        val_out.extend(iter);
    } else {
        let hasher = NtHasher::<true>::new(k);
        let runner = simd_minimizers::canonical_minimizers(k, w)
            .hasher(&hasher)
            .run(AsciiSeq(seq), pos_out);
        let iter = runner.values_u64();
        val_out.extend(iter);
    }
}
