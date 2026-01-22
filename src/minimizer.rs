use packed_seq::AsciiSeq;
use simd_minimizers::seq_hash::NtHasher;

use crate::cfg::MinimizerCfg;

/// Return (positions, minimizer_values) using simd_minimizers Builder API.
///
/// - forward_only=true  -> NtHasher::<false>
/// - forward_only=false -> canonical minimizers (strand-agnostic) using NtHasher::<true>
pub fn sampled_minimizers(seq: &[u8], mm: &MinimizerCfg) -> (Vec<u32>, Vec<u64>) {
    let k = mm.k;
    let w = mm.w;

    let mut pos: Vec<u32> = Vec::new();

    if mm.forward_only {
        let hasher = NtHasher::<false>::new(k);
        let vals: Vec<u64> = simd_minimizers::minimizers(k, w)
            .hasher(&hasher)
            .run(AsciiSeq(seq), &mut pos)
            .values_u64()
            .collect();
        (pos, vals)
    } else {
        // Canonical minimizers; still return values from simd_minimizers so matching is consistent.
        let hasher = NtHasher::<true>::new(k);
        let vals: Vec<u64> = simd_minimizers::canonical_minimizers(k, w)
            .hasher(&hasher)
            .run(AsciiSeq(seq), &mut pos)
            .values_u64()
            .collect();
        (pos, vals)
    }
}
