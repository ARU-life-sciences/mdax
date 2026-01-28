use gxhash::HashMapExt;

use crate::foldback::BinStat;

pub struct FoldScratch {
    // ...
    pub rc: Vec<u8>,
    // minimizer buffers
    pub pos_f: Vec<u32>,
    pub val_f: Vec<u64>,
    pub pos_rc: Vec<u32>,
    pub val_rc: Vec<u64>,

    // ...
    pub idx_f: gxhash::HashMap<u64, Vec<i32>>,
    pub repetitive: gxhash::HashMap<u64, ()>,
    pub stats: gxhash::HashMap<i32, BinStat>,
    pub best_pts: Vec<(i32, i32)>,
}

impl FoldScratch {
    pub fn new() -> Self {
        Self {
            rc: Vec::new(),
            pos_f: Vec::new(),
            val_f: Vec::new(),
            pos_rc: Vec::new(),
            val_rc: Vec::new(),
            idx_f: gxhash::HashMap::new(),
            repetitive: gxhash::HashMap::new(),
            stats: gxhash::HashMap::new(),
            best_pts: Vec::new(),
        }
    }

    #[inline]
    pub fn clear(&mut self) {
        self.idx_f.clear();
        self.repetitive.clear();
        self.stats.clear();
        self.best_pts.clear();
    }
}

#[derive(Debug, Default)]
pub struct SigScratch {
    pub right_rc: Vec<u8>,

    // minimizer output buffers
    pub pos_l: Vec<u32>,
    pub val_l: Vec<u64>,
    pub pos_r: Vec<u32>,
    pub val_r: Vec<u64>,

    // combined values (for hashing)
    pub combined: Vec<u64>,
}
