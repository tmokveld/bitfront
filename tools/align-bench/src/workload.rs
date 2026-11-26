use crate::Mode;

const DNA: &[u8; 4] = b"ACGT";
const GLOBAL_RATES: &[(usize, &[u32])] = &[
    (64, &[0, 100, 500, 2_000]),
    (256, &[0, 100, 500, 2_000]),
    (1_000, &[0, 100, 500, 2_000]),
    (10_000, &[0, 10, 100, 500, 1_000]),
    (100_000, &[0, 10, 100]),
    (1_000_000, &[0, 1]),
];
const BOUNDED_LENGTH: usize = 256;
const BOUNDED_RATE_BP: u32 = 1_250;
const BOUNDED_EDITS: usize = edit_count(BOUNDED_LENGTH, BOUNDED_RATE_BP);

pub(crate) struct Case {
    pub(crate) mode: Mode,
    pub(crate) length: usize,
    pub(crate) rate_bp: u32,
    pub(crate) flank: usize,
}

pub(crate) fn cases(smoke: bool) -> Vec<Case> {
    let mut cases = Vec::new();
    let mut add = |mode, length, rate_bp, flank| {
        cases.push(Case {
            mode,
            length,
            rate_bp,
            flank,
        });
    };
    for &(length, rates) in GLOBAL_RATES {
        for &rate_bp in rates {
            if !smoke
                || matches!(
                    (length, rate_bp),
                    (64, 0 | 100 | 2_000) | (1_000, 500) | (1_000_000, 1)
                )
            {
                add(Mode::Global, length, rate_bp, 0);
            }
        }
    }

    // Only the query contains N. Each of its N bytes costs exactly one edit.
    for k in [BOUNDED_EDITS, BOUNDED_EDITS - 1] {
        add(Mode::AtMost(k), BOUNDED_LENGTH, BOUNDED_RATE_BP, 0);
    }

    // Also score one identical flanked panel both globally and with the tail free.
    for (length, rate_bp) in [(256, 500), (10_000, 100)] {
        if !smoke || length == 256 {
            if length == 256 {
                add(Mode::Global, length, rate_bp, 3 * length);
            }
            for mode in [Mode::Prefix, Mode::Suffix] {
                add(mode, length, rate_bp, 3 * length);
            }
        }
    }
    cases
}

impl Case {
    pub(crate) fn pair(&self, index: usize) -> (Vec<u8>, Vec<u8>) {
        // Rate and mode are excluded so rows of the same length share a source.
        let seed = 0x243f_6a88_85a3_08d3 ^ ((self.length as u64) << 32) ^ index as u64;
        let mut seeds = SplitMix64(seed);
        let mut source_rng = SplitMix64(seeds.next());
        let mut flank_rng = SplitMix64(seeds.next());
        let mut edit_rng = SplitMix64(seeds.next());

        let leading = if self.mode == Mode::Suffix {
            self.flank
        } else {
            0
        };
        let trailing = self.flank - leading;
        let source_end = leading + self.length;
        let mut target = Vec::with_capacity(source_end + trailing);
        target.extend((0..leading).map(|_| flank_rng.base()));
        target.extend((0..self.length).map(|_| source_rng.base()));
        target.extend((0..trailing).map(|_| flank_rng.base()));
        let source = &target[leading..source_end];

        let edits = edit_count(self.length, self.rate_bp);
        if edits == 0 {
            return (source.to_vec(), target);
        }

        let bounded = matches!(self.mode, Mode::AtMost(_));
        let mut edits_left = edits;
        let mut query = Vec::with_capacity(source.len() + edits);
        for (i, &base) in source.iter().enumerate() {
            // Sample exactly `edits` distinct sites without a position buffer.
            if edits_left != 0 && edit_rng.below(source.len() - i) < edits_left {
                edits_left -= 1;
                if bounded {
                    query.push(b'N');
                    continue;
                }
                match edit_rng.below(3) {
                    0 => query.push(DNA[(base_index(base) + 1 + edit_rng.below(3)) % DNA.len()]),
                    1 => {
                        query.push(edit_rng.base());
                        query.push(base);
                    }
                    _ => {} // Delete the source base.
                }
            } else {
                query.push(base);
            }
        }
        (query, target)
    }
}

const fn edit_count(length: usize, rate_bp: u32) -> usize {
    if rate_bp == 0 {
        return 0;
    }
    let rounded = ((length as u64 * rate_bp as u64 + 5_000) / 10_000) as usize;
    if rounded == 0 { 1 } else { rounded }
}

fn base_index(base: u8) -> usize {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => unreachable!(),
    }
}

struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn below(&mut self, limit: usize) -> usize {
        (self.next() % limit as u64) as usize
    }

    fn base(&mut self) -> u8 {
        DNA[self.below(DNA.len())]
    }
}
