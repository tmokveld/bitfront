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
pub(crate) const HELD_OUT: [usize; 4] = [0x28ba, 0x7319, 0xab71, 0xed09];
const CONTROL_LENGTH: usize = 10_000;
const SPARSE_EDITS: usize = 8;

#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum Scenario {
    Seeded(Option<usize>),
    MixedBelow,
    MixedAt,
    SparseLoose,
}

impl Scenario {
    pub(crate) fn name(self) -> Option<&'static str> {
        match self {
            Self::Seeded(_) => None,
            Self::MixedBelow => Some("mixed-below"),
            Self::MixedAt => Some("mixed-at"),
            Self::SparseLoose => Some("sparse-loose"),
        }
    }

    pub(crate) fn cutoff(self, distance: usize) -> Option<usize> {
        match self {
            Self::Seeded(k) => k,
            Self::MixedBelow => Some(distance.checked_sub(1).expect("boundary miss needs a nonzero distance")),
            Self::MixedAt => Some(distance),
            Self::SparseLoose => Some(4_096),
        }
    }
}

pub(crate) const MIXED_AT: Case = Case::control(Scenario::MixedAt, 1000);
pub(crate) const SPARSE_LOOSE: Case =
    Case::control(Scenario::SparseLoose, (SPARSE_EDITS * 10_000 / CONTROL_LENGTH) as u32);
pub(crate) const REGRESSIONS: [Case; 3] = [
    Case::control(Scenario::MixedBelow, 1000),
    MIXED_AT,
    SPARSE_LOOSE,
];

pub(crate) struct Case {
    pub(crate) mode: Mode,
    pub(crate) length: usize,
    pub(crate) rate_bp: u32,
    pub(crate) flank: usize,
    pub(crate) scenario: Scenario,
    pub(crate) swap: bool,
}

pub(crate) fn cases(smoke: bool) -> Vec<Case> {
    let mut cases = Vec::new();
    for &(length, rates) in GLOBAL_RATES {
        for &rate_bp in rates {
            if !smoke
                || matches!(
                    (length, rate_bp),
                    (64, 0 | 100 | 2_000) | (1_000, 500) | (1_000_000, 1)
                )
            {
                cases.push(Case::seeded(Mode::Global, length, rate_bp, 0, None));
            }
        }
    }

    // Only the query contains N. Each of its N bytes costs exactly one edit.
    for mode in [Mode::Global, Mode::Prefix, Mode::Suffix] {
        let flank = if mode == Mode::Global { 0 } else { 3 * BOUNDED_LENGTH };
        for k in [BOUNDED_EDITS, BOUNDED_EDITS - 1] {
            cases.push(Case::seeded(mode, BOUNDED_LENGTH, BOUNDED_RATE_BP, flank, Some(k)));
        }
    }
    if !smoke {
        let (length, rate_bp) = (10_000, 1_000);
        let edits = edit_count(length, rate_bp);
        for k in [edits - 1, edits] {
            cases.push(Case::seeded(Mode::Global, length, rate_bp, 0, Some(k)));
        }
        cases.extend(REGRESSIONS);
    }

    // Also score one identical flanked panel both globally and with the tail free.
    for (length, rate_bp) in [(256, 500), (10_000, 100)] {
        if !smoke || length == 256 {
            if length == 256 {
                cases.push(Case::seeded(Mode::Global, length, rate_bp, 3 * length, None));
            }
            for mode in [Mode::Prefix, Mode::Suffix] {
                cases.push(Case::seeded(mode, length, rate_bp, 3 * length, None));
            }
        }
    }
    cases
}

impl Case {
    pub(crate) const fn seeded(mode: Mode, length: usize, rate_bp: u32, flank: usize, cutoff: Option<usize>) -> Self {
        Self { mode, length, rate_bp, flank, scenario: Scenario::Seeded(cutoff), swap: false }
    }

    const fn control(scenario: Scenario, rate_bp: u32) -> Self {
        Self { mode: Mode::Global, length: CONTROL_LENGTH, rate_bp, flank: 0, scenario, swap: false }
    }

    pub(crate) fn pair(&self, index: usize) -> (Vec<u8>, Vec<u8>) {
        let (index, marked) = match self.scenario {
            Scenario::Seeded(cutoff) => (index, cutoff.is_some()),
            Scenario::MixedBelow | Scenario::MixedAt => (HELD_OUT[index], false),
            Scenario::SparseLoose => {
                let mut rng = SplitMix64::held_out(index);
                let mut target = rng.sequence(self.length, DNA);
                let mut query = target.clone();
                // Eight T→A substitutions, including both ends; both bytes also occur unchanged.
                for j in 0..SPARSE_EDITS {
                    let site = j * (self.length - 1) / (SPARSE_EDITS - 1);
                    target[site] = b'T';
                    query[site] = b'A';
                }
                return (query, target);
            }
        };
        // Rate and mode are excluded so rows of the same length share a source.
        let seed = 0x243f_6a88_85a3_08d3 ^ ((self.length as u64) << 32) ^ index as u64;
        let mut seeds = SplitMix64(seed);
        let mut source_rng = SplitMix64(seeds.next());
        let mut flank_rng = SplitMix64(seeds.next());
        let edit_rng = SplitMix64(seeds.next());

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
        let query = self.edit_query(&target[leading..source_end], edit_rng, marked);
        if self.swap { (target, query) } else { (query, target) }
    }

    fn edit_query(&self, source: &[u8], mut edit_rng: SplitMix64, marked: bool) -> Vec<u8> {
        let edits = edit_count(self.length, self.rate_bp);
        if edits == 0 {
            return source.to_vec();
        }

        let mut edits_left = edits;
        let capacity = if marked { source.len() } else { source.len() + edits };
        let mut query = Vec::with_capacity(capacity);
        for (i, &base) in source.iter().enumerate() {
            // Sample exactly `edits` distinct sites without a position buffer.
            if edits_left != 0 && edit_rng.below(source.len() - i) < edits_left {
                edits_left -= 1;
                if marked {
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
        query
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

pub(crate) struct SplitMix64(pub(crate) u64);

impl SplitMix64 {
    pub(crate) fn held_out(index: usize) -> Self {
        Self(0x4e6a_1d35_0cb7_84f0 ^ HELD_OUT[index] as u64)
    }

    pub(crate) fn sequence(&mut self, length: usize, alphabet: &[u8]) -> Vec<u8> {
        (0..length).map(|_| alphabet[self.below(alphabet.len())]).collect()
    }

    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    pub(crate) fn below(&mut self, limit: usize) -> usize {
        (self.next() % limit as u64) as usize
    }

    fn base(&mut self) -> u8 {
        DNA[self.below(DNA.len())]
    }
}
