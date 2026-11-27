// Global unit-edit recurrence and forward word extension

use crate::word;
use rust_edlib::AlignError;

#[inline]
fn chunk<const N: usize, const REVERSE: bool>(s: &[u8]) -> Option<(&[u8; N], &[u8])> {
    if REVERSE {
        s.split_last_chunk().map(|(rest, last)| (last, rest))
    } else {
        s.split_first_chunk()
    }
}

/// Skip whole equal blocks. Inspect individual words only in the first
/// unequal block or in the trailing partial block.
#[inline]
fn common_end<const REVERSE: bool>(mut a: &[u8], mut b: &[u8]) -> usize {
    let mut matched = 0;
    while let (Some((x, rest_a)), Some((y, rest_b))) = (chunk::<32, REVERSE>(a), chunk::<32, REVERSE>(b)) {
        if x != y {
            return matched + matching_words::<REVERSE>(x, y);
        }
        matched += 32;
        (a, b) = (rest_a, rest_b);
    }
    matched + matching_words::<REVERSE>(a, b)
}

/// Count equal bytes a word at a time, working from the back for suffix searches.
#[inline]
pub(crate) fn matching_words<const REVERSE: bool>(a: &[u8], b: &[u8]) -> usize {
    #[inline]
    fn scan<'a, 'b, const REV: bool>(
        words: impl Iterator<Item = (&'a [u8], &'b [u8])>,
        tail: impl Iterator<Item = (&'a u8, &'b u8)>,
    ) -> usize {
        let mut matched = 0;
        for (x, y) in words {
            let xor = u64::from_le_bytes(x.try_into().unwrap())
                ^ u64::from_le_bytes(y.try_into().unwrap());
            if xor != 0 {
                let bits = if REV { xor.leading_zeros() } else { xor.trailing_zeros() };
                return matched + bits as usize / 8;
            }
            matched += 8;
        }
        matched + tail.take_while(|(x, y)| x == y).count()
    }
    let n = a.len().min(b.len());
    if REVERSE {
        let x = a[a.len() - n..].rchunks_exact(8);
        let y = b[b.len() - n..].rchunks_exact(8);
        let tail = x.remainder().iter().rev().zip(y.remainder().iter().rev());
        scan::<true>(x.zip(y), tail)
    } else {
        let x = a[..n].chunks_exact(8);
        let y = b[..n].chunks_exact(8);
        let tail = x.remainder().iter().zip(y.remainder());
        scan::<false>(x.zip(y), tail)
    }
}

pub(crate) fn trim<'a, 'b>(a: &'a [u8], b: &'b [u8]) -> (&'a [u8], &'b [u8]) {
    let prefix = common_end::<false>(a, b);
    let (a, b) = (&a[prefix..], &b[prefix..]);
    if a.is_empty() || b.is_empty() {
        return (a, b);
    }
    let suffix = common_end::<true>(a, b);
    (&a[..a.len() - suffix], &b[..b.len() - suffix])
}

/// Budget by shorter *trimmed* length.
pub(crate) fn cap(n: usize, m: usize) -> usize {
    match n.min(m) {
        0..=64 => 3,
        65..=1024 => 16,
        _ => 64,
    }
}

pub(crate) fn coordinates_fit(n: usize, m: usize) -> bool {
    n.saturating_add(m) < i32::MAX as usize
}

/// Unbounded dispatch on nonempty trimmed slices. The initial budget only
/// selects a scoring path; a projection never supplies the returned distance.
pub(crate) fn hybrid(a: &[u8], b: &[u8], k: usize) -> Result<usize, AlignError> {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if short.len() == 1 {
        return Ok(long.len() - usize::from(long.contains(&short[0])));
    }
    let gap = long.len() - short.len();
    if gap > k && long.len() <= isize::MAX as usize / 64 {
        if is_subsequence(short, long) {
            return Ok(gap);
        }
        // With at least as many insertions as query bytes, a few other edits
        // can still be cheaper to certify than a multiword Myers pass.
        if gap >= short.len() && short.len() > 64 {
            if let Some(extra) = near_subsequence(short, long) {
                return Ok(gap + extra);
            }
        }
    }
    if short.len() <= 64 {
        if gap <= k {
            if let Some(d) = within_three(short, long) {
                return Ok(d);
            }
        }
        return Ok(word::upper_bound(short, long));
    }
    // If one doubling could cover the gap, the first checkpoint is still a
    // useful density sample even when it cannot yet reach the endpoint.
    if gap <= 2 * k && coordinates_fit(a.len(), b.len()) {
        if let Ok(Some(d)) = bidirectional(a, b, k) {
            return Ok(d);
        }
    }
    rust_edlib::distance(a, b)
}

/// Both ends mismatch after trimming. A script of two or three edits must
/// edit each end, with at most one edit between the two resulting match runs.
/// All three possibilities for an end fit in registers; no wavefront rows.
fn within_three(a: &[u8], b: &[u8]) -> Option<usize> {
    fn reaches<const REV: bool>(a: &[u8], b: &[u8], gap: i32) -> [usize; 3] {
        let (x, y) = if REV {
            (&a[..a.len() - 1], &b[..b.len() - 1])
        } else {
            (&a[1..], &b[1..])
        };
        // A deletion needs gap + 1 insertions; a substitution still needs
        // gap insertions. Don't extend edits that cannot fit in three steps.
        let del = if gap <= 1 { matching_words::<REV>(x, b) } else { 0 };
        let sub = if gap <= 2 { 1 + matching_words::<REV>(x, y) } else { 0 };
        [del, sub, 1 + matching_words::<REV>(a, y)]
    }
    let gap = (b.len() - a.len()) as i32;
    let [ld, ls, li] = reaches::<false>(a, b, gap);
    let [rd, rs, ri] = reaches::<true>(a, b, gap);
    // Maximum target bytes covered by two ends with net diagonal z.
    let span = |z| match z {
        -1 => (ld + rs).max(ls + rd),
        0 => (ld + ri).max(ls + rs).max(li + rd),
        1 => (ls + ri).max(li + rs),
        2 => li + ri,
        _ => 0,
    };
    let two = span(gap);
    if two >= b.len() {
        return Some(2);
    }
    // A middle substitution or insertion consumes a target byte; deletion does not.
    (two.max(span(gap - 1)) + 1 >= b.len() || span(gap + 1) >= b.len()).then_some(3)
}

fn is_subsequence(short: &[u8], mut long: &[u8]) -> bool {
    let mut spare = long.len() - short.len();
    for &byte in short {
        let Some(skip) = long[..=spare].iter().position(|&b| b == byte) else {
            return false;
        };
        long = &long[skip + 1..];
        spare -= skip;
    }
    true
}

/// With the mandatory gap subtracted, insertions are free, substitutions cost
/// one and deletions cost two. Keep the earliest target end for each budget:
/// an earlier end can always skip forward, so it dominates every later end.
fn near_subsequence(short: &[u8], long: &[u8]) -> Option<usize> {
    let miss = long.len() + 1;
    let mut ends = [0; 3];
    for &byte in short {
        let [p0, p1, p2] = ends;
        let sub = |p: usize| (p + 1).min(miss);
        let advance = |p: usize, alternative: usize| {
            // A match must finish strictly before the substitution or deletion
            // to improve it. With no alternative, search to the target end.
            let stop = alternative.saturating_sub(1);
            if p < stop {
                long[p..stop].iter().position(|&b| b == byte).map_or(alternative, |i| p + i + 1)
            } else {
                alternative
            }
        };
        ends = [advance(p0, miss), advance(p1, sub(p0)), advance(p2, sub(p1).min(p0))];
        if ends[2] == miss {
            return None;
        }
    }
    ends.iter().position(|&p| p != miss)
}

const UNREACHABLE: i32 = -2; // Adding one to a sentinel must still be negative.
const INLINE_WIDTH: usize = 2 * 16 + 1;

struct Row {
    inline: [i32; INLINE_WIDTH + 4],
    spill: Vec<i32>,
}

impl Row {
    fn new() -> Self {
        Self { inline: [UNREACHABLE; INLINE_WIDTH + 4], spill: Vec::new() }
    }

    fn stored(&self, width: usize) -> &[i32] {
        if width <= INLINE_WIDTH { &self.inline } else { &self.spill }
    }

    // Two guards on either side for a one-diagonal change of the row bounds.
    // Live cells always start at 2, so only the right guards can become stale
    // when a narrowed row reuses backing that previously held a wider row.
    fn slots(&mut self, width: usize) -> Result<&mut [i32], AlignError> {
        let len = width + 4;
        if width <= INLINE_WIDTH {
            return Ok(&mut self.inline[..len]);
        }
        if len > self.spill.len() {
            self.spill.try_reserve(len - self.spill.len()).map_err(|_| AlignError::InputTooLong)?;
            self.spill.resize(len, UNREACHABLE);
        }
        Ok(&mut self.spill)
    }
}

// Insertion and substitution advance h; deletion advances only v.
#[inline]
fn next_offset(left: i32, center: i32, right: i32) -> i32 {
    right.max(left.max(center) + 1)
}

#[inline]
fn word_at<const REVERSE: bool>(bytes: &[u8], pos: i32) -> Option<u64> {
    // In reverse, the word ends at len - pos instead of starting at pos.
    // Negative starts (including words crossing the left edge) fail get.
    let start = if REVERSE { bytes.len() as i32 - pos - 8 } else { pos };
    let i = start as u32 as usize;
    Some(u64::from_le_bytes(bytes.get(i..i.checked_add(8)?)?.try_into().unwrap()))
}

#[inline]
fn advance<const ANCHORED: bool, const REVERSE: bool>(
    a: &[u8], b: &[u8], lo: i32, source: &[i32], dest: &mut [i32], mut upper: Option<&mut i32>,
) -> bool {
    // An edit plus a mismatch in the first word advances by at most eight
    // bytes. Only the slow path can exceed this bound on the whole row.
    if let Some(bound) = upper.as_deref_mut() {
        *bound = bound.saturating_add(8);
    }
    // A short global row costs less to seed and extend in one pass. On
    // larger rows the separate seed loop can vectorize across diagonals.
    if !ANCHORED && dest.len() <= 16 {
        for ((z, slot), p) in (lo..).zip(dest).zip(source.windows(3)) {
            *slot = next_offset(p[0], p[1], p[2]);
            extend_diagonal::<false, REVERSE>(a, b, z, slot, upper.as_deref_mut());
        }
        return false;
    }
    for (p, slot) in source.windows(3).zip(dest.iter_mut()) {
        *slot = next_offset(p[0], p[1], p[2]);
    }
    for (z, slot) in (lo..).zip(dest) {
        if extend_diagonal::<ANCHORED, REVERSE>(a, b, z, slot, upper.as_deref_mut()) {
            return true;
        }
    }
    false
}

// Called once per live diagonal; keep the seed in a register on the in-place path.
#[inline(always)]
fn extend_diagonal<const ANCHORED: bool, const REVERSE: bool>(
    a: &[u8], b: &[u8], z: i32, slot: &mut i32, upper: Option<&mut i32>,
) -> bool {
    let mut h = *slot;
    // A mismatch within the first word cannot reach the query end.
    if let (Some(x), Some(y)) = (word_at::<REVERSE>(a, h - z), word_at::<REVERSE>(b, h)) {
        let xor = x ^ y;
        if xor != 0 {
            let bits = if REVERSE { xor.leading_zeros() } else { xor.trailing_zeros() };
            *slot = h + (bits / 8) as i32;
            return false;
        }
        h += 8;
    }
    // Saturate a step just outside the rectangle. Removing one byte from both
    // prefixes cannot increase their edit distance, so the clipped cell is reachable.
    let edge = (b.len() as i32).min(a.len() as i32 + z);
    *slot = if h < 0 {
        UNREACHABLE
    } else if h >= edge {
        edge
    } else {
        let v = (h - z) as usize;
        let p = h as usize;
        let (x, y) = if REVERSE {
            (&a[..a.len() - v], &b[..b.len() - p])
        } else {
            (&a[v..], &b[p..])
        };
        let matched = if ANCHORED {
            matching_words::<REVERSE>(x, y)
        } else {
            common_end::<REVERSE>(x, y)
        };
        h + matched as i32
    };
    if let Some(bound) = upper {
        *bound = (*bound).max(*slot);
    }
    ANCHORED && *slot - z == a.len() as i32
}

struct Front {
    rows: [Row; 2],
    odd: bool,
    lo: i32,
    width: usize,
}

impl Front {
    fn new() -> Self {
        let mut rows = [Row::new(), Row::new()];
        rows[0].inline[2] = 0; // The caller already trimmed the matching end.
        Self { rows, odd: false, lo: 0, width: 1 }
    }

    fn offsets(&self) -> &[i32] {
        &self.rows[self.odd as usize].stored(self.width)[2..2 + self.width]
    }

    fn furthest(&self) -> i32 {
        *self.offsets().iter().max().unwrap()
    }

    #[inline]
    fn step<const ANCHORED: bool, const REVERSE: bool>(
        &mut self, a: &[u8], b: &[u8], lo: i32, hi: i32, upper: Option<&mut i32>,
    ) -> Result<bool, AlignError> {
        let [even, odd] = &mut self.rows;
        let (prev, curr) = if self.odd { (odd, even) } else { (even, odd) };
        let width = (hi - lo + 1) as usize;
        let first = (lo - self.lo + 1) as usize;
        let source = &prev.stored(self.width)[first..first + width + 2];
        let slots = curr.slots(width)?;
        // The same row storage held the frontier two scores ago.
        if width <= self.width {
            slots[2 + width..4 + width].fill(UNREACHABLE);
        }
        self.lo = lo;
        self.width = width;
        self.odd = !self.odd;
        Ok(advance::<ANCHORED, REVERSE>(a, b, lo, source, &mut slots[2..2 + width], upper))
    }

    #[inline]
    fn expand<const ANCHORED: bool, const REVERSE: bool>(
        &mut self, a: &[u8], b: &[u8], score: i32, upper: Option<&mut i32>,
    ) -> Result<bool, AlignError> {
        let lo = -score.min(a.len() as i32);
        let hi = score.min(b.len() as i32);
        self.step::<ANCHORED, REVERSE>(a, b, lo, hi, upper)
    }
}

/// Reverse diagonal r is forward diagonal `end - r`. Prefix distance
/// is nondecreasing along a diagonal: overlapping reaches on that diagonal
/// therefore certify a path costing at most the sum of the two budgets.
/// Conversely, splitting any optimal path at these edit budgets gives an
/// overlap. Testing successive budget sums therefore finds the minimum.
fn overlap(f: &Front, r: &Front, end: i32, m: i32) -> bool {
    let hi_f = f.lo + (f.width as i32 - 1);
    let hi_r = r.lo + (r.width as i32 - 1);
    let lo = f.lo.max(end - hi_r);
    let hi = hi_f.min(end - r.lo);
    if hi < lo {
        return false;
    }
    let forward = &f.offsets()[(lo - f.lo) as usize..=(hi - f.lo) as usize];
    let reverse = &r.offsets()[(end - hi - r.lo) as usize..=(end - lo - r.lo) as usize];
    forward.iter().zip(reverse.iter().rev()).any(|(&x, &y)| x >= m - y)
}

// On dense pairs, stop spending quadratic wavefront work and switch to
// bit-parallel Myers. Progress is sampled at doubling scores, never used
// as a certificate or as a position from which to restart.
const MIN_BYTES_PER_EDIT: u64 = 12;

/// A one-word band can still certify the answer if we haven't ruled it out.
/// Otherwise use frontier progress to choose a first Myers band, never a
/// restart position. Grow a missed band rather than retrying disproved widths.
fn finish_dense(a: &[u8], b: &[u8], score: u64, reached: u64) -> Result<Option<usize>, AlignError> {
    let (short, long) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    let gap = long.len() - short.len();
    let estimate = score * b.len() as u64 / reached.max(1);
    let can_certify = score < word::BAND as u64 && gap <= word::BAND;
    let first = if can_certify && estimate <= word::BAND as u64 {
        let upper = word::upper_bound(short, long);
        if upper <= word::BAND {
            return Ok(Some(upper));
        }
        upper as u64
    } else {
        // Skip a full one-word pass when it is unlikely to certify the score.
        // Shorter probes have less progress data; allow more headroom.
        // A miss still widens the band and retries from the origin.
        let headroom = if can_certify { estimate / 2 } else { estimate / 4 };
        estimate + headroom
    };
    // A displaced or very similar middle can make either hint much too
    // wide. Never speculate by more than a constant times the disproved budget.
    let hint = first.min(16 * score).min(long.len() as u64) as usize;
    if hint <= gap + 64 {
        // The ordinary score-only call already starts with at least this band.
        return Ok(None);
    }
    let mut slack = hint - gap;
    loop {
        if let Some(d) = rust_edlib::distance_at_most(short, long, gap + slack)? {
            return Ok(Some(d));
        }
        // At slack == short.len(), deleting/substituting/inserting always
        // fits the longer length, so the bounded call succeeds by then.
        slack = (2 * slack).min(short.len());
    }
}

fn bidirectional(a: &[u8], b: &[u8], mut checkpoint: usize) -> Result<Option<usize>, AlignError> {
    // Both matching ends have already been removed by trim.
    let (mut f, mut r) = (Front::new(), Front::new());
    let (mut left, mut right) = (0, 0);
    let m = b.len() as i32;
    let end = m - a.len() as i32;
    // At d <= BAND, projected wavefront work d²/2 can only exceed a
    // two-cells-per-byte scan for targets shorter than BAND²/4. Sample these
    // early for a one-word certificate, not for a speculative wide restart.
    let first_wide = checkpoint;
    if b.len() < word::BAND * word::BAND / 4 {
        checkpoint /= 2;
    }
    for score in 1..=a.len().max(b.len()) {
        if score % 2 == 1 {
            f.expand::<false, false>(a, b, (score / 2 + 1) as i32, Some(&mut left))?;
        } else {
            r.expand::<false, true>(a, b, (score / 2) as i32, Some(&mut right))?;
        }
        if left >= m - right && overlap(&f, &r, end, m) {
            return Ok(Some(score));
        }
        if score == checkpoint {
            checkpoint *= 2;
            let s = score as u64;
            let consider_wide = score >= first_wide;
            let reached = f.furthest() as u64 + r.furthest() as u64;
            if reached < s * MIN_BYTES_PER_EDIT
                && (consider_wide || s * m as u64 <= word::BAND as u64 * reached)
            {
                // On short pairs a few more wavefront rows can be cheaper
                // than even a linear word-band pass. Compare projected d²/2
                // with two cells per byte, where d ≈ s*m/reached.
                if score >= INLINE_WIDTH || s * s * m as u64 > 4 * reached * reached {
                    return finish_dense(a, b, s, reached);
                }
            }
        }
    }
    Ok(None)
}

pub(crate) enum AnchorProbe {
    Hit(usize),
    Miss,
    Myers,
}

const TIGHT_ANCHOR_WIDTH: usize = 2 * 32 + 1;

/// The cutoff guarantees that even the widest row fits on the stack.
fn anchored_inline<const REVERSE: bool, const WIDTH: usize>(
    a: &[u8], b: &[u8], limit: i32, end: i32,
) -> AnchorProbe {
    let mut row = [UNREACHABLE; WIDTH];
    row[0] = 0; // The pinned match was already trimmed.
    // At score s, lo = -s and hi = min(s, end + limit - s). The width
    // never shrinks, so newly exposed slots are still unreachable.
    let max_width = (end + limit + 1) as usize;
    for score in 1..=limit {
        let width = (2 * score as usize + 1).min(max_width);
        let (mut left, mut center) = (UNREACHABLE, UNREACHABLE);
        for (z, slot) in (-score..).zip(&mut row[..width]) {
            // lo moved left by one: this slot holds the old right neighbor.
            // Carry the two overwritten neighbors instead of building a second row.
            let right = *slot;
            *slot = next_offset(left, center, right);
            (left, center) = (center, right);
            if extend_diagonal::<true, REVERSE>(a, b, z, slot, None) {
                return AnchorProbe::Hit(score as usize);
            }
        }
    }
    AnchorProbe::Miss
}

/// Both remainders begin at a mismatch in traversal order and the one-sided
/// length gap fits k. Only exhausting the eligible diagonals proves a miss.
pub(crate) fn anchored<const REVERSE: bool>(
    a: &[u8], b: &[u8], k: usize, checkpoint: usize,
) -> AnchorProbe {
    if checkpoint == 0 || !coordinates_fit(a.len(), b.len()) {
        return AnchorProbe::Myers;
    }
    let limit = k.min(a.len());
    let end = b.len() as i32 - a.len() as i32;
    if limit <= INLINE_WIDTH / 2 {
        return anchored_inline::<REVERSE, INLINE_WIDTH>(a, b, limit as i32, end);
    }
    if limit <= TIGHT_ANCHOR_WIDTH / 2 {
        // Finishing scores s+1..=k visits at most (k-s)(k+s+2) diagonals.
        // Compare with a full Myers pass in 64-row blocks. This estimate
        // only chooses a scorer; the wavefront still has to prove hit or miss.
        let s = checkpoint.min(limit);
        let remaining = (limit - s) * (limit + s + 2);
        let blocks = a.len().div_ceil(u64::BITS as usize) as u64;
        if remaining as u64 <= blocks * b.len() as u64 {
            return anchored_inline::<REVERSE, TIGHT_ANCHOR_WIDTH>(a, b, limit as i32, end);
        }
    }
    let gap = a.len().saturating_sub(b.len());
    // For a tight deletion budget the strip can stay narrow even at high
    // scores. Otherwise skip a probe that cannot close the gap by its checkpoint.
    if gap > checkpoint && limit - gap >= INLINE_WIDTH {
        return AnchorProbe::Myers;
    }
    let mut through = checkpoint.min(limit) as i32;
    let limit = limit as i32;
    let mut front = Front::new();
    let mut first = 1;
    loop {
        // Only sample density between score intervals.
        for score in first..=through {
            // Every endpoint has z <= end; each remaining edit can reduce z by at most one.
            let hi = score.min(end + (limit - score));
            match front.step::<true, REVERSE>(a, b, -score, hi, None) {
                Ok(true) => return AnchorProbe::Hit(score as usize),
                Ok(false) => {}
                Err(_) => return AnchorProbe::Myers, // Myers has independent storage.
            }
        }
        if through == limit {
            return AnchorProbe::Miss;
        }
        if front.width >= INLINE_WIDTH && (front.furthest() as u64) < through as u64 * MIN_BYTES_PER_EDIT {
            return AnchorProbe::Myers;
        }
        first = through + 1;
        through = through.saturating_mul(2).min(limit);
    }
}

/// Global probe on remainders with both common ends removed.
pub(crate) fn distance_at_most(a: &[u8], b: &[u8], k: usize) -> Result<Option<usize>, AlignError> {
    let gap = a.len().abs_diff(b.len());
    if gap > k {
        return Ok(None);
    }
    if a.is_empty() || b.is_empty() {
        return Ok(Some(gap));
    }
    if !coordinates_fit(a.len(), b.len()) {
        return rust_edlib::distance_at_most(a, b, k);
    }
    let upper = a.len().max(b.len());
    if k < upper {
        bounded::<true>(a, b, k as i32)
    } else {
        bounded::<false>(a, b, upper as i32)
    }
}

fn bounded<const NARROW: bool>(a: &[u8], b: &[u8], limit: i32) -> Result<Option<usize>, AlignError> {
    let (n, m) = (a.len() as i32, b.len() as i32);
    let end = m - n;
    let mut front = Front::new();
    // Before this score, |z| <= score implies |z - end| <= limit - score.
    let full_through = if NARROW { (limit - end.abs()) / 2 } else { limit };
    let mut checkpoint: i32 = 64;
    for score in 1..=limit {
        let (mut lo, mut hi) = (-score.min(n), score.min(m));
        if NARROW && score > full_through {
            // Each remaining edit can change the diagonal by at most one.
            // Compute in i64: end +/- remaining can exceed i32 even when the
            // coordinates fit. Intersecting with the row before casting fits.
            let remaining = i64::from(limit - score);
            lo = i64::from(lo).max(i64::from(end) - remaining) as i32;
            hi = i64::from(hi).min(i64::from(end) + remaining) as i32;
        }
        front.step::<false, false>(a, b, lo, hi, None)?;
        if (lo..=hi).contains(&end) && front.offsets()[(end - lo) as usize] == m {
            return Ok(Some(score as usize));
        }
        if score == checkpoint {
            checkpoint = checkpoint.saturating_mul(2);
            // Keep tight cutoffs on the wavefront. For a wider cutoff, compare
            // the current diagonal work with a Myers pass over the whole band.
            let blocks = (limit.min(n).min(m) as u64).div_ceil(u64::BITS.into());
            if checkpoint <= limit
                && score as u64 * front.width as u64 >= blocks * n.max(m) as u64
                && (front.furthest() as u64) < score as u64 * MIN_BYTES_PER_EDIT
            {
                drop(front);
                // Restart from the origin with the full cutoff, not the unused edits.
                return rust_edlib::distance_at_most(a, b, limit as usize);
            }
        }
    }
    Ok(None)
}
