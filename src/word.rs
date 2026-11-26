//! Score a global alignment in a moving window of 64 query rows. Every
//! boundary value represents a real edit script, so the result is an upper
//! bound even when the optimal path leaves the window.

const WIDTH: usize = u64::BITS as usize;
const STEP: usize = 8;
pub(super) const BAND: usize = WIDTH - STEP;

struct Column {
    positive: u64,
    negative: u64,
}

impl Column {
    const INITIAL: Self = Self { positive: u64::MAX, negative: 0 };

    #[inline]
    fn advance(&mut self, eq: u64) {
        // Myers vertical differences. The top boundary advances by insertion.
        let (pv, mv) = (self.positive, self.negative);
        let xv = eq | mv;
        let sum = (eq & pv).wrapping_add(pv);
        // Keep the shifted horizontal positive mask complemented: its low
        // bit is zero, so inserting the top boundary needs no extra OR.
        let not_ph = ((sum | eq | pv) & !mv) << 1;
        let mh = pv & (eq | !sum);
        self.positive = (mh << 1) | (not_ph & !xv);
        self.negative = xv & !not_ph;
    }

    #[inline]
    fn value(&self, top: usize, rows: usize) -> usize {
        let used = u64::MAX.checked_shr((WIDTH - rows) as u32).unwrap_or(0);
        top + (self.positive & used).count_ones() as usize
            - (self.negative & used).count_ones() as usize
    }
}

/// Every global path crosses the middle target column at some query row i.
/// Minimize prefix(i) + suffix(rows-i). The two half-column recurrences
/// are independent, cutting the serial dependency through the target in half.
fn meet(masks: &[u64; 256], rows: usize, target: &[u8]) -> usize {
    let (prefix, suffix) = target.split_at(target.len() / 2);
    let reverse = |bits: u64| bits.reverse_bits() >> (WIDTH - rows);
    let (mut left, mut right) = (Column::INITIAL, Column::INITIAL);
    let mut back = suffix.iter().rev();
    for (&a, &b) in prefix.iter().zip(back.by_ref()) {
        left.advance(masks[a as usize]);
        right.advance(reverse(masks[b as usize]));
    }
    if let Some(&b) = back.next() {
        right.advance(reverse(masks[b as usize]));
    }

    let cost = |i| left.value(prefix.len(), i) + right.value(suffix.len(), rows - i);
    let gap = target.len() - rows;
    let slack = (cost(rows / 2) - gap) / 2;
    // The midpoint split supplies a feasible script. Length gaps alone
    // exclude rows outside this interval from improving it.
    let lo = prefix.len().saturating_sub(gap + slack);
    let hi = rows.min(prefix.len() + slack);
    let mut score = cost(lo);
    let mut best = score;
    let (mut lp, mut lm) = (left.positive >> lo, left.negative >> lo);
    let (mut rp, mut rm) = (reverse(right.positive) >> lo, reverse(right.negative) >> lo);
    for _ in lo..hi {
        // Moving the split down adds a left row and removes a right row.
        score = score + ((lp & 1) + (rm & 1)) as usize - ((lm & 1) + (rp & 1)) as usize;
        best = best.min(score);
        (lp, lm, rp, rm) = (lp >> 1, lm >> 1, rp >> 1, rm >> 1);
    }
    best
}

/// Exact for a one-word query. Otherwise `gap <= BAND` is required, and
/// a result at most `BAND` is exact: on diagonal z = j-i, any path costs at
/// least |z| + |gap-z| to pass through it and still reach the endpoint.
pub(super) fn upper_bound(query: &[u8], target: &[u8]) -> usize {
    debug_assert!(!query.is_empty() && query.len() <= target.len());
    let mut masks = [0u64; 256];
    for (i, &byte) in query.iter().take(WIDTH).enumerate() {
        masks[byte as usize] |= 1 << i;
    }
    let mut col = Column::INITIAL;
    if query.len() <= WIDTH {
        return meet(&masks, query.len(), target);
    }

    let gap = target.len() - query.len();
    debug_assert!(gap <= BAND);
    // Finish each tile when the lowest viable diagonal reaches its bottom
    // row. After the slide, its top boundary is above every viable row of
    // the next column. No path costing <= BAND needs the discarded rows.
    let below = (BAND - gap) / 2;
    let (mut base, mut start, mut top) = (0, 0, 0);
    loop {
        let remaining = query.len() - base;
        let end = if remaining == WIDTH { target.len() } else { base + WIDTH - below };
        let rotation = (base % WIDTH) as u32;
        for &byte in &target[start..end] {
            col.advance(masks[byte as usize].rotate_right(rotation));
        }
        top += end - start;
        if remaining == WIDTH {
            return col.value(top, WIDTH);
        }
        let shift = STEP.min(remaining - WIDTH);
        top = col.value(top, shift);
        // Carry the overlapping rows forward; initialize new bottom rows by
        // deletion from the old bottom. Mask slots wrap in absolute query order.
        col.positive = !(!col.positive >> shift);
        col.negative >>= shift;
        let mut bit = 1u64.rotate_left(rotation);
        for (&old, &new) in query[base..base + shift].iter()
            .zip(&query[base + WIDTH..base + WIDTH + shift])
        {
            masks[old as usize] &= !bit;
            masks[new as usize] |= bit;
            bit = bit.rotate_left(1);
        }
        base += shift;
        start = end;
    }
}
