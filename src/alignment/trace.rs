use std::mem::size_of;

use crate::kernel::{MIN_BYTES_PER_EDIT, UNREACHABLE, advance, coordinates_fit};

use super::{Alignment, Error, Frame, HISTORY_BYTES, Op, Run, push_run, trace_edit};

pub(super) enum Probe {
    Hit(Alignment),
    Miss,
    Switch,
}

const GUARDS: usize = 4;
const INLINE_SCORE: i32 = 32;
// Rows 0..=s have at most (s+1)^2 live cells, plus four guards per row.
const INLINE_REACHES: usize = (INLINE_SCORE as usize + 1) * (INLINE_SCORE as usize + 1 + GUARDS);
const SPILL_REACHES: usize = HISTORY_BYTES / size_of::<i32>() - INLINE_REACHES;
const GLOBAL_LEAF_LIMIT: i32 = 16;

#[derive(Clone, Copy)]
struct Shape {
    n: i32,
    m: i32,
    limit: i32,
}

impl Shape {
    #[inline]
    fn band<const ANCHORED: bool>(self, score: i32) -> (i32, i32) {
        let end = i64::from(self.m - self.n);
        let remaining = i64::from(self.limit - score);
        let mut lo = -score.min(self.n);
        if !ANCHORED {
            lo = i64::from(lo).max(end - remaining) as i32;
        }
        // Each remaining edit changes the diagonal by at most one. An
        // anchored endpoint can have any diagonal at or below `end`.
        let hi = i64::from(score.min(self.m)).min(end + remaining) as i32;
        (lo, hi)
    }
}

fn grow_history(spill: &mut Vec<i32>, prefix: &[i32], end: usize) -> bool {
    if end > SPILL_REACHES {
        return false;
    }
    if end > spill.capacity() {
        let next = (spill.capacity() * 2).max(INLINE_REACHES * 2).max(end).min(SPILL_REACHES);
        if spill.try_reserve_exact(next - spill.len()).is_err() || spill.capacity() > SPILL_REACHES {
            return false;
        }
    }
    if spill.is_empty() {
        spill.extend_from_slice(prefix);
    }
    spill.resize(end, UNREACHABLE);
    true
}

fn global_history_fits(limit: i32, gap: i32) -> bool {
    if limit <= INLINE_SCORE {
        return true;
    }
    let limit = limit as u64;
    let slack = limit - gap as u64;
    // The endpoint cone is bounded by min(2*s, 2*(limit-s), slack) + 1.
    // Start with its widest row everywhere and subtract both end triangles.
    let taper = slack.div_ceil(2);
    let slots = (limit + 1) * (slack + 1 + GUARDS as u64) - 2 * taper * (slack + 1 - taper);
    slots <= INLINE_REACHES as u64
}

fn traceback<const ANCHORED: bool>(
    mut history: &[i32],
    shape: Shape,
    score: i32,
    mut h: i32,
    mut v: i32,
    runs: &mut Vec<Run>,
) -> Result<(), Error> {
    // Bounds depend only on the score, so no row headers are needed.
    for s in (0..score).rev() {
        if h == 0 || v == 0 {
            break;
        }
        let (lo, hi) = shape.band::<ANCHORED>(s);
        let (earlier, row) = history.split_at(history.len() - ((hi - lo) as usize + 1 + GUARDS));
        history = earlier;
        let index = (h - v - lo + 2) as usize;
        let [left, center, right] = row[index - 1..index + 2] else { unreachable!() };
        trace_edit(runs, &mut v, &mut h, left, center, right)?;
    }
    debug_assert!(h == 0 || v == 0);
    push_run(runs, Op::Insert, v as usize)?;
    push_run(runs, Op::Delete, h as usize)?;
    Ok(())
}

#[inline]
pub(super) fn probe<const ANCHORED: bool, const REVERSE: bool>(
    q: &[u8],
    t: &[u8],
    limit: usize,
    frame: Frame,
) -> Result<Probe, Error> {
    if !coordinates_fit(q.len(), t.len()) {
        return Ok(Probe::Switch);
    }
    let shape = Shape { n: q.len() as i32, m: t.len() as i32, limit: limit as i32 };
    wavefront::<ANCHORED, REVERSE>(q, t, shape, frame)
}

#[inline(never)]
fn wavefront<const ANCHORED: bool, const REVERSE: bool>(
    q: &[u8], t: &[u8], mut shape: Shape, frame: Frame,
) -> Result<Probe, Error> {
    let end = shape.m - shape.n;
    let original_limit = shape.limit;
    let gap = if ANCHORED { (-end).max(0) } else { end.abs() };
    if gap > original_limit {
        return Ok(Probe::Miss);
    }
    if !ANCHORED && !global_history_fits(original_limit, gap) {
        // Narrow a loose global cone to a speculative leaf bound. A hit is
        // exact; exhausting a reduced bound only requests a new backend.
        shape.limit = original_limit.min(GLOBAL_LEAF_LIMIT);
        if gap > shape.limit {
            return Ok(Probe::Switch);
        }
    }
    let mut inline = [UNREACHABLE; INLINE_REACHES];
    let mut spill = Vec::new();
    inline[2] = 0; // The first core bytes differ in traversal order.
    let (mut prev_start, mut prev_lo, mut start) = (0, 0, 1 + GUARDS);
    let mut checkpoint = 2 * INLINE_SCORE;

    for score in 1..=shape.limit {
        let (lo, hi) = shape.band::<ANCHORED>(score);
        let width = (hi - lo) as usize + 1;
        let stop = start + width + GUARDS;
        let history = if !ANCHORED || stop <= INLINE_REACHES {
            &mut inline[..stop]
        } else {
            if !grow_history(&mut spill, &inline[..start.min(INLINE_REACHES)], stop) {
                return Ok(Probe::Switch);
            }
            &mut spill[..]
        };
        let (prev, curr) = history.split_at_mut(start);
        // A row boundary moves by at most one diagonal per edit. Two
        // sentinel cells on each side keep every three-way seed in bounds.
        let first = prev_start + (lo - prev_lo + 1) as usize;
        let source = &prev[first..first + width + 2];
        let row = &mut curr[2..2 + width];
        let reached = advance::<ANCHORED, REVERSE>(q, t, lo, source, row, None);
        let h = if ANCHORED {
            reached
        } else {
            ((lo..=hi).contains(&end) && row[(end - lo) as usize] == shape.m).then_some(shape.m)
        };
        if let Some(h) = h {
            let runs = frame.trace_runs(score as usize, shape.n.min(h) as usize)?;
            return frame.finish::<REVERSE>(score as usize, h as usize, !REVERSE, runs, |runs| {
                traceback::<ANCHORED>(prev, shape, score, h, shape.n, runs)
            }).map(Probe::Hit);
        }
        if ANCHORED && score == checkpoint && score < shape.limit {
            let furthest = row.iter().copied().max().unwrap().max(0) as u64;
            if width >= INLINE_SCORE as usize && furthest < score as u64 * MIN_BYTES_PER_EDIT {
                return Ok(Probe::Switch);
            }
            checkpoint = checkpoint.saturating_mul(2);
        }
        (prev_start, prev_lo, start) = (start, lo, stop);
    }
    Ok(if shape.limit == original_limit { Probe::Miss } else { Probe::Switch })
}
