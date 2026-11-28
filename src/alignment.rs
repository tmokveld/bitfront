use crate::{Error, flank, kernel};
use rust_edlib::{AlignmentMode, GlobalProbe};
use std::ops::Range;

mod checkpoint;
mod trace;

/// Unit operations in read/reference order: `Insert` consumes only a query byte,
/// `Delete` only a target byte, and `Match` and `Mismatch` consume one of each.
pub use rust_edlib::EditOperation as Op;

/// A positive number of consecutive operations of the same kind.
pub use rust_edlib::Run;

/// One optimal alignment of the entire query. Tied paths and target spans are
/// not canonical. An empty target span can still have query-only operations.
pub use rust_edlib::RunAlignment as Alignment;

/// Aligns the entire query to the entire target. For a score alone, use [`crate::distance`].
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn align(query: &[u8], target: &[u8]) -> Result<Alignment, Error> {
    unbounded::<false, false>(query, target)
}

/// Returns a global alignment if its exact score is at most `k` (inclusive).
/// `Ok(None)` proves that the distance exceeds `k`.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn align_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Alignment>, Error> {
    search::<false, false>(query, target, Some(k))
}

/// Aligns the entire query to an optimal target prefix, possibly empty. Trailing
/// target bytes are free. For a score alone, use [`crate::prefix_distance`].
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
///
/// ```
/// use bitfront::{Op, Run, prefix_align};
/// let hit = prefix_align(b"ACGT", b"ATCGTzz")?;
/// assert_eq!((hit.score, hit.target), (1, 0..5));
/// assert_eq!(hit.runs, [Run { op: Op::Match, len: 1 },
///                      Run { op: Op::Delete, len: 1 }, // target-only T
///                      Run { op: Op::Match, len: 3 }]);
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn prefix_align(query: &[u8], target: &[u8]) -> Result<Alignment, Error> {
    unbounded::<true, false>(query, target)
}

/// Returns a [`prefix_align`] hit at or below `k` (inclusive), or `None` if every prefix costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn prefix_align_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Alignment>, Error> {
    search::<true, false>(query, target, Some(k))
}

/// Aligns the entire query to an optimal target suffix, possibly empty. Both
/// range and runs are in forward coordinates; the leading target flank is free.
/// For a score alone, use [`crate::suffix_distance`].
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn suffix_align(query: &[u8], target: &[u8]) -> Result<Alignment, Error> {
    unbounded::<true, true>(query, target)
}

/// Returns a [`suffix_align`] hit at or below `k` (inclusive), or `None` if every suffix costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn suffix_align_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Alignment>, Error> {
    search::<true, true>(query, target, Some(k))
}

/// Aligns the entire query to an optimal target substring, possibly empty.
/// Both target flanks are free. The half-open range uses original target
/// coordinates; runs are coalesced and in forward order. Tied spans are not
/// canonical. For the score and range without runs, use [`crate::infix_locate`].
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn infix_align(query: &[u8], target: &[u8]) -> Result<Alignment, Error> {
    Ok(rust_edlib::align_runs(query, target, AlignmentMode::Infix, None, &[])?
        .expect("an unbounded alignment has a path"))
}

/// Returns an optimal [`infix_align`] path if its score is at most `k`
/// (inclusive), or `None` if every target substring costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
///
/// ```
/// let hit = bitfront::infix_align_at_most(b"ab", b"_ax_ab_", 1)?.unwrap();
/// assert_eq!((hit.score, hit.target), (0, 4..6));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn infix_align_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Alignment>, Error> {
    rust_edlib::align_runs(query, target, AlignmentMode::Infix, Some(k), &[])
}

fn unbounded<const ANCHORED: bool, const REVERSE: bool>(query: &[u8], target: &[u8]) -> Result<Alignment, Error> {
    Ok(search::<ANCHORED, REVERSE>(query, target, None)?.expect("an unbounded alignment has a path"))
}

const HISTORY_BYTES: usize = 256 * 1024;

fn push_run(runs: &mut Vec<Run>, op: Op, len: usize) -> Result<(), Error> {
    if len != 0 {
        if let Some(last) = runs.last_mut().filter(|r| r.op == op) {
            last.len += len;
        } else {
            runs.try_reserve(1).map_err(|_| Error::InputTooLong)?;
            runs.push(Run { op, len });
        }
    }
    Ok(())
}

/// Undo one edit at an exact-score cell with both coordinates positive.
/// Some optimal last edit ends no later than the largest seed, and matches
/// from there to j. If the seed overshoots j, clip its predecessor along
/// the diagonal: shorter prefixes cost no more, and a cost below s-1 would
/// contradict optimality. Thus scoring's clamped reaches suffice for trace.
#[inline]
fn trace_edit(
    runs: &mut Vec<Run>, i: &mut i32, j: &mut i32, left: i32, center: i32, right: i32,
) -> Result<(), Error> {
    let (sub, del) = (center + 1, left + 1);
    let best = sub.max(del).max(right);
    let matched = (*j - best).max(0);
    debug_assert!(matched <= (*i).min(*j));
    let op = if sub == best { Op::Mismatch } else if right == best { Op::Insert } else { Op::Delete };
    push_run(runs, Op::Match, matched as usize)?;
    push_run(runs, op, 1)?;
    *i -= matched + i32::from(op != Op::Delete);
    *j -= matched + i32::from(op != Op::Insert);
    Ok(())
}

// Both backends trace only the trimmed core. Keep its framing in original
// coordinates so a suffix's free prefix can never leak into a run or range.
#[derive(Clone, Copy)]
struct Frame {
    prefix: usize,
    suffix: usize,
    target_len: usize,
}

impl Frame {
    fn span<const REVERSE: bool>(self, used: usize) -> Range<usize> {
        let length = self.prefix + used + self.suffix;
        if REVERSE { self.target_len - length..self.target_len } else { 0..length }
    }

    fn trace_runs(self, score: usize, shorter_core: usize) -> Result<Vec<Run>, Error> {
        const MAX_INITIAL_RUNS: usize = 32;
        // d edits separate at most d+1 match runs. Likewise, m bytes of the
        // shorter aligned span separate at most m+1 runs of longer-side gaps.
        let separators = score.min(self.prefix + shorter_core + self.suffix).min(MAX_INITIAL_RUNS / 2);
        let mut runs = Vec::new();
        runs.try_reserve_exact((2 * separators + 1).min(MAX_INITIAL_RUNS)).map_err(|_| Error::InputTooLong)?;
        Ok(runs)
    }

    fn finish<const REVERSE: bool>(
        self, score: usize, used: usize, backwards: bool, mut runs: Vec<Run>,
        emit: impl FnOnce(&mut Vec<Run>) -> Result<(), Error>,
    ) -> Result<Alignment, Error> {
        let (first, last) = if backwards { (self.suffix, self.prefix) } else { (self.prefix, self.suffix) };
        push_run(&mut runs, Op::Match, first)?;
        emit(&mut runs)?;
        push_run(&mut runs, Op::Match, last)?;
        if backwards {
            runs.reverse();
        }
        Ok(Alignment { score, target: self.span::<REVERSE>(used), runs })
    }

    fn single<const REVERSE: bool>(self, op: Op, len: usize) -> Result<Alignment, Error> {
        let score = if op == Op::Match { 0 } else { len };
        let used = if op == Op::Insert { 0 } else { len };
        self.finish::<REVERSE>(score, used, false, Vec::new(), |runs| push_run(runs, op, len))
    }

    fn restore<const REVERSE: bool>(self, hit: Option<Alignment>) -> Result<Option<Alignment>, Error> {
        let Some(mut hit) = hit else { return Ok(None) };
        // The backend span is in the trimmed window, not the original target.
        hit.target = self.span::<REVERSE>(hit.target.len());
        // Each trimmed end borders a mismatch in the core, so neither added
        // match run can adjoin another match. Keep the backend's vector.
        let runs = &mut hit.runs;
        let extra = usize::from(self.prefix != 0) + usize::from(self.suffix != 0);
        runs.try_reserve_exact(extra).map_err(|_| Error::InputTooLong)?;
        if self.prefix != 0 {
            runs.insert(0, Run { op: Op::Match, len: self.prefix });
        }
        if self.suffix != 0 {
            runs.push(Run { op: Op::Match, len: self.suffix });
        }
        Ok(Some(hit))
    }
}

fn search<const ANCHORED: bool, const REVERSE: bool>(
    query: &[u8], target: &[u8], cutoff: Option<usize>,
) -> Result<Option<Alignment>, Error> {
    let k = cutoff.unwrap_or(usize::MAX);
    let gap = if ANCHORED { query.len().saturating_sub(target.len()) } else { query.len().abs_diff(target.len()) };
    if gap > k {
        return Ok(None);
    }
    let mut frame = Frame { prefix: 0, suffix: 0, target_len: target.len() };
    if k == 0 {
        let exact = if !ANCHORED { target == query } else if REVERSE { target.ends_with(query) } else { target.starts_with(query) };
        return exact.then(|| frame.single::<REVERSE>(Op::Match, query.len())).transpose();
    }
    if ANCHORED && k == gap {
        // At the positive length gap only the whole target is long enough.
        return rust_edlib::align_runs(query, target, AlignmentMode::Global, Some(k), &[]);
    }
    let (query, target) = if ANCHORED {
        let (q, t) = flank::trimmed_window::<REVERSE>(query, target, k);
        let pinned = query.len() - q.len();
        if REVERSE { frame.suffix = pinned } else { frame.prefix = pinned }
        (q, t)
    } else {
        (frame.prefix, frame.suffix) = kernel::matching_ends(query, target);
        (&query[frame.prefix..query.len() - frame.suffix], &target[frame.prefix..target.len() - frame.suffix])
    };
    // A depleted flank query has an empty window. Any other empty side costs
    // precisely the mandatory gap already checked above.
    let single = match (query.len(), target.len()) {
        (0, n) => Some((Op::Delete, n)),
        (n, 0) => Some((Op::Insert, n)),
        (1, 1) => Some((Op::Mismatch, 1)), // The pinned match was trimmed.
        _ => None,
    };
    if let Some((op, len)) = single {
        return frame.single::<REVERSE>(op, len).map(Some);
    }
    let upper = if ANCHORED { query.len() } else { query.len().max(target.len()) };
    let limit = k.min(upper);
    let mut native = trace::probe::<ANCHORED, REVERSE>(query, target, limit, frame)?;
    let mut global = None;
    if !ANCHORED && matches!(native, trace::Probe::Switch) {
        // Keep the cheap inline leaf above, especially for tight cutoffs.
        // Wider skews go directly to the backend rather than building both
        // checkpoint fronts just to cover the mandatory gap.
        let skewed = gap > 32 && gap >= query.len().min(target.len()).div_ceil(2);
        if !skewed && kernel::coordinates_fit(query.len(), target.len()) {
            let wide = kernel::wide_cutoff(query.len(), target.len(), limit);
            global = wide.then(|| GlobalProbe::new(query, target));
            if let Some(probe) = global.as_ref().filter(|probe| probe.diagonal_score().is_some()) {
                return frame.restore::<REVERSE>(probe.align_runs(Some(limit))?);
            }
            native = checkpoint::probe(query, target, limit, wide, frame)?;
        }
    }
    match native {
        trace::Probe::Hit(hit) => Ok(Some(hit)),
        trace::Probe::Miss => Ok(None),
        // Neither a partial wavefront nor a failed internal budget certifies
        // a miss. Restart the trimmed core with the caller's original cutoff.
        trace::Probe::Switch => {
            let hit = if let Some(probe) = global {
                probe.align_runs(cutoff)?
            } else {
                let mode = if REVERSE { AlignmentMode::Suffix } else if ANCHORED { AlignmentMode::Prefix } else { AlignmentMode::Global };
                rust_edlib::align_runs(query, target, mode, cutoff, &[])?
            };
            frame.restore::<REVERSE>(hit)
        }
    }
}
