//! Global traceback beyond the short trace. Keep the cheap, narrow rows;
//! thereafter save every sixteenth row from each end. Once the fronts meet,
//! only a cone of at most 2 * STRIDE + 1 diagonals needs replay per checkpoint.

use super::trace::Probe;
use super::{Alignment, Error, Frame, HISTORY_BYTES, Op, Run, push_run, trace_edit};
use crate::kernel::{Front, MIN_BYTES_PER_EDIT, UNREACHABLE, advance, meeting, uneven_progress};

const STRIDE: usize = 16;
const EARLY: usize = 2 * STRIDE;
type Small = [i32; EARLY * (EARLY + 4)];
const WIDTH: usize = 2 * STRIDE + 3; // One guard at either end.
type Tile = [[i32; WIDTH]; STRIDE];
const PER_END: usize = (HISTORY_BYTES - size_of::<Tile>() - 2 * size_of::<Small>()) / (2 * size_of::<i32>());

// Score d has at most 2*d+1 offsets and four guard slots. The first d
// guarded rows occupy d*(d+4) slots, even when a rectangle clips their ends.
fn small_start(d: usize) -> usize {
    d * (d + 4) + 2
}

struct History {
    small: Small,
    saved: Vec<i32>, // Concatenated checkpoints; replay consumes these from the tail.
}

struct Trail {
    front: Front,
    history: History,
    depth: i32,
    upper: i32,
}

impl Trail {
    fn new() -> Self {
        let history = History { small: [UNREACHABLE; EARLY * (EARLY + 4)], saved: Vec::new() };
        Self { front: Front::new(), history, depth: 0, upper: 0 }
    }

    fn step<const REVERSE: bool>(&mut self, q: &[u8], t: &[u8]) -> Result<bool, Error> {
        // Save on departure: the last saved depth is strictly below the
        // eventual split depth, including when we meet exactly on a checkpoint.
        let depth = self.depth as usize;
        if depth < EARLY {
            let row = self.front.offsets();
            let start = small_start(depth);
            self.history.small[start..start + row.len()].copy_from_slice(row);
        } else if depth % STRIDE == 0 {
            let row = self.front.offsets();
            let saved = &mut self.history.saved;
            let len = saved.len();
            let need = len + row.len();
            if need > PER_END {
                return Ok(false);
            }
            if need > saved.capacity() {
                let capacity = need.max(2 * saved.capacity()).min(PER_END);
                if saved.try_reserve_exact(capacity - len).is_err() || saved.capacity() > PER_END {
                    return Ok(false);
                }
            }
            saved.extend_from_slice(row);
        }
        self.depth += 1;
        self.front.expand::<REVERSE>(q, t, self.depth, Some(&mut self.upper))?;
        Ok(true)
    }
}

// (i,j) is reachable at exactly `depth`, but need not be the furthest cell on
// its diagonal (in particular on the reverse side of the meeting point).
fn replay<const REVERSE: bool>(
    q: &[u8], t: &[u8], history: &History, mut depth: i32, mut i: i32, mut j: i32,
    runs: &mut Vec<Run>,
) -> Result<(), Error> {
    let (n, m) = (q.len() as i32, t.len() as i32);
    let mut saved = &history.saved[..];
    while depth > EARLY as i32 {
        let base = (depth - 1) / STRIDE as i32 * STRIDE as i32;
        let radius = depth - base;
        let center = j - i;
        let origin = center - radius - 1; // Local slot zero is a guard.
        let mut tile = [[UNREACHABLE; WIDTH]; STRIDE];
        let (left, right) = (-base.min(n), base.min(m));
        let width = (right - left + 1) as usize;
        let (rest, row) = saved.split_at(saved.len() - width);
        saved = rest;
        let (lo, hi) = (left.max(center - radius), right.min(center + radius));
        let start = (lo - origin) as usize;
        let width = (hi - lo + 1) as usize;
        tile[0][start..start + width].copy_from_slice(&row[(lo - left) as usize..=(hi - left) as usize]);
        for d in 1..radius {
            let s = base + d;
            let r = radius - d;
            let (lo, hi) = ((center - r).max(-s.min(n)), (center + r).min(s.min(m)));
            let start = (lo - origin) as usize;
            let end = (hi - origin + 1) as usize;
            let (earlier, here) = tile.split_at_mut(d as usize);
            let _ = advance::<false, REVERSE>(q, t, lo, &earlier[d as usize - 1][start - 1..end + 1],
                &mut here[0][start..end], None);
        }
        for row in tile[..radius as usize].iter().rev() {
            if i == 0 || j == 0 {
                debug_assert_eq!(depth, i + j);
                return push_run(runs, if i == 0 { Op::Delete } else { Op::Insert }, depth as usize);
            }
            let z = (j - i - origin) as usize;
            trace_edit(runs, &mut i, &mut j, row[z - 1], row[z], row[z + 1])?;
            depth -= 1;
        }
    }
    while depth != 0 {
        if i == 0 || j == 0 {
            debug_assert_eq!(depth, i + j);
            return push_run(runs, if i == 0 { Op::Delete } else { Op::Insert }, depth as usize);
        }
        depth -= 1;
        let z = (small_start(depth as usize) as i32 + (j - i) + depth.min(n)) as usize;
        let row = &history.small;
        trace_edit(runs, &mut i, &mut j, row[z - 1], row[z], row[z + 1])?;
    }
    debug_assert_eq!((i, j), (0, 0));
    Ok(())
}

pub(super) fn probe(q: &[u8], t: &[u8], limit: usize, mut wide: bool, frame: Frame) -> Result<Probe, Error> {
    let (n, m) = (q.len() as i32, t.len() as i32);
    let (mut f, mut r) = (Trail::new(), Trail::new());
    let mut inspect = 64usize;
    for score in 1..=limit {
        let room = if score % 2 == 1 { f.step::<false>(q, t)? } else { r.step::<true>(q, t)? };
        if !room {
            return Ok(Probe::Switch);
        }
        if f.upper >= m - r.upper {
            if let Some((z, j)) = meeting(&f.front, &r.front, m - n, m) {
                let i = j - z;
                drop((f.front, r.front));
                let mut runs = frame.trace_runs(score, q.len().min(t.len()))?;
                replay::<false>(q, t, &f.history, f.depth, i, j, &mut runs)?;
                push_run(&mut runs, Op::Match, frame.prefix)?;
                runs.reverse();
                replay::<true>(q, t, &r.history, r.depth, n - i, m - j, &mut runs)?;
                push_run(&mut runs, Op::Match, frame.suffix)?;
                return Ok(Probe::Hit(Alignment { score, target: 0..frame.target_len, runs }));
            }
        }
        if score == inspect && score < limit {
            let reached = f.front.furthest() as u64 + r.front.furthest() as u64;
            let s = score as u64;
            // Unlike score-only, a near-bound hit benefits from retained history.
            // Switch early only if current progress projects over twice the limit.
            // The full original cutoff goes to Myers; a projection never proves a miss.
            let far_outside = wide && score <= limit / 2 && reached < s * MIN_BYTES_PER_EDIT
                && reached * (limit as u64) < s * m as u64 / 2;
            if reached < 4 * s || (far_outside && !uneven_progress(q, t, s, reached)) {
                return Ok(Probe::Switch);
            }
            if far_outside { wide = false; }
            inspect = inspect.saturating_mul(2);
        }
    }
    Ok(Probe::Miss)
}

#[cfg(test)]
mod tests {
    use super::{Frame, Probe, probe};

    #[test]
    fn exhausted_bound_is_a_miss_even_at_a_density_checkpoint() {
        let (q, t) = ([b'A'; 1024], [b'C'; 1024]);
        let frame = Frame { prefix: 0, suffix: 0, target_len: t.len() };
        assert!(matches!(probe(&q, &t, 64, false, frame), Ok(Probe::Miss)));
    }
}
