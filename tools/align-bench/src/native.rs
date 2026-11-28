use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, Heuristics, MemoryModel, WFAligner};
use std::ops::Range;

use bitfront::Op;

use crate::{Mode, Task};

impl Mode {
    pub(crate) fn rust_edlib(self, query: &[u8], target: &[u8], cutoff: Option<usize>) -> Option<usize> {
        match (self, cutoff) {
            (Self::Global, None) => rust_edlib::distance(query, target).map(Some),
            (Self::Global, Some(k)) => rust_edlib::distance_at_most(query, target, k),
            (Self::Prefix, None) => rust_edlib::prefix_distance(query, target).map(Some),
            (Self::Prefix, Some(k)) => rust_edlib::prefix_distance_at_most(query, target, k),
            (Self::Suffix, None) => rust_edlib::suffix_distance(query, target).map(Some),
            (Self::Suffix, Some(k)) => rust_edlib::suffix_distance_at_most(query, target, k),
        }.expect("rust-edlib failed to score benchmark input")
    }
}

pub(crate) struct Wfa2 {
    aligner: WFAligner,
    mode: Mode,
    cutoff: Option<usize>,
}

pub(crate) struct WfaPath {
    pub(crate) score: usize,
    pub(crate) target: Range<usize>,
    cigar: Vec<u32>,
    first_run: usize,
}

impl WfaPath {
    pub(crate) fn runs(&self) -> impl Iterator<Item = (Op, usize)> + '_ {
        self.cigar[self.first_run..].iter().copied().map(decode_run)
    }
}

fn decode_run(word: u32) -> (Op, usize) {
    let op = match word & 0xf {
        1 => Op::Insert,
        2 => Op::Delete,
        7 => Op::Match,
        8 => Op::Mismatch,
        _ => panic!("unexpected WFA2 packed CIGAR operation: {word:#x}"),
    };
    (op, (word >> 4) as usize)
}

// WFA2 checks its step limit before extending the last wavefront.
fn max_steps(k: usize) -> i32 {
    k.checked_add(1).and_then(|steps| i32::try_from(steps).ok()).expect("WFA2 bound is too large")
}

impl Wfa2 {
    pub(crate) fn new(mode: Mode, cutoff: Option<usize>, task: Task) -> Self {
        let scope = match task {
            Task::Score => AlignmentScope::Score,
            Task::Path => AlignmentScope::Alignment,
        };
        let builder = WFAligner::builder(scope, MemoryModel::MemoryHigh)
            .edit()
            .with_heuristics(Heuristics::none())
            .with_max_num_threads(1);
        let builder = if let Some(k) = cutoff {
            builder.with_max_alignment_steps(max_steps(k))
        } else {
            builder
        };
        let mut aligner = builder.build().expect("could not construct WFA2 aligner");
        if cutoff.is_some() {
            // The Rust wrapper snapshots native memory_used even on a step
            // limit. Only normal cleanup initializes that field, so complete
            // a zero-cost alignment here before any timed call can miss.
            let initial = aligner.align_end_to_end(b"A", b"A");
            assert_eq!(initial.status, AlignmentStatus::StatusAlgCompleted);
        }
        Self { aligner, mode, cutoff }
    }

    pub(crate) fn score(&mut self, query: &[u8], target: &[u8], cutoff: Option<usize>) -> Option<usize> {
        // One warm aligner per row. Only the mixed boundary rows change k between pairs.
        if self.cutoff != cutoff {
            let (Some(_), Some(k)) = (self.cutoff, cutoff) else {
                panic!("a benchmark row cannot mix bounded and unbounded calls");
            };
            self.aligner.set_max_alignment_steps(max_steps(k));
            self.cutoff = cutoff;
        }
        let _ = i32::try_from(query.len()).expect("WFA2 query is too long");
        let target_len = i32::try_from(target.len()).expect("WFA2 target is too long");
        // A fresh native call resets the frontiers even after a max-step
        // pause, so keep using the same aligner and its allocated buffers.
        let result = match self.mode {
            Mode::Global => self.aligner.align_end_to_end(query, target),
            Mode::Prefix => self.aligner.align_ends_free(query, 0, 0, target, 0, target_len),
            Mode::Suffix => self.aligner.align_ends_free(query, 0, 0, target, target_len, 0),
        };
        match result.status {
            // The completed unit-edit wavefront score is the nonnegative edit distance.
            AlignmentStatus::StatusAlgCompleted => {
                Some(usize::try_from(result.score).expect("WFA2 returned an invalid distance"))
            }
            AlignmentStatus::StatusMaxStepsReached if self.cutoff.is_some() => None,
            _ => panic!("WFA2 {:?} failed: {result:?}", self.mode),
        }
    }

    pub(crate) fn path(&mut self, query: &[u8], target: &[u8], cutoff: Option<usize>) -> Option<WfaPath> {
        let score = self.score(query, target, cutoff)?;
        // Alignment scope has already run native traceback. Pack its CIGAR
        // (including mismatches) and own the resulting run vector here, in the
        // timed call. SAM orientation is already correct for query=pattern
        // and target=text; replay decodes each run without expanding it.
        let mut cigar = self.aligner.sam_packed_cigar(true);
        let mut target = 0..target.len();
        let mut first_run = 0;
        // WFA2 pads the traceback to both input lengths. With the entire
        // target edge free, a target-only run there is necessarily free:
        // charging any part of that run would contradict optimality.
        match self.mode {
            Mode::Global => {}
            Mode::Prefix => {
                if let Some((Op::Delete, len)) = cigar.last().copied().map(decode_run) {
                    let _ = cigar.pop();
                    target.end -= len;
                }
            }
            Mode::Suffix => {
                if let Some((Op::Delete, len)) = cigar.first().copied().map(decode_run) {
                    first_run = 1;
                    target.start = len;
                }
            }
        }
        Some(WfaPath { score, target, cigar, first_run })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_scores_are_inclusive_after_step_limit_misses() {
        let mut wfa2 = Wfa2::new(Mode::Global, Some(1), Task::Score);
        let pairs: &[(&[u8], &[u8], usize, Option<usize>)] = &[
            (b"AA", b"CC", 1, None),
            (b"AA", b"CC", 2, Some(2)),
            (b"AA", b"CC", 0, None),
            (b"AC", b"AG", 1, Some(1)),
            (b"AC", b"ACG", 1, Some(1)),
            (b"AC", b"", 1, None),
            (b"TT", b"TT", 1, Some(0)),
        ];
        for &(query, target, k, expected) in pairs {
            assert_eq!(Mode::Global.rust_edlib(query, target, Some(k)), expected);
            assert_eq!(wfa2.score(query, target, Some(k)), expected);
        }
    }

    #[test]
    fn anchored_scores_consume_the_query_and_free_only_the_other_end() {
        let query = b"CAG";
        for (mode, target, distance) in [
            (Mode::Prefix, b"CAGTTT", 0),
            (Mode::Suffix, b"CAGTTT", 3),
            (Mode::Suffix, b"TTTCAG", 0),
        ] {
            assert_eq!(mode.rust_edlib(query, target, None), Some(distance));
            assert_eq!(Wfa2::new(mode, None, Task::Score).score(query, target, None), Some(distance));
        }
    }

    #[test]
    fn bounded_anchors_can_hit_at_the_limit_after_a_wrong_end_miss() {
        let (prefix, suffix) = (b"CATTTT", b"TTTCAT");
        for (mode, hit, miss) in [(Mode::Prefix, prefix, suffix), (Mode::Suffix, suffix, prefix)] {
            let mut wfa2 = Wfa2::new(mode, Some(1), Task::Score);
            for (target, expected) in [(miss, None), (hit, Some(1))] {
                assert_eq!(mode.rust_edlib(b"CAG", target, Some(1)), expected);
                assert_eq!(wfa2.score(b"CAG", target, Some(1)), expected);
            }
        }
    }

    #[test]
    fn bounded_paths_return_owned_tracebacks_at_the_limit_after_misses() {
        let (prefix, suffix) = (b"ACGTCCC", b"CCCACGT");
        for (mode, hit, miss, span) in [
            (Mode::Prefix, prefix, suffix, 0..4),
            (Mode::Suffix, suffix, prefix, 3..7),
        ] {
            let mut wfa2 = Wfa2::new(mode, Some(1), Task::Path);
            assert!(wfa2.path(b"ACAT", miss, Some(1)).is_none());
            let path = wfa2.path(b"ACAT", hit, Some(1)).unwrap();
            let exact = wfa2.path(b"ACGT", hit, Some(1)).unwrap();
            assert_eq!((exact.score, exact.target.clone()), (0, span.clone()));
            assert_eq!(exact.runs().collect::<Vec<_>>(), [(Op::Match, 4)]);
            assert_eq!((path.score, path.target.clone()), (1, span));
            assert_eq!(
                path.runs().collect::<Vec<_>>(),
                [(Op::Match, 2), (Op::Mismatch, 1), (Op::Match, 1)]
            );
        }
    }

    #[test]
    fn anchored_paths_trim_only_free_target_gaps() {
        let cases: [(Mode, &[u8], &[u8], Range<usize>, [(Op, usize); 2]); 4] = [
            (Mode::Prefix, b"ACGT", b"TACGTCCC", 0..5, [(Op::Delete, 1), (Op::Match, 4)]),
            (Mode::Suffix, b"ACGT", b"CCCACGTA", 3..8, [(Op::Match, 4), (Op::Delete, 1)]),
            (Mode::Prefix, b"TACGT", b"ACGTCCC", 0..4, [(Op::Insert, 1), (Op::Match, 4)]),
            (Mode::Suffix, b"ACGTA", b"CCCACGT", 3..7, [(Op::Match, 4), (Op::Insert, 1)]),
        ];
        for (mode, query, target, span, runs) in cases {
            let path = Wfa2::new(mode, None, Task::Path).path(query, target, None).unwrap();
            assert_eq!(path.score, 1);
            assert_eq!(path.target, span);
            assert_eq!(path.runs().collect::<Vec<_>>(), runs);
        }
    }

    #[test]
    fn paths_keep_anchor_coordinates_with_empty_inputs() {
        for (mode, span) in [(Mode::Prefix, 0..0), (Mode::Suffix, 3..3)] {
            let mut wfa2 = Wfa2::new(mode, None, Task::Path);
            let free = wfa2.path(b"", b"ACG", None).unwrap();
            assert_eq!(free.score, 0);
            assert_eq!(free.target, span);
            assert_eq!(free.runs().next(), None);
            let gap = wfa2.path(b"ACG", b"", None).unwrap();
            assert_eq!(gap.score, 3);
            assert_eq!(gap.target, 0..0);
            assert_eq!(gap.runs().collect::<Vec<_>>(), [(Op::Insert, 3)]);
            let empty = wfa2.path(b"", b"", None).unwrap();
            assert_eq!(empty.score, 0);
            assert_eq!(empty.target, 0..0);
            assert_eq!(empty.runs().next(), None);
        }
    }
}
