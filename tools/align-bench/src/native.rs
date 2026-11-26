use rust_wfa2::aligner::{AlignmentScope, AlignmentStatus, Heuristics, MemoryModel, WFAligner};

use crate::Mode;

impl Mode {
    pub(crate) fn rust_edlib(self, query: &[u8], target: &[u8]) -> Option<usize> {
        match self {
            Self::Global => rust_edlib::distance(query, target).map(Some),
            Self::Prefix => rust_edlib::prefix_distance(query, target).map(Some),
            Self::Suffix => rust_edlib::suffix_distance(query, target).map(Some),
            Self::AtMost(k) => rust_edlib::align(query, target, &rust_edlib::AlignConfig {
                max_distance: Some(k), ..Default::default()
            }).map(|result| result.edit_distance()),
        }.expect("rust-edlib failed to score benchmark input")
    }
}

pub(crate) struct Wfa2 {
    aligner: WFAligner,
    mode: Mode,
}

impl Wfa2 {
    pub(crate) fn new(mode: Mode) -> Self {
        let builder = WFAligner::builder(AlignmentScope::Score, MemoryModel::MemoryHigh)
            .edit()
            .with_heuristics(Heuristics::none())
            .with_max_num_threads(1);
        let builder = if let Mode::AtMost(k) = mode {
            // The limit is checked after computing a wavefront, but before
            // extending it. Compute k+1 to allow the entire score-k step.
            let max_steps = k
                .checked_add(1)
                .and_then(|steps| i32::try_from(steps).ok())
                .expect("WFA2 bound is too large");
            builder.with_max_alignment_steps(max_steps)
        } else {
            builder
        };
        let mut aligner = builder.build().expect("could not construct WFA2 aligner");
        if matches!(mode, Mode::AtMost(_)) {
            // The Rust wrapper snapshots native memory_used even on a step
            // limit. Only normal cleanup initializes that field, so complete
            // a zero-cost alignment here before any timed call can miss.
            let initial = aligner.align_end_to_end(b"A", b"A");
            assert_eq!(initial.status, AlignmentStatus::StatusAlgCompleted);
        }
        Self { aligner, mode }
    }

    pub(crate) fn score(&mut self, query: &[u8], target: &[u8]) -> Option<usize> {
        let _ = i32::try_from(query.len()).expect("WFA2 query is too long");
        let target_len = i32::try_from(target.len()).expect("WFA2 target is too long");
        // A fresh native call resets the frontiers even after a max-step
        // pause, so keep using the same aligner and its allocated buffers.
        let result = match self.mode {
            Mode::Global | Mode::AtMost(_) => self.aligner.align_end_to_end(query, target),
            Mode::Prefix => self.aligner.align_ends_free(query, 0, 0, target, 0, target_len),
            Mode::Suffix => self.aligner.align_ends_free(query, 0, 0, target, target_len, 0),
        };
        match result.status {
            // The completed unit-edit wavefront score is the nonnegative
            // edit distance. There is no need to request a CIGAR.
            AlignmentStatus::StatusAlgCompleted => {
                Some(usize::try_from(result.score).expect("WFA2 returned an invalid distance"))
            }
            AlignmentStatus::StatusMaxStepsReached if matches!(self.mode, Mode::AtMost(_)) => None,
            _ => panic!("WFA2 {:?} failed: {result:?}", self.mode),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounded_scores_are_inclusive_after_step_limit_misses() {
        let mut wfa2 = Wfa2::new(Mode::AtMost(1));
        let pairs: &[(&[u8], &[u8], Option<usize>)] = &[
            (b"AA", b"CC", None),
            (b"AC", b"AG", Some(1)),
            (b"AC", b"ACG", Some(1)),
            (b"AC", b"", None),
            (b"TT", b"TT", Some(0)),
        ];
        for &(query, target, expected) in pairs {
            assert_eq!(Mode::AtMost(1).rust_edlib(query, target), expected);
            assert_eq!(wfa2.score(query, target), expected);
        }
    }

    #[test]
    fn anchored_scores_consume_the_query_and_free_only_the_other_end() {
        let query = b"CAG";
        let target = b"CAGTTT";
        assert_eq!(Mode::Prefix.rust_edlib(query, target), Some(0));
        assert_eq!(Wfa2::new(Mode::Prefix).score(query, target), Some(0));
        assert_eq!(Mode::Suffix.rust_edlib(query, target), Some(3));
        assert_eq!(Wfa2::new(Mode::Suffix).score(query, target), Some(3));
        assert_eq!(Mode::Suffix.rust_edlib(query, b"TTTCAG"), Some(0));
        assert_eq!(Wfa2::new(Mode::Suffix).score(query, b"TTTCAG"), Some(0));
    }
}
