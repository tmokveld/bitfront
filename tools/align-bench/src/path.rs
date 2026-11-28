use crate::{Mode, native::Wfa2, timing::Pair};
use bitfront::Op;
use rust_edlib::{AlignmentMode, RunAlignment};
use std::ops::Range;

impl Mode {
    pub(super) fn bitfront_path(self, q: &[u8], t: &[u8], k: Option<usize>) -> Option<bitfront::Alignment> {
        match (self, k) {
            (Self::Global, None) => bitfront::align(q, t).map(Some),
            (Self::Prefix, None) => bitfront::prefix_align(q, t).map(Some),
            (Self::Suffix, None) => bitfront::suffix_align(q, t).map(Some),
            (Self::Global, Some(k)) => bitfront::align_at_most(q, t, k),
            (Self::Prefix, Some(k)) => bitfront::prefix_align_at_most(q, t, k),
            (Self::Suffix, Some(k)) => bitfront::suffix_align_at_most(q, t, k),
        }.expect("Bitfront traceback failed")
    }

    pub(super) fn edlib_path(self, q: &[u8], t: &[u8], k: Option<usize>) -> Option<RunAlignment> {
        let mode = match self {
            Self::Global => AlignmentMode::Global,
            Self::Prefix => AlignmentMode::Prefix,
            Self::Suffix => AlignmentMode::Suffix,
        };
        rust_edlib::align_runs(q, t, mode, k, &[]).expect("rust-edlib traceback failed")
    }
}

// Check each native representation without expanding runs or requiring tied
// paths (or endpoints) to be identical. This is never in the timed closure.
fn replay(mode: Mode, pair: &Pair, span: Range<usize>, ops: impl Iterator<Item = (Op, usize)>) {
    if mode != Mode::Suffix { assert_eq!(span.start, 0); }
    if mode != Mode::Prefix { assert_eq!(span.end, pair.target.len()); }
    let q = &pair.query;
    let t = &pair.target[span];
    let (mut i, mut j, mut edits) = (0, 0, 0);
    for (op, n) in ops {
        assert!(n > 0);
        match op {
            Op::Match => { assert_eq!(&q[i..i+n], &t[j..j+n]); i += n; j += n; }
            Op::Mismatch => {
                assert!(q[i..i+n].iter().zip(&t[j..j+n]).all(|(a, b)| a != b));
                i += n; j += n; edits += n;
            }
            Op::Insert => { i += n; edits += n; }
            Op::Delete => { j += n; edits += n; }
        }
    }
    assert_eq!((i, j, edits), (q.len(), t.len(), pair.distance), "{mode:?} replay");
}

pub(super) fn validate(mode: Mode, pair: &Pair, expected: Option<usize>, wfa2: &mut Wfa2) {
    let k = pair.cutoff;
    let (q, t) = (&pair.query, &pair.target);
    let check_runs = |engine: &str, hit: Option<RunAlignment>| {
        assert_eq!(hit.as_ref().map(|a| a.score), expected, "{engine} {mode:?}, k={k:?}");
        if let Some(hit) = hit {
            assert!(hit.runs.windows(2).all(|r| r[0].op != r[1].op), "{engine}: unmerged runs");
            replay(mode, pair, hit.target, hit.runs.into_iter().map(|r| (r.op, r.len)));
        }
    };
    check_runs("Bitfront", mode.bitfront_path(q, t, k));
    check_runs("rust-edlib", mode.edlib_path(q, t, k));
    let w = wfa2.path(q, t, k);
    assert_eq!(w.as_ref().map(|a| a.score), expected, "WFA2 {mode:?}, k={k:?}");
    if let Some(hit) = w {
        replay(mode, pair, hit.target.clone(), hit.runs());
    }
}
