use super::{scalar, scalar_row};
use crate::{
    Alignment, Error, Op, align, align_at_most, infix_align, infix_align_at_most, infix_distance,
    infix_distance_at_most, infix_locate, infix_locate_at_most, prefix_align, prefix_align_at_most,
    suffix_align, suffix_align_at_most,
};

#[derive(Clone, Copy, Debug)]
enum Anchor {
    Global,
    Prefix,
    Suffix,
    Infix,
}

type Unbounded = fn(&[u8], &[u8]) -> Result<Alignment, Error>;
type Bounded = fn(&[u8], &[u8], usize) -> Result<Option<Alignment>, Error>;

fn take<'a>(rest: &mut &'a [u8], len: usize) -> &'a [u8] {
    let (head, tail) = rest.split_at_checked(len).expect("run overruns its input");
    *rest = tail;
    head
}

fn replay(anchor: Anchor, mut query: &[u8], target: &[u8], path: &Alignment, optimum: usize) {
    let mut span = target.get(path.target.clone()).expect("invalid target range");
    match anchor {
        Anchor::Global => assert_eq!(path.target, 0..target.len()),
        Anchor::Prefix => assert_eq!(path.target.start, 0),
        Anchor::Suffix => assert_eq!(path.target.end, target.len()),
        Anchor::Infix => {} // Neither endpoint is pinned.
    }

    let mut cost = 0;
    let mut previous = None;
    for run in &path.runs {
        let (op, len) = (run.op, run.len);
        assert!(len > 0, "empty {op:?} run");
        assert_ne!(previous, Some(op), "adjacent {op:?} runs are not coalesced");
        previous = Some(op);
        match op {
            Op::Match => assert_eq!(take(&mut query, len), take(&mut span, len)),
            Op::Mismatch => assert!(
                take(&mut query, len)
                    .iter()
                    .zip(take(&mut span, len))
                    .all(|(q, t)| q != t),
                "mismatch run contains equal bytes"
            ),
            Op::Insert => {
                take(&mut query, len);
            }
            Op::Delete => {
                take(&mut span, len);
            }
        }
        if op != Op::Match {
            cost += len;
        }
    }
    assert!(query.is_empty(), "unconsumed query bytes");
    assert!(span.is_empty(), "unconsumed bytes in the aligned target range");
    assert_eq!(cost, path.score, "script disagrees with the reported score");
    assert_eq!(cost, optimum, "{anchor:?} script is not optimal");
}

fn check_paths(anchor: Anchor, query: &[u8], target: &[u8], score: usize, cutoffs: &[usize]) {
    let (unbounded, bounded): (Unbounded, Bounded) = match anchor {
        Anchor::Global => (align, align_at_most),
        Anchor::Prefix => (prefix_align, prefix_align_at_most),
        Anchor::Suffix => (suffix_align, suffix_align_at_most),
        Anchor::Infix => (infix_align, infix_align_at_most),
    };
    replay(anchor, query, target, &unbounded(query, target).unwrap(), score);
    for &k in cutoffs {
        let path = bounded(query, target, k).unwrap();
        assert_eq!(path.is_some(), score <= k, "{anchor:?}, optimum={score}, k={k}");
        if let Some(path) = path {
            replay(anchor, query, target, &path, score);
        }
    }
}

#[test]
fn exhaustive_binary_alignment_paths_through_four_bytes() {
    let mut words = Vec::new();
    for len in 0..=4 {
        for bits in 0..(1 << len) {
            let word: Vec<u8> = (0..len)
                .map(|i| if bits & (1 << i) == 0 { 0 } else { 255 })
                .collect();
            let reversed = word.iter().rev().copied().collect::<Vec<_>>();
            words.push((word, reversed));
        }
    }
    // Includes empty queries and spans. A tied endpoint or script can be any optimum.
    for (query, rev_query) in &words {
        for (target, rev_target) in &words {
            let row = scalar_row(query, target);
            let global = row[target.len()];
            let prefix = row.into_iter().min().unwrap();
            let suffix = scalar_row(rev_query, rev_target).into_iter().min().unwrap();
            // `prefix` covers start zero. Each remaining row scores all prefixes
            // of that suffix, including the empty one.
            let infix = (1..=target.len())
                .flat_map(|start| scalar_row(query, &target[start..]))
                .fold(prefix, usize::min);
            assert_eq!(infix_distance(query, target).unwrap(), infix);
            let check_location = |hit: crate::Location| {
                assert_eq!(hit.score, infix);
                assert_eq!(scalar(query, &target[hit.target]), infix);
            };
            check_location(infix_locate(query, target).unwrap());
            for k in 0..=query.len() {
                assert_eq!(infix_distance_at_most(query, target, k).unwrap(), (infix <= k).then_some(infix));
                let hit = infix_locate_at_most(query, target, k).unwrap();
                assert_eq!(hit.is_some(), infix <= k);
                if let Some(hit) = hit { check_location(hit); }
            }
            for (anchor, score) in [
                (Anchor::Global, global),
                (Anchor::Prefix, prefix),
                (Anchor::Suffix, suffix),
                (Anchor::Infix, infix),
            ] {
                check_paths(anchor, query, target, score, &[0, score.saturating_sub(1), score, usize::MAX]);
            }
        }
    }
}

#[test]
fn bounded_infix_prefers_the_later_exact_match() {
    let query = b"ab";
    let target = b"_ax_ab_";
    assert_eq!(infix_distance_at_most(query, target, 1).unwrap(), Some(0));
    let hit = infix_align_at_most(query, target, 1).unwrap().unwrap();
    assert_eq!(hit.target, 4..6);
    replay(Anchor::Infix, query, target, &hit, 0);
}

#[test]
fn run_cigar_formats_only_the_selected_span() {
    use crate::{Cigar, CigarFormat};

    let hit = infix_align(b"abcXef", b"__abcdef__").unwrap();
    assert_eq!(hit.target, 2..8);
    assert_eq!(hit.cigar(CigarFormat::Extended).to_string(), "3=1X2=");
    assert_eq!(Cigar::from_runs(&hit.runs, CigarFormat::Standard).to_string(), "6M");
}

#[test]
fn infix_locations_across_word_and_cutoff_boundaries() {
    for len in [63, 64, 65, 127, 128, 129, 257] {
        let query: Vec<u8> = (0..len).map(|i| (i % 251) as u8).collect();
        for stride in [1, 3, 37] {
            let mut core = query.clone();
            for b in core.iter_mut().step_by(stride) { *b = 255; }
            core.insert(len / 2, 254);
            core.remove(len / 4);
            // A poorer earlier match, then two tied copies of the edited core.
            let target = [[255; 31].as_slice(), &query[..len / 3], &[254; 67], &core, &[255; 19], &core].concat();
            for (q, t) in [(query.as_slice(), target.as_slice()), (target.as_slice(), query.as_slice())] {
                let d = infix_distance(q, t).unwrap();
                let check = |hit: crate::Location| {
                    assert_eq!(hit.score, d);
                    assert_eq!(scalar(q, &t[hit.target]), d);
                };
                check(infix_locate(q, t).unwrap());
                for k in [d.saturating_sub(1), d, usize::MAX] {
                    let hit = infix_locate_at_most(q, t, k).unwrap();
                    assert_eq!(hit.is_some(), d <= k);
                    if let Some(hit) = hit { check(hit); }
                }
            }
        }
    }
}

#[test]
fn skewed_global_paths_replay_in_both_directions_at_the_length_bound() {
    // Matching trim on both ends, several separated gap runs, and tied A's.
    let target = [b"head".as_slice(), &[b'x'; 241], b"AABA", &[b'y'; 193], b"BCz", b"tail"].concat();
    for (query, extra) in [(b"headAABCtail".as_slice(), 0), (b"headAAB\xfftail", 1)] {
        // The absent byte is reached late; it requires one substitution beyond the gap.
        let gap = target.len() - query.len();
        let score = gap + extra;
        assert_eq!(scalar(query, &target), score);
        for (q, t) in [(query, target.as_slice()), (target.as_slice(), query)] {
            check_paths(Anchor::Global, q, t, score, &[gap - 1, gap, gap + 1, usize::MAX]);
        }
    }
}

#[test]
fn alignment_trims_only_the_ends_it_pins() {
    // Both ends match globally, but retaining the match at the free end of
    // either anchored alignment would require consuming the intervening bytes.
    let query = b"BA";
    let target = b"BXBXXA";
    check_paths(Anchor::Global, query, target, 4, &[3, 4]);
    check_paths(Anchor::Prefix, query, target, 1, &[0, 1]);
    check_paths(Anchor::Suffix, query, target, 1, &[0, 1]);
}

#[test]
fn suffix_indels_replay_in_forward_order() {
    // The optimum drops Q from the query and X from the target at distinct places.
    let query = b"abcQdefgh";
    let target = b"ZZZZabcdeXfgh";
    let score = (0..=target.len()).map(|i| scalar(query, &target[i..])).min().unwrap();
    assert_eq!(score, 2);
    check_paths(Anchor::Suffix, query, target, score, &[score - 1, score]);
}

fn check_with_free_flanks(query: &[u8], core: &[u8], score: usize) {
    let flank = [b'G'; 4096];
    let target = [&flank[..], core, &flank].concat();
    let cutoffs = [score - 1, score];
    check_paths(Anchor::Global, query, core, score, &cutoffs);
    check_paths(Anchor::Prefix, query, &target[flank.len()..], score, &cutoffs);
    check_paths(Anchor::Suffix, query, &target[..flank.len() + core.len()], score, &cutoffs);
}

#[test]
fn sparse_indels_replay_across_long_free_flanks() {
    let edits = 40;
    let mut core = vec![b'A'; 128];
    for _ in 0..edits {
        core.extend_from_slice(&[b'A'; 24]);
        core.push(b'C');
    }
    core.extend_from_slice(&[b'A'; 128]);
    let query = vec![b'A'; core.len() - edits];
    // Every C needs an edit. Dropping even one C at the free end would
    // shorten the span by more than 2 * edits, leaving it too short to tie.
    check_with_free_flanks(&query, &core, edits);
}

#[test]
fn dense_edits_still_give_exact_paths_at_the_cutoff() {
    let query = [b'A'; 1665];
    let mut core = query;
    // Sparse progress before a dense patch, with shared bytes at both ends.
    for (i, byte) in core[64..query.len() - 64].iter_mut().enumerate() {
        if i % 24 == 0 || (384..1152).contains(&i) {
            *byte = b'C';
        }
    }
    // Each C requires an edit. Any omitted C also shortens the span by one
    // relative to the query, costing at least as much as consuming it.
    let score = core.iter().filter(|&&byte| byte == b'C').count();
    check_with_free_flanks(&query, &core, score);
}

#[test]
fn many_spaced_substitutions_still_give_an_exact_path() {
    // This many edits also exceeds the bidirectional checkpoint budget.
    // Handoff must still find the inclusive hit, not report a native miss.
    const EDITS: usize = 1537;
    const SPACING: usize = 8;
    let query = vec![b'A'; SPACING * (EDITS - 1) + 1];
    let mut target = query.clone();
    // Spacing keeps progress sparse; mismatches at both ends prevent trimming.
    for byte in target.iter_mut().step_by(SPACING) {
        *byte = b'C';
    }
    // Each of the EDITS non-A bytes needs an edit, and substitutions suffice.
    check_paths(Anchor::Global, &query, &target, EDITS, &[EDITS - 1, EDITS]);
}

#[test]
fn bounded_global_substitutions_and_a_cheaper_indel_path() {
    // No globally novel byte. At changed sites {0, 65} and {130, 255} are
    // disjoint, then a rotation introduces overlapping symbols and cheaper gaps.
    let q: Vec<_> = [0u8, 65, 130, 255].into_iter().cycle().take(1536).collect();
    let mut t = q.clone();
    for i in (16..1520).step_by(8) {
        t[i..i + 2].copy_from_slice(&[130, 255]);
    }
    let check = |a: &[u8], b: &[u8], d: usize| {
        for (a, b) in [(a, b), (b, a)] {
            for k in [d - 1, d] {
                assert_eq!(crate::distance_at_most(a, b, k).unwrap(), (k >= d).then_some(d));
                let path = align_at_most(a, b, k).unwrap();
                assert_eq!(path.is_some(), k >= d);
                if let Some(hit) = path { replay(Anchor::Global, a, b, &hit, d); }
            }
        }
    };
    let substitutions = q.iter().zip(&t).filter(|(x, y)| x != y).count();
    check(&q, &t, substitutions);
    t[512..1312].rotate_left(1);
    let d = scalar(&q, &t);
    assert!(d < q.iter().zip(&t).filter(|(x, y)| x != y).count());
    check(&q, &t, d);
}

#[test]
fn inconclusive_global_proof_can_fall_back_to_an_exact_path() {
    let mut query = [b'A'; 1536];
    let mut target = query;
    target[..240].fill(b'C');
    target[1296..].fill(b'C');
    // Both mismatch directions occur immediately, and dense edits at both
    // ends leave too little progress for native checkpoint traceback.
    query[0] = b'C';
    target[0] = b'A';
    let query = [b"head".as_slice(), &query, b"tail"].concat();
    let target = [b"head".as_slice(), &target, b"tail"].concat();
    let score = scalar(&query, &target);
    check_paths(Anchor::Global, &query, &target, score, &[score - 1, score]);
}

#[test]
fn mixed_indels_replay_with_slack_at_either_anchor() {
    let mut query = Vec::new();
    let mut core = Vec::new();
    let mut state = 0x0866_8df1_8276_773b_u64;
    for i in 0..288 {
        // The first and second halves drift in opposite directions. Each
        // edit is separated by a variable-length, pseudorandom match stretch.
        match (i < 144, i % 4) {
            (true, 0 | 1) | (false, 0) => query.push(250),
            (true, 2) | (false, 1 | 2) => core.push(251),
            _ => {
                query.push(250);
                core.push(251);
            }
        }
        if i == 287 {
            break;
        }
        for _ in 0..(15 + i % 5) {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let byte = (state >> 57) as u8;
            query.push(byte);
            core.push(byte);
        }
    }

    let row = scalar_row(&query, &core);
    let global = row[core.len()];
    let prefix = row.iter().copied().min().unwrap();
    // Ending before or after the last target marker ties, with no preference
    // required of the public path API.
    assert_eq!((row[core.len() - 1], global), (prefix, prefix));
    let rev_query: Vec<_> = query.iter().rev().copied().collect();
    let rev_core: Vec<_> = core.iter().rev().copied().collect();
    let suffix = scalar_row(&rev_query, &rev_core).into_iter().min().unwrap();
    // No flank byte can match. Removing a free-end flank column either
    // discards a deletion or exchanges a substitution for an insertion,
    // so the scalar row over the core suffices for the flanked oracle.
    let flank = [255; 8192];
    let prefix_target = [&core[..], &flank].concat();
    let suffix_target = [&flank[..], &core].concat();
    for (anchor, target, score) in [
        (Anchor::Global, &core[..], global),
        (Anchor::Prefix, &prefix_target, prefix),
        (Anchor::Suffix, &suffix_target, suffix),
    ] {
        assert!(score > 256, "{anchor:?}: {score}");
        check_paths(anchor, &query, target, score, &[score - 1, score, score + 37, usize::MAX]);
    }
}

#[test]
fn anchored_hit_can_use_the_last_byte_allowed_by_the_cutoff() {
    let query = b"0123456789abcdefGHIJKLMNOPQRSTUV";
    let core = b"0123456789abcdefXYZGHIJKLMNOPQRSTUV";
    let flank = [255; 4096];
    let prefix_target = [&core[..], &flank].concat();
    let suffix_target = [&flank[..], core].concat();
    // The three nonmatching bytes occur too far from either edge to avoid
    // with fewer edits. At that score the whole core is the unique span.
    let score = core.len() - query.len();
    let rev_query: Vec<_> = query.iter().rev().copied().collect();
    let rev_core: Vec<_> = core.iter().rev().copied().collect();
    for (anchor, target, row) in [
        (Anchor::Prefix, &prefix_target, scalar_row(query, core)),
        (Anchor::Suffix, &suffix_target, scalar_row(&rev_query, &rev_core)),
    ] {
        assert_eq!(row[core.len()], score);
        assert!(row[..core.len()].iter().all(|&earlier| earlier > score));
        check_paths(anchor, query, target, score, &[score - 1, score, score + 1, usize::MAX]);
    }
}

#[test]
fn one_edit_above_the_length_gap_can_choose_a_shorter_span() {
    // The full target needs four edits; dropping the free-end C needs just three.
    check_paths(Anchor::Prefix, b"xxxAB", b"ABC", 3, &[2, 3]);
    check_paths(Anchor::Suffix, b"BAxxx", b"CBA", 3, &[2, 3]);
}

#[test]
fn high_score_two_run_anchored_path_does_not_reserve_per_edit() {
    const EDITS: usize = 512;
    let mut query = [b'C'; EDITS + 1];
    query[EDITS] = b'A';
    // At the exact gap only the whole target fits; one extra edit also admits
    // the empty span. Both searches must return the compact two-run optimum.
    for k in [EDITS, EDITS + 1] {
        let hit = prefix_align_at_most(&query, b"A", k).unwrap().unwrap();
        replay(Anchor::Prefix, &query, b"A", &hit, EDITS);
        assert_eq!(hit.runs.len(), 2);
        assert!(hit.runs.capacity() < EDITS / 8, "retained capacity should follow runs, not edits");
    }
}
