use super::{
    distance, distance_at_most, flank, kernel, prefix_distance, prefix_distance_at_most,
    suffix_distance, suffix_distance_at_most,
};

// Wagner–Fischer prefix DP, independent of both the diagonal frontier and Myers.
fn scalar_row(a: &[u8], b: &[u8]) -> Vec<usize> {
    let mut row: Vec<_> = (0..=b.len()).collect();
    for (i, x) in a.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, y) in b.iter().enumerate() {
            let above = row[j + 1];
            row[j + 1] = (diagonal + usize::from(x != y)).min(above.min(row[j]) + 1);
            diagonal = above;
        }
    }
    row
}

fn scalar(a: &[u8], b: &[u8]) -> usize {
    scalar_row(a, b)[b.len()]
}

fn check(a: &[u8], b: &[u8], d: usize) {
    assert_eq!(distance(a, b), Ok(d));
    if let Some(k) = d.checked_sub(1) {
        assert_eq!(distance_at_most(a, b, k), Ok(None));
    }
    assert_eq!(distance_at_most(a, b, d), Ok(Some(d)));
    assert_eq!(distance_at_most(a, b, d + 1), Ok(Some(d)));
}

fn check_anchor<const REVERSE: bool>(query: &[u8], target: &[u8], d: usize) {
    let public = if REVERSE { suffix_distance } else { prefix_distance };
    assert_eq!(public(query, target), Ok(d));
    for cap in [0, d.saturating_sub(1), d, d + 1] {
        assert_eq!(flank::with_cap::<REVERSE>(query, target, cap), Ok(d));
    }
    let at_most = if REVERSE { suffix_distance_at_most } else { prefix_distance_at_most };
    for k in [0, d.saturating_sub(1), d, d + 1, usize::MAX] {
        assert_eq!(at_most(query, target, k), Ok((d <= k).then_some(d)), "reverse={REVERSE}, k={k}");
    }
}

#[test]
fn dense_flanks_keep_the_full_cutoff_after_the_probe() {
    let query = [b'A'; 129];
    let bad = [b'C'; 24];
    let good = [b'A'; 105];
    let free = [b'G'; 80];
    // Only 105 query bytes can match. Substituting the other 24 attains the
    // bound, even though the first 16 edit levels cannot consume the query.
    let left = [bad.as_slice(), &good, &free].concat();
    let right = [free.as_slice(), &good, &bad].concat();
    check_anchor::<false>(&query, &left, 24);
    check_anchor::<true>(&query, &right, 24);
}

#[test]
fn bounded_flanks_cover_full_and_clipped_inline_rows() {
    for d in [15, 16, 17] {
        let mut query = [b'A'; 96];
        query[..d].fill(b'C');
        let mut reverse = query;
        reverse.reverse();
        // Each C must be edited. Vary substitutions vs. deletions to cover
        // a single diagonal, a width that plateaus, and the full 2k+1 strip.
        for substituted in [0, 1, d / 2, d] {
            let mut target = vec![b'G'; substituted];
            target.extend_from_slice(&query[d..]);
            if substituted == d {
                target.extend_from_slice(&[b'T'; 32]);
            }
            check_anchor::<false>(&query, &target, d);
            target.reverse();
            check_anchor::<true>(&reverse, &target, d);
        }
    }
}

fn check_anchors(query: &[u8], target: &[u8]) {
    // Minimum in the last DP row includes the empty target prefix. Reverse
    // copies here, in the oracle only, make its suffix implementation simple.
    let prefix = scalar_row(query, target).into_iter().min().unwrap();
    let rev = |s: &[u8]| s.iter().rev().copied().collect::<Vec<_>>();
    let suffix = scalar_row(&rev(query), &rev(target)).into_iter().min().unwrap();
    check_anchor::<false>(query, target, prefix);
    check_anchor::<true>(query, target, suffix);
}

#[test]
fn exhaustive_binary_strings_through_six_bytes() {
    let mut words = Vec::new();
    for len in 0..=6 {
        for bits in 0..(1 << len) {
            let word = (0..len)
                .map(|i| if bits & (1 << i) == 0 { 0 } else { 255 })
                .collect::<Vec<_>>();
            words.push(word);
        }
    }
    for a in &words {
        for b in &words {
            check(a, b, scalar(a, b));
            check_anchors(a, b);
        }
    }
}

#[test]
fn word_boundaries_and_unequal_patterns_in_both_orders() {
    for m in [63, 64, 65] {
        let a: Vec<_> = (0..m).map(|i| (17 * i) as u8).collect();
        for n in [63, 64, 65] {
            let mut b: Vec<_> = (0..n).map(|i| (17 * i) as u8).collect();
            for i in (0..n - 1).step_by(3) {
                b[i] ^= 0xff;
            }
            b[n - 1] ^= 0xff; // Keep the word boundary intact after trimming.
            let d = scalar(&a, &b);
            check(&a, &b, d);
            check(&b, &a, d);
        }
    }
}

#[test]
fn global_extension_after_a_gap_hits_each_word_lane_and_tail() {
    let a: Vec<_> = (0..33).map(|i| (37 * i) as u8).collect();
    for mismatch in 0..a.len() - 1 {
        let mut b = vec![0xff];
        b.extend_from_slice(&a);
        b[mismatch + 1] ^= 0x80;
        b[a.len()] ^= 1; // Keep the right edge in the wavefront, not in trim.
        let d = scalar(&a, &b);
        check(&a, &b, d);
        check(&b, &a, d);
    }
}

#[test]
fn trim_shared_ends_without_overlap() {
    let flank: Vec<u8> = (0..=255).cycle().take(4097).collect();
    let pad = |middle: &[u8]| [flank.as_slice(), middle, flank.as_slice()].concat();
    let a = pad(b"\0abababab\xff");
    let b = pad(b"\xffbabababa\0");
    let d = scalar(b"\0abababab\xff", b"\xffbabababa\0");
    check(&a, &b, d);
    check(&a, &a, 0);
    check(&a, &a[1..], 1);
    check(b"aaaa", b"aaaaaaaaa", 5);
    assert_eq!(distance_at_most(&a, &b, usize::MAX), Ok(Some(d)));

    // Unequal remainders, including a partial reverse word on either side.
    let (x, y) = (&b"\x80"[..], &b"\x81\xfe"[..]);
    for tail in 0..=16 {
        let padded = |middle: &[u8]| [&flank[..3], middle, &flank[..tail]].concat();
        let (a, b) = (padded(x), padded(y));
        assert_eq!(kernel::trim(&a, &b), (x, y));
        assert_eq!(kernel::trim(&b, &a), (y, x));
        check(&a, &b, 2);
        check(&b, &a, 2);
    }
}

#[test]
fn public_answers_after_hit_miss_and_length_gap_skip() {
    let a: Vec<_> = (0..=255).collect();
    let mut near = a.clone();
    near[0] ^= 1;
    near[255] ^= 1;
    let k = kernel::cap(a.len(), near.len());
    assert_eq!(distance_at_most(&a, &near, k), Ok(Some(2)));
    assert_eq!(distance(&a, &near), Ok(2));

    // Interior matches are tempting false restart points on the losing probe.
    let mut far = a.clone();
    far[..96].fill(255);
    far[160..].fill(0);
    assert_eq!(distance_at_most(&a, &far, k), Ok(None));
    check(&a, &far, scalar(&a, &far));

    let short = [0xff; 65];
    let long = [0; 192];
    assert!(long.len() - short.len() > kernel::cap(short.len(), long.len()));
    check(&short, &long, long.len());
    check(&long, &short, long.len());
}

#[test]
fn large_gap_needs_an_ordered_subsequence() {
    let short: Vec<u8> = (0..80).map(|i| 255 - i).collect();
    let mut long: Vec<u8> = short.iter().flat_map(|&b| [0, b, 1]).collect();
    for swapped in [false, true] {
        if swapped {
            // Same byte counts, but these two unique bytes no longer occur in order.
            long.swap(1 + 3 * 39, 1 + 3 * 40);
        }
        let d = scalar(&short, &long);
        assert_eq!(d == long.len() - short.len(), !swapped);
        check(&short, &long, d);
        check(&long, &short, d);
    }
}

#[test]
fn near_miss_path_leaves_the_first_endpoint_band() {
    let a: Vec<u8> = (0..192).collect();
    let mut b = [&a[..80], &[255; 9], &a[80..160], &a[169..]].concat();
    b[0] = 254;
    b[191] = 253;
    // The path crosses diagonal 9 after score 9. It is outside the
    // cap-16 band there, but competitive in the larger search.
    let d = scalar(&a, &b);
    assert_eq!(d, 20);
    check(&a, &b, d);
    check(&b, &a, d);
}

#[test]
fn bounded_scores_for_candidate_overlaps() {
    // Successive suffix/prefix windows, as in TRGT's overlap search. A miss
    // rules out the current window, not a shorter one with different ends.
    let (left, right) = (b"ACGT", b"CTAA");
    let score = |len, k| distance_at_most(&left[left.len() - len..], &right[..len], k);
    assert_eq!(score(4, 3), Ok(None));
    assert_eq!(score(3, 3), Ok(Some(2)));
    assert_eq!(score(2, 1), Ok(Some(1)));
}

#[test]
fn bounded_global_rows_expand_spill_and_shrink() {
    let a: Vec<_> = (0..128).map(|i| (29 * i) as u8).collect();
    for gap in [0, 7, 31] {
        let mut b = a.clone();
        for x in b.iter_mut().step_by(3) {
            *x ^= 0x80;
        }
        b[127] ^= 0x40;
        b.splice(57..57, (0..gap).map(|i| (157 * i) as u8));
        let d = scalar(&a, &b);
        // A tight bound spills both ping-pong rows at scores 17 and 18,
        // then puts both back inline before reaching the endpoint.
        assert!(d - gap >= 36);
        for (x, y) in [(&a, &b), (&b, &a)] {
            check(x, y, d);
            assert_eq!(distance_at_most(x, y, usize::MAX), Ok(Some(d)));
            // The corresponding anchored optima can end away from `m - n`.
            check_anchors(x, y);
        }
    }
}

#[test]
fn wavefront_coordinate_limit_without_giant_allocations() {
    let max = i32::MAX as usize;
    assert!(kernel::coordinates_fit(max - 2, 1));
    assert!(!kernel::coordinates_fit(max - 1, 1));
    assert!(!kernel::coordinates_fit(usize::MAX, 1));
}

#[test]
fn only_the_pinned_end_can_be_trimmed() {
    // Stripping the shared byte at the *free* end would change 1 to 0.
    assert_eq!(prefix_distance(b"BA", b"BXBA"), Ok(1));
    assert_eq!(suffix_distance(b"AB", b"ABXB"), Ok(1));
    check_anchors(b"ABC", b"AXBCxxxxxxxx");
    check_anchors(b"xxxxxABCxxxxxxxx", b"AXBC");
    let repeat: Vec<u8> = [0, 255, 1].into_iter().cycle().take(4097).collect();
    let observation = [&repeat[..], b"xy"].concat();
    let allele = [&repeat[..], b"xBy", &repeat[..]].concat();
    assert_eq!(prefix_distance(&observation, &allele), Ok(1));
    let reverse = |s: &[u8]| s.iter().rev().copied().collect::<Vec<_>>();
    assert_eq!(suffix_distance(&reverse(&observation), &reverse(&allele)), Ok(1));
}

#[test]
fn flank_optima_cross_words_and_probe_budgets() {
    for n in [63, 64, 65, 127, 128, 129, 256] {
        let query: Vec<u8> = (0..n).map(|i| (29 * i) as u8).collect();
        let mut target = query.clone();
        target.insert(n / 2, 0xff);
        target[0] ^= 0x80;
        target[n] ^= 0x40;
        let near = [&target[..], b"unobserved"].concat();
        check_anchors(&query, &near);
        for x in target.iter_mut().step_by(2) { *x ^= 0xff; }
        check_anchors(&query, &target);
        check_anchors(&query, &target[..n / 2]);
    }
}

#[test]
fn global_fronts_meet_across_diagonals_and_score_parities() {
    // No shared ends to trim. Long insertions and deletions on opposite
    // sides put the join away from both end diagonals, in spilled rows.
    let a: Vec<u8> = (0..600).map(|i| ((i * 37 + i / 7) % 251) as u8).collect();
    for inserted in [30, 31, 37] {
        let mut b = [&a[..140], &[255; 37][..inserted], &a[140..430], &a[461..]].concat();
        b[0] = 253;
        *b.last_mut().unwrap() = 254;
        let d = scalar(&a, &b);
        assert_eq!(d, 2 + inserted + 31);
        check(&a, &b, d);
        check(&b, &a, d);
    }
}

#[test]
fn sparse_anchored_indels_reach_a_free_endpoint_after_the_initial_probe() {
    let query: Vec<u8> = (0..600).map(|i| (i % 251) as u8).collect();
    let mut core = query.clone();
    for i in (0..600).step_by(21).rev() {
        core.insert(i, 255);
    }
    let edge = [254; 83];
    let prefix = [core.as_slice(), &edge].concat();
    let suffix = [&edge, core.as_slice()].concat();
    check_anchors(&query, &prefix);
    check_anchors(&query, &suffix);
}

#[test]
fn one_word_midpoint_can_be_at_either_query_edge() {
    for n in [2, 31, 32, 63, 64] {
        let query: Vec<_> = (0..n).map(|i| (199 * i) as u8).collect();
        for extra in [n + 2, n + 3] {
            let pad = vec![255; extra];
            for mut target in [[&pad[..], &query].concat(), [&query[..], &pad].concat()] {
                let pos = target.iter().position(|&b| b == query[n / 2]).unwrap();
                target[pos] = 254;
                let d = scalar(&query, &target);
                assert_eq!(super::word::upper_bound(&query, &target), d);
                assert_eq!(distance(&target, &query), Ok(d));
            }
        }
    }
}

#[test]
fn sliding_word_band_certifies_its_edges_but_not_a_large_excursion() {
    let a: Vec<u8> = (0..265).map(|i| ((37 * i + i / 7) % 251) as u8).collect();
    for n in [65, 72, 127, 193, 265] {
        let a = &a[..n];
        for gap in [0, 1, 7, 28, 55, 56] {
            // An optimal path can spend all its slack on an excursion to
            // either extreme diagonal, including in the partial last tile.
            let deleted = (super::word::BAND - gap) / 2;
            let inserted = gap + deleted;
            for near_start in [false, true] {
                let mut b = a.to_vec();
                if near_start {
                    b.drain(1..1 + deleted);
                    b.splice(b.len() - 1..b.len() - 1, std::iter::repeat_n(255, inserted));
                } else {
                    b.drain(n - 1 - deleted..n - 1);
                    b.splice(1..1, std::iter::repeat_n(255, inserted));
                }
                let d = scalar(a, &b);
                assert_eq!(super::word::upper_bound(a, &b), d);
                assert_eq!(distance(&b, a), Ok(d));
            }
        }
    }

    let mut b = [&a[..36], &[255; 48], &a[36..198], &a[246..]].concat();
    b[..32].fill(254);
    b[240..].fill(253);
    let d = scalar(&a, &b);
    assert!(super::word::upper_bound(&a, &b) > d);
    assert_eq!(distance(&a, &b), Ok(d));
    assert_eq!(distance(&b, &a), Ok(d));
}

#[test]
fn dense_handoff_does_not_treat_a_missed_prediction_as_a_score() {
    let a = [0; 2048];
    let mut b = a;
    for (i, byte) in b.iter_mut().enumerate() {
        if (320..1728).contains(&i) || i % 8 == 0 || i == 2047 {
            *byte = 255;
        }
    }
    // Every nonzero byte requires one edit, and substitutions suffice.
    let d = b.iter().filter(|&&byte| byte != 0).count();
    assert_eq!(distance(&a, &b), Ok(d));
    assert_eq!(distance(&b, &a), Ok(d));
}
