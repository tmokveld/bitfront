use crate::{Error, Location, kernel};

// First progress check, not a hard edit limit.
pub(crate) const CAP: usize = 16;

pub(crate) fn trimmed_window<'q, 't, const REVERSE: bool>(
    query: &'q [u8], target: &'t [u8], k: usize,
) -> (&'q [u8], &'t [u8]) {
    let matched = kernel::matching_words::<REVERSE>(query, target);
    let n = query.len() - matched;
    let m = target.len() - matched;
    // Deleting the remaining query costs n. Beyond n + min(k, n) target
    // bytes, the length difference alone exceeds k or that deletion cost.
    let extent = m.min(n.saturating_add(k.min(n)));
    if REVERSE {
        (&query[..n], &target[m - extent..m])
    } else {
        (&query[matched..], &target[matched..matched + extent])
    }
}

/// Score the whole query against a prefix (or, in reverse, suffix). Never swap
/// query and target: only the unpinned target end is free. A positive cap
/// starts the density checks; zero goes directly to Myers.
pub(crate) fn with_cap<const REVERSE: bool>(query: &[u8], target: &[u8], cap: usize) -> Result<usize, Error> {
    Ok(at_most::<REVERSE>(query, target, query.len(), cap)?
        .expect("the empty flank is within the query-length bound"))
}

pub(crate) fn at_most<const REVERSE: bool>(
    query: &[u8], target: &[u8], k: usize, checkpoint: usize,
) -> Result<Option<usize>, Error> {
    search::<REVERSE, _>(query, target, k, checkpoint, |score, _| score, |q, t, _| {
        if REVERSE {
            rust_edlib::suffix_distance_at_most(q, t, k)
        } else {
            rust_edlib::prefix_distance_at_most(q, t, k)
        }
    })
}

pub(crate) fn locate_at_most<const REVERSE: bool>(
    query: &[u8], target: &[u8], k: usize,
) -> Result<Option<Location>, Error> {
    let n = target.len();
    let location = |score, consumed| Location {
        score,
        target: if REVERSE { n - consumed..n } else { 0..consumed },
    };
    search::<REVERSE, _>(query, target, k, CAP, location, |q, t, pinned| {
        let hit = if REVERSE {
            rust_edlib::suffix_locate_at_most(q, t, k)
        } else {
            rust_edlib::prefix_locate_at_most(q, t, k)
        }?;
        Ok(hit.map(|h| location(h.score, pinned + h.target.len())))
    })
}

fn search<const REVERSE: bool, T>(
    query: &[u8], target: &[u8], k: usize, checkpoint: usize,
    hit: impl FnOnce(usize, usize) -> T,
    myers: impl FnOnce(&[u8], &[u8], usize) -> Result<Option<T>, Error>,
) -> Result<Option<T>, Error> {
    let gap = query.len().saturating_sub(target.len());
    if gap > k {
        return Ok(None);
    }
    if gap != 0 && k == gap {
        // A shorter span would already exceed k on length alone.
        return Ok(rust_edlib::distance_at_most(query, target, k)?.map(|d| hit(d, target.len())));
    }
    let (q, t) = trimmed_window::<REVERSE>(query, target, k);
    let pinned = query.len() - q.len();
    if q.is_empty() || t.is_empty() {
        return Ok(Some(hit(q.len(), pinned)));
    }
    if k == 0 {
        return Ok(None); // The first remaining bytes at the pinned end differ.
    }
    match kernel::anchored::<REVERSE>(q, t, k, checkpoint) {
        kernel::AnchorProbe::Hit { score, consumed } => Ok(Some(hit(score as usize, pinned + consumed as usize))),
        kernel::AnchorProbe::Miss => Ok(None),
        // A speculative probe has not certified a starting cell or spent any of the cutoff.
        kernel::AnchorProbe::Myers => myers(q, t, pinned),
    }
}
