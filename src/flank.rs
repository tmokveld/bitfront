use crate::{Error, kernel};

// First progress check, not a hard edit limit.
pub(crate) const CAP: usize = 16;

fn trimmed_window<'q, 't, const REVERSE: bool>(
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
    // Excess target bytes are free; excess query bytes must be deleted.
    if query.len().saturating_sub(target.len()) > k {
        return Ok(None);
    }
    let (query, target) = trimmed_window::<REVERSE>(query, target, k);
    if query.is_empty() || target.is_empty() {
        return Ok(Some(query.len()));
    }
    if k == 0 {
        return Ok(None); // The first remaining bytes at the pinned end differ.
    }
    match kernel::anchored::<REVERSE>(query, target, k, checkpoint) {
        kernel::AnchorProbe::Hit(d) => Ok(Some(d)),
        kernel::AnchorProbe::Miss => Ok(None),
        // Restart at the pinned end with the caller's cutoff, not k minus the
        // probe cost. No particular path or endpoint was certified by the probe.
        kernel::AnchorProbe::Myers if REVERSE => rust_edlib::suffix_distance_at_most(query, target, k),
        kernel::AnchorProbe::Myers => rust_edlib::prefix_distance_at_most(query, target, k),
    }
}
