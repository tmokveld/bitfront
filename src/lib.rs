//! Exact unit-cost Levenshtein scores on bytes. Global distance consumes both
//! inputs; prefix and suffix distances consume the query and one anchored end
//! of the target. No alignment path or endpoint list is constructed.
#![forbid(unsafe_code)]

mod flank;
mod kernel;
mod word;

/// Required score dimensions or buffers cannot be accommodated.
pub use rust_edlib::AlignError as Error;

/// Returns the minimum number of single-byte edits transforming `a` into `b`.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required Myers dimensions or buffer
/// sizes cannot be represented. Exceeding the wavefront budget is not an error.
///
/// ```
/// assert_eq!(bitfront::distance(b"kitten", b"sitting")?, 3);
/// assert_eq!(bitfront::distance(b"\x00\xff", b"\xff")?, 1);
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn distance(a: &[u8], b: &[u8]) -> Result<usize, Error> {
    let (a, b) = kernel::trim(a, b);
    kernel::hybrid(a, b, kernel::cap(a.len(), b.len()))
}

/// Returns the exact global edit distance if it is at most `k` (inclusive).
/// `Ok(None)` proves the distance exceeds `k`. Unlike [`distance`], this
/// search stops at the bound and never restarts an unbounded scorer.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the remaining wavefront coordinates or
/// row storage cannot be accommodated. This is not a bound miss.
///
/// ```
/// assert_eq!(bitfront::distance_at_most(b"kitten", b"sitting", 2)?, None);
/// assert_eq!(bitfront::distance_at_most(b"kitten", b"sitting", 3)?, Some(3));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn distance_at_most(a: &[u8], b: &[u8], k: usize) -> Result<Option<usize>, Error> {
    if k == 0 {
        return Ok((a == b).then_some(0));
    }
    // Trimming cannot change the length gap. Reject before scanning either end.
    if a.len().abs_diff(b.len()) > k {
        return Ok(None);
    }
    let (a, b) = kernel::trim(a, b);
    kernel::distance_at_most(a, b, k)
}

/// Minimum global edit distance from the entire `query` to any prefix of
/// `target`, including its empty prefix. The unconsumed target suffix is free.
/// Swapping the arguments changes the question.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required Myers dimensions or buffers cannot be accommodated.
///
/// ```
/// assert_eq!(bitfront::prefix_distance(b"CAGCAG", b"CAGCAGCAA")?, 0);
/// assert_eq!(bitfront::prefix_distance(b"CAGCAGCAA", b"CAGCAG")?, 3);
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn prefix_distance(query: &[u8], target: &[u8]) -> Result<usize, Error> {
    flank::with_cap::<false>(query, target, flank::CAP)
}

/// Returns the exact [`prefix_distance`] if it is at most `k` (inclusive).
/// `Ok(None)` means every target prefix, including the empty one, costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required scoring dimensions or
/// buffers cannot be accommodated.
///
/// ```
/// assert_eq!(bitfront::prefix_distance_at_most(b"ACCT", b"ACGTTAGC", 1)?, Some(1));
/// assert_eq!(bitfront::prefix_distance_at_most(b"ACCT", b"ACGTTAGC", 0)?, None);
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn prefix_distance_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<usize>, Error> {
    flank::at_most::<false>(query, target, k, flank::CAP)
}

/// Minimum global edit distance from the entire `query` to any suffix of
/// `target`, including its empty suffix. The unconsumed target prefix is free.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required Myers dimensions or buffers cannot be accommodated.
///
/// ```
/// assert_eq!(bitfront::suffix_distance(b"CAA", b"CAGCAGCAA")?, 0);
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn suffix_distance(query: &[u8], target: &[u8]) -> Result<usize, Error> {
    flank::with_cap::<true>(query, target, flank::CAP)
}

/// Returns the exact [`suffix_distance`] if it is at most `k` (inclusive).
/// `Ok(None)` means every target suffix, including the empty one, costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required scoring dimensions or
/// buffers cannot be accommodated.
///
/// ```
/// assert_eq!(bitfront::suffix_distance_at_most(b"TAGC", b"ACGTTAGC", 0)?, Some(0));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn suffix_distance_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<usize>, Error> {
    flank::at_most::<true>(query, target, k, flank::CAP)
}

#[cfg(test)]
mod tests;
