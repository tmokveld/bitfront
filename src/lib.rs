//! Exact unit-cost Levenshtein scores and alignments on bytes. All modes consume
//! the entire query. Global alignment consumes the entire target, prefix and
//! suffix alignment pin one end, and infix alignment allows any contiguous
//! target substring, including an empty one. Use [`distance`], [`prefix_distance`],
//! [`suffix_distance`], or [`infix_distance`] for a score without constructing a
//! path or endpoint list. [`prefix_locate`], [`suffix_locate`], and [`infix_locate`]
//! also return one optimal target range without constructing an edit path.
#![forbid(unsafe_code)]

mod alignment;
mod flank;
mod kernel;

pub use alignment::{
    Alignment, Op, Run, align, align_at_most, infix_align, infix_align_at_most, prefix_align,
    prefix_align_at_most, suffix_align, suffix_align_at_most,
};

/// Required alignment dimensions or buffers cannot be accommodated.
pub use rust_edlib::AlignError as Error;

pub use rust_edlib::{
    Cigar, CigarFormat, Location, infix_distance, infix_distance_at_most, infix_locate,
    infix_locate_at_most,
};

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
    kernel::hybrid(a, b)
}

/// Returns the exact global edit distance if it is at most `k` (inclusive).
/// `Ok(None)` proves the distance exceeds `k`. At the length gap, checking
/// whether the shorter input is a subsequence decides both hits and misses.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if the required scoring coordinates or
/// buffers cannot be accommodated. This is not a bound miss.
///
/// ```
/// assert_eq!(bitfront::distance_at_most(b"kitten", b"sitting", 2)?, None);
/// assert_eq!(bitfront::distance_at_most(b"kitten", b"sitting", 3)?, Some(3));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn distance_at_most(a: &[u8], b: &[u8], k: usize) -> Result<Option<usize>, Error> {
    let gap = a.len().abs_diff(b.len());
    if k < gap {
        return Ok(None);
    }
    if k == gap {
        return rust_edlib::distance_at_most(a, b, k);
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

/// Returns the score and one optimal target prefix as a half-open range,
/// without traceback. The range can be empty; tied ranges are not canonical.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn prefix_locate(query: &[u8], target: &[u8]) -> Result<Location, Error> {
    Ok(prefix_locate_at_most(query, target, query.len())?.expect("the empty prefix is within the bound"))
}

/// Returns one [`prefix_locate`] hit at or below `k` (inclusive), or `None` if every prefix costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
///
/// ```
/// let hit = bitfront::prefix_locate_at_most(b"ACGT", b"ATCGTzz", 1)?.unwrap();
/// assert_eq!((hit.score, hit.target), (1, 0..5));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn prefix_locate_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Location>, Error> {
    flank::locate_at_most::<false>(query, target, k)
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

/// Returns the score and one optimal target suffix in forward, half-open byte
/// coordinates, without traceback. The range can be empty; ties are not canonical.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
pub fn suffix_locate(query: &[u8], target: &[u8]) -> Result<Location, Error> {
    Ok(suffix_locate_at_most(query, target, query.len())?.expect("the empty suffix is within the bound"))
}

/// Returns one [`suffix_locate`] hit at or below `k` (inclusive), or `None` if every suffix costs more.
///
/// # Errors
/// Returns [`Error::InputTooLong`] if required dimensions or buffers cannot be accommodated.
///
/// ```
/// let hit = bitfront::suffix_locate_at_most(b"TAGC", b"ACGTTAGC", 0)?.unwrap();
/// assert_eq!((hit.score, hit.target), (0, 4..8));
/// # Ok::<(), bitfront::Error>(())
/// ```
pub fn suffix_locate_at_most(query: &[u8], target: &[u8], k: usize) -> Result<Option<Location>, Error> {
    flank::locate_at_most::<true>(query, target, k)
}

#[cfg(test)]
mod tests;
