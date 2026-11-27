use bitfront_bench::{Mode, run_with_revision};

fn main() {
    run_with_revision(env!("BITFRONT_PREVIOUS_COMMIT"), |mode, query, target, cutoff| {
        match (mode, cutoff) {
            (Mode::Global, None) => previous::distance(query, target).map(Some),
            (Mode::Global, Some(k)) => previous::distance_at_most(query, target, k),
            (Mode::Prefix, None) => previous::prefix_distance(query, target).map(Some),
            (Mode::Prefix, Some(k)) => previous::prefix_distance_at_most(query, target, k),
            (Mode::Suffix, None) => previous::suffix_distance(query, target).map(Some),
            (Mode::Suffix, Some(k)) => previous::suffix_distance_at_most(query, target, k),
        }.expect("previous Bitfront failed to score benchmark input")
    });
}
