mod heap;
mod memory;
mod native;
mod options;
mod path;
mod timing;
mod workload;

use criterion::{Criterion, Throughput};
use native::Wfa2;
use options::{Options, USAGE};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::PathBuf;
use std::time::Duration;
use timing::{Pair, measure, median_per_pair};
use workload::{Case, Scenario};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Task {
    Score,
    Path,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Global,
    Prefix,
    Suffix,
}

impl Mode {
    fn name(self) -> &'static str {
        match self { Self::Global => "global", Self::Prefix => "prefix", Self::Suffix => "suffix" }
    }

    fn bitfront(self, query: &[u8], target: &[u8], cutoff: Option<usize>) -> Option<usize> {
        match (self, cutoff) {
            (Self::Global, None) => bitfront::distance(query, target).map(Some),
            (Self::Global, Some(k)) => bitfront::distance_at_most(query, target, k),
            (Self::Prefix, None) => bitfront::prefix_distance(query, target).map(Some),
            (Self::Prefix, Some(k)) => bitfront::prefix_distance_at_most(query, target, k),
            (Self::Suffix, None) => bitfront::suffix_distance(query, target).map(Some),
            (Self::Suffix, Some(k)) => bitfront::suffix_distance_at_most(query, target, k),
        }
        .expect("bitfront failed to score benchmark input")
    }
}

struct Markdown {
    home: PathBuf,
    out: BufWriter<File>,
}

fn case_pairs(case: &Case, count: usize) -> Vec<Pair> {
    (0..count).map(|i| {
        let (query, target) = case.pair(i);
        let distance = case.mode.rust_edlib(&query, &target, None).expect("an unbounded alignment has a score");
        Pair { query, target, distance, cutoff: case.scenario.cutoff(distance) }
    }).collect()
}

const CASE_COLUMNS: &str = " Search | Template bp | Query bp | Target bp | Nominal edit rate | Exact distance | Global extra edits |";
const CASE_ALIGNMENT: &str = "---|---:|---:|---:|---:|---:|---:|";

// Timing and memory tables must identify the same panel, including hits at an inclusive cutoff.
fn write_case(out: &mut impl Write, case: &Case, pairs: &[Pair]) -> io::Result<()> {
    let [mut query, mut target, mut distance, mut extra, mut cutoff] = [(usize::MAX, 0); 5];
    let mut hits = 0;
    let global = case.mode == Mode::Global;
    let include = |(low, high): &mut (usize, usize), value: usize| {
        *low = (*low).min(value);
        *high = (*high).max(value);
    };
    for pair in pairs {
        let (q, t, d) = (pair.query.len(), pair.target.len(), pair.distance);
        include(&mut query, q);
        include(&mut target, t);
        include(&mut distance, d);
        if global {
            include(&mut extra, d - q.abs_diff(t));
        }
        if let Some(k) = pair.cutoff {
            include(&mut cutoff, k);
            hits += usize::from(d <= k);
        }
    }
    write!(out, "| {}", case.mode.name())?;
    if let Some(name) = case.scenario.name() {
        write!(out, " ({name})")?;
    }
    if case.swap {
        write!(out, " (swapped)")?;
    }
    if cutoff.0 <= cutoff.1 {
        write!(out, " ≤ ")?;
        write_bounds(out, cutoff)?;
        write!(out, " ({hits}/{} hit)", pairs.len())?;
    }
    write!(out, " | {} |", case.length)?;
    write_range(out, query)?;
    write_range(out, target)?;
    write!(out, " {:.2}% |", f64::from(case.rate_bp) / 100.)?;
    write_range(out, distance)?;
    if global { write_range(out, extra) } else { write!(out, " — |") }
}

fn write_range(out: &mut impl Write, bounds: (usize, usize)) -> io::Result<()> {
    write!(out, " ")?;
    write_bounds(out, bounds)?;
    write!(out, " |")
}

fn write_bounds(out: &mut impl Write, (low, high): (usize, usize)) -> io::Result<()> {
    if low == high {
        write!(out, "{low}")
    } else {
        write!(out, "{low}–{high}")
    }
}

fn compare(
    criterion: &mut Criterion, report: Option<&mut Markdown>, case: &Case, pairs: &[Pair], task: Task,
    previous: Option<(&str, &impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>)>,
) -> io::Result<()> {
    let count = pairs.len();
    let mut wfa2 = Wfa2::new(case.mode, pairs[0].cutoff, task);
    for (i, pair) in pairs.iter().enumerate() {
        let k = pair.cutoff;
        let expected = k.is_none_or(|k| pair.distance <= k).then_some(pair.distance);
        match task {
            Task::Score => {
                let (q, t) = (&pair.query, &pair.target);
                assert_eq!(case.mode.bitfront(q, t, k), expected, "Bitfront, pair {i}");
                assert_eq!(wfa2.score(q, t, k), expected, "WFA2, pair {i}");
                if k.is_some() {
                    assert_eq!(case.mode.rust_edlib(q, t, k), expected, "rust-edlib, pair {i}");
                }
                if let Some((_, score)) = previous {
                    assert_eq!(score(case.mode, q, t, k), expected, "previous Bitfront, pair {i}");
                }
            }
            Task::Path => path::validate(case.mode, pair, expected, &mut wfa2),
        }
    }

    let mode = case.mode.name();
    // Keep the original template-first ID even when swapping puts the edited,
    // possibly variable-length input in the target position.
    let target_len = case.length + case.flank;
    let mut id = format!("{mode}-{}x{target_len}-r{}-p{count}", case.length, case.rate_bp);
    if let Scenario::Seeded(Some(k)) = case.scenario {
        use std::fmt::Write as _;
        write!(id, "-k{k}").unwrap();
    } else if let Some(name) = case.scenario.name() {
        id.push('-');
        id.push_str(name);
    }
    if case.swap {
        id.push_str("-swap");
    }
    let (label, edlib_id) = match task {
        Task::Score => ("score", "rust-edlib"),
        // Keep the old unit-path measurements out of the run-path baseline.
        Task::Path => { id.push_str("-path"); ("path", "rust-edlib-runs") }
    };
    let directory = report.as_ref().map(|r| r.home.join(&id));
    let dir = directory.as_deref();
    let mut group = criterion.benchmark_group(id);
    group.throughput(Throughput::Elements(count as u64));
    match task {
        Task::Score => {
            measure(&mut group, dir, "bitfront", pairs, |q, t, k| case.mode.bitfront(q, t, k))?;
            measure(&mut group, dir, edlib_id, pairs, |q, t, k| case.mode.rust_edlib(q, t, k))?;
            measure(&mut group, dir, "wfa2", pairs, |q, t, k| wfa2.score(q, t, k))?;
            if let Some((rev, score)) = previous {
                measure(&mut group, dir, rev, pairs, |q, t, k| score(case.mode, q, t, k))?;
            }
        }
        Task::Path => {
            measure(&mut group, dir, "bitfront", pairs, |q, t, k| case.mode.bitfront_path(q, t, k))?;
            measure(&mut group, dir, edlib_id, pairs, |q, t, k| case.mode.edlib_path(q, t, k))?;
            measure(&mut group, dir, "wfa2", pairs, |q, t, k| wfa2.path(q, t, k))?;
        }
    }
    group.finish();

    if let Some((report, dir)) = report.zip(dir) {
        let ns = |engine| median_per_pair(dir, engine, count);
        let (b, e, w) = (ns("bitfront")?, ns(edlib_id)?, ns("wfa2")?);
        let out = &mut report.out;
        write!(out, "| {label} ")?;
        write_case(out, case, pairs)?;
        write!(out, " {b:.0} | {e:.0} | {w:.0} | {:.2}× | {:.2}× |", e / b, w / b)?;
        if let Some((rev, _)) = previous {
            let p = ns(rev)?;
            write!(out, " {p:.0} | {:.2}× |", p / b)?;
        }
        writeln!(out)?;
    }
    Ok(())
}

pub fn run() {
    type Score = fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>;
    run_panel(None::<(&str, Score)>).expect("alignment benchmark failed");
}

pub fn run_with_revision(
    revision: &str, score: impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>,
) {
    run_panel(Some((revision, score))).expect("revision benchmark failed");
}

fn run_panel(
    previous: Option<(&str, impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>)>,
) -> io::Result<()> {
    let summary = std::env::var_os("BITFRONT_BENCH_SUMMARY");
    let options = if let Some(path) = &summary {
        match Options::parse(std::env::args().skip(1)) {
            Ok(Some(options)) => options,
            Ok(None) => return writeln!(File::create(path)?, "{USAGE}"),
            Err(message) => {
                eprintln!("{message}\n{USAGE}");
                std::process::exit(2);
            }
        }
    } else {
        Options::default()
    };
    let tasks = options.tasks.unwrap_or(if previous.is_some() { &[Task::Score] } else { &[Task::Score, Task::Path] });
    if previous.is_some() && tasks.contains(&Task::Path) {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "--compare supports score rows only"));
    }
    let samples = options.samples();
    let (pairs, warm_up_ms, measurement_ms) = if options.smoke { (1, 5, 20) } else { (4, 100, 250) };
    let cases = options.case.map_or_else(|| workload::cases(options.smoke), |case| vec![case]);
    let previous = previous.as_ref().map(|(rev, score)| (*rev, score));
    if heap::in_worker() {
        return memory::report(&mut io::stdout().lock(), &cases, pairs, tasks, previous);
    }
    let mut report = if let Some(path) = summary {
        let home = std::env::var_os("CRITERION_HOME").expect("run.sh sets CRITERION_HOME");
        Some(Markdown { home: home.into(), out: BufWriter::new(File::create(path)?) })
    } else {
        None
    };
    let mut criterion = Criterion::default()
        .sample_size(samples)
        .warm_up_time(Duration::from_millis(warm_up_ms))
        .measurement_time(Duration::from_millis(measurement_ms))
        .nresamples(10_000);
    if let Some(report) = &mut report {
        criterion = criterion.without_plots();
        let out = &mut report.out;
        writeln!(out, "### Exact alignments: Bitfront vs. rust-edlib and WFA2-lib\n")?;
        writeln!(out, "{pairs} seeded pair(s) per row. Criterion median ns/pair over {samples} samples; \
            ratios above 1 favor Bitfront. WFA2: MemoryHigh, no heuristics. Path rows include owned runs \
            for Bitfront and rust-edlib, and packed CIGAR for WFA2; replay is outside the timer. \
            Global extra edits subtract the mandatory length gap from each pair's exact distance \
            before taking the range. A cutoff range gives the minimum and maximum per-pair cutoffs.\n")?;
        let (extra_columns, extra_alignment) = if let Some((rev, _)) = &previous {
            writeln!(out, "Bitfront is the current working tree; previous is `{rev}`.\n")?;
            (" previous ns/pair | previous / Bitfront |", "---:|---:|")
        } else {
            ("", "")
        };
        writeln!(out, "| Task |{CASE_COLUMNS} Bitfront ns/pair | rust-edlib ns/pair | WFA2 ns/pair \
            | rust-edlib / Bitfront | WFA2 / Bitfront |{extra_columns}")?;
        writeln!(out, "|---|{CASE_ALIGNMENT}---:|---:|---:|---:|---:|{extra_alignment}")?;
    } else {
        // Direct cargo bench uses Criterion's native filters and reporting.
        criterion = criterion.configure_from_args();
    }
    for case in cases {
        let inputs = case_pairs(&case, pairs);
        for &task in tasks {
            compare(&mut criterion, report.as_mut(), &case, &inputs, task, previous)?;
        }
    }
    criterion.final_summary();
    if let Some(report) = report {
        heap::append_tables(report.out.into_inner().map_err(|err| err.into_error())?)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{Case, Mode, Pair, Scenario, case_pairs, workload, write_case};

    #[test]
    fn held_out_cutoffs_bracket_identical_indel_pairs_and_stay_loose_for_sparse_hits() {
        let standard = workload::cases(false);
        let pairs = |scenario| {
            let case = standard.iter().find(|c| c.scenario == scenario).unwrap();
            case_pairs(case, workload::HELD_OUT.len())
        };
        assert!(workload::cases(true).iter().all(|c| matches!(c.scenario, Scenario::Seeded(_))));
        let hits = pairs(Scenario::MixedAt);
        let distances = || hits.iter().map(|p| p.distance);
        assert_eq!((distances().min(), distances().max()), (Some(955), Some(965)));
        for (i, (miss, hit)) in pairs(Scenario::MixedBelow).iter().zip(&hits).enumerate() {
            assert!(miss.query == hit.query && miss.target == hit.target, "mixed inputs differ for seed {i}");
            assert_eq!(miss.cutoff, Some(hit.distance - 1));
            assert_eq!(hit.cutoff, Some(hit.distance));
        }
        for sparse in pairs(Scenario::SparseLoose) {
            assert_eq!((sparse.distance, sparse.cutoff), (8, Some(4096)));
        }
    }

    #[test]
    fn case_row_reports_actual_lengths_and_extra_distance_per_pair() {
        let pairs = [
            Pair { query: b"AAAAA".to_vec(), target: b"A".to_vec(), distance: 4, cutoff: None },
            Pair { query: b"AAAAA".to_vec(), target: b"CCC".to_vec(), distance: 5, cutoff: None },
        ];
        let case = Case { swap: true, ..Case::seeded(Mode::Global, 3, 10_000, 2, None) };
        let mut row = Vec::new();
        write_case(&mut row, &case, &pairs).unwrap();
        assert_eq!(
            std::str::from_utf8(&row).unwrap(),
            "| global (swapped) | 3 | 5 | 1–3 | 100.00% | 4–5 | 0–3 |",
        );
    }
}
