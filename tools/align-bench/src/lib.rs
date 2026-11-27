mod native;
mod options;
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
use workload::Case;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Global,
    Prefix,
    Suffix,
}

impl Mode {
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

fn compare(
    criterion: &mut Criterion, report: Option<&mut Markdown>, case: Case, count: usize,
    previous: Option<(&str, &impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>)>,
) -> io::Result<()> {
    let pairs: Vec<Pair> = (0..count).map(|i| case.pair(i)).collect();
    let mut wfa2 = Wfa2::new(case.mode, case.cutoff);
    let (mut low, mut high, mut hits) = (usize::MAX, 0, 0);
    for (i, (q, t)) in pairs.iter().enumerate() {
        let e = case.mode.rust_edlib(q, t, case.cutoff);
        let distance = if case.cutoff.is_some() {
            case.mode.rust_edlib(q, t, None)
        } else {
            e
        }.expect("unbounded rust-edlib returned no distance");
        let expected = case.cutoff.is_none_or(|k| distance <= k).then_some(distance);
        assert_eq!(
            [case.mode.bitfront(q, t, case.cutoff), e, wfa2.score(q, t)],
            [expected; 3],
            "[Bitfront, rust-edlib, WFA2] disagree with rust-edlib's unbounded distance {distance}: \
             {:?}, cutoff {:?}, source {} bp, edit rate {}/10000, pair {i}",
            case.mode, case.cutoff, case.length, case.rate_bp,
        );
        if let Some((_, score)) = previous {
            assert_eq!(score(case.mode, q, t, case.cutoff), expected,
                "previous Bitfront disagrees at {:?}, cutoff {:?}, source {} bp, edit rate {}/10000, pair {i}",
                case.mode, case.cutoff, case.length, case.rate_bp);
        }
        low = low.min(distance);
        high = high.max(distance);
        hits += usize::from(expected.is_some());
    }

    let mode = match case.mode {
        Mode::Global => "global",
        Mode::Prefix => "prefix",
        Mode::Suffix => "suffix",
    };
    let target_len = pairs[0].1.len();
    let mut id = format!("{mode}-{}x{target_len}-r{}-p{count}", case.length, case.rate_bp);
    if let Some(k) = case.cutoff {
        use std::fmt::Write as _;
        write!(id, "-k{k}").unwrap();
    }
    let directory = report.as_ref().map(|r| r.home.join(&id));
    let dir = directory.as_deref();
    let mut group = criterion.benchmark_group(id);
    group.throughput(Throughput::Elements(count as u64));
    measure(&mut group, dir, "bitfront", &pairs, |q, t| case.mode.bitfront(q, t, case.cutoff))?;
    measure(&mut group, dir, "rust-edlib", &pairs, |q, t| case.mode.rust_edlib(q, t, case.cutoff))?;
    measure(&mut group, dir, "wfa2", &pairs, |q, t| wfa2.score(q, t))?;
    if let Some((rev, score)) = previous {
        measure(&mut group, dir, rev, &pairs, |q, t| score(case.mode, q, t, case.cutoff))?;
    }
    group.finish();

    if let Some((report, dir)) = report.zip(dir) {
        let ns = |engine| median_per_pair(dir, engine, count);
        let (b, e, w) = (ns("bitfront")?, ns("rust-edlib")?, ns("wfa2")?);
        let out = &mut report.out;
        write!(out, "| {mode}")?;
        if let Some(k) = case.cutoff {
            write!(out, " ≤ {k} ({hits}/{count} hit)")?;
        }
        write!(out, " | {} | {target_len} | {:.2}% | {low}", case.length, f64::from(case.rate_bp) / 100.)?;
        if low != high {
            write!(out, "–{high}")?;
        }
        write!(out, " | {b:.0} | {e:.0} | {w:.0} | {:.2}× | {:.2}× |", e / b, w / b)?;
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
    let mut report = None;
    let options = if let Some(path) = std::env::var_os("BITFRONT_BENCH_SUMMARY") {
        let mut out = BufWriter::new(File::create(path)?);
        let options = match Options::parse(std::env::args().skip(1)) {
            Ok(Some(options)) => options,
            Ok(None) => {
                writeln!(out, "{USAGE}")?;
                return out.flush();
            }
            Err(message) => {
                eprintln!("{message}\n{USAGE}");
                std::process::exit(2);
            }
        };
        let home = std::env::var_os("CRITERION_HOME").expect("run.sh sets CRITERION_HOME");
        report = Some(Markdown { home: home.into(), out });
        options
    } else {
        Options::default()
    };
    let samples = options.samples();
    let (pairs, warm_up_ms, measurement_ms) = if options.smoke { (1, 5, 20) } else { (4, 100, 250) };
    let mut criterion = Criterion::default()
        .sample_size(samples)
        .warm_up_time(Duration::from_millis(warm_up_ms))
        .measurement_time(Duration::from_millis(measurement_ms))
        .nresamples(10_000);
    if let Some(report) = &mut report {
        criterion = criterion.without_plots();
        let out = &mut report.out;
        writeln!(out, "### Exact edit distance: Bitfront vs. rust-edlib and WFA2-lib\n")?;
        writeln!(out, "{pairs} seeded pair(s) per row. Criterion median ns/pair over {samples} samples; \
            ratios above 1 favor Bitfront. WFA2: score only, MemoryHigh, no heuristics.\n")?;
        let (extra_columns, extra_alignment) = if let Some((rev, _)) = &previous {
            writeln!(out, "Bitfront is the current working tree; previous is `{rev}`.\n")?;
            (" previous ns/pair | previous / Bitfront |", "---:|---:|")
        } else {
            ("", "")
        };
        writeln!(out, "| Search | Source bp | Target bp | Nominal edit rate | Exact distance | Bitfront ns/pair \
            | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |{extra_columns}")?;
        writeln!(out, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|{extra_alignment}")?;
    } else {
        // Direct cargo bench uses Criterion's native filters and reporting.
        criterion = criterion.configure_from_args();
    }
    let cases = match options.case {
        Some(case) => vec![case],
        None => workload::cases(options.smoke),
    };
    for case in cases {
        compare(&mut criterion, report.as_mut(), case, pairs, previous.as_ref().map(|(r, f)| (*r, f)))?;
    }
    criterion.final_summary();
    if let Some(mut report) = report {
        report.out.flush()?;
    }
    Ok(())
}
