mod native;
mod workload;

use native::Wfa2;
use std::hint::black_box;
use std::time::{Duration, Instant};
use workload::Case;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Global,
    AtMost(usize),
    Prefix,
    Suffix,
}

impl Mode {
    fn bitfront(self, query: &[u8], target: &[u8]) -> Option<usize> {
        match self {
            Self::Global => bitfront::distance(query, target).map(Some),
            Self::AtMost(k) => bitfront::distance_at_most(query, target, k),
            Self::Prefix => bitfront::prefix_distance(query, target).map(Some),
            Self::Suffix => bitfront::suffix_distance(query, target).map(Some),
        }
        .expect("bitfront failed to score benchmark input")
    }
}

type Pair = (Vec<u8>, Vec<u8>);
const SAMPLES: usize = 5;

fn batch(
    pairs: &[Pair],
    passes: usize,
    mut score: impl FnMut(&[u8], &[u8]) -> Option<usize>,
) -> Duration {
    let start = Instant::now();
    for _ in 0..passes {
        for (query, target) in pairs {
            let _ = black_box(score(black_box(query), black_box(target)));
        }
    }
    start.elapsed()
}

fn medians(
    pairs: usize,
    first: usize,
    min_time: Duration,
    mut run: impl FnMut(usize, usize) -> Duration,
) -> [f64; 3] {
    // Calibrate whole-panel passes per engine: even the slowest sample sees every pair.
    let mut passes = [1; 3];
    for (engine, n) in passes.iter_mut().enumerate() {
        while run(engine, *n) < min_time {
            *n *= 2;
        }
    }
    let mut ns = [[0.; SAMPLES]; 3];
    for sample in 0..SAMPLES {
        for offset in 0..3 {
            let engine = (first + sample + offset) % 3;
            let elapsed = run(engine, passes[engine]);
            ns[engine][sample] = elapsed.as_secs_f64() * 1e9 / (passes[engine] * pairs) as f64;
        }
    }
    ns.map(|mut times| {
        times.sort_unstable_by(f64::total_cmp);
        times[SAMPLES / 2]
    })
}

fn compare(case: Case, ordinal: usize, count: usize, min_time: Duration) {
    let pairs: Vec<Pair> = (0..count).map(|i| case.pair(i)).collect();
    let mut wfa2 = Wfa2::new(case.mode);
    let (mut low, mut high, mut hits) = (usize::MAX, 0, 0);
    for (i, (q, t)) in pairs.iter().enumerate() {
        let e = case.mode.rust_edlib(q, t);
        let distance = match case.mode {
            Mode::AtMost(_) => Mode::Global.rust_edlib(q, t),
            _ => e,
        }.expect("unbounded rust-edlib returned no distance");
        let expected = match case.mode {
            Mode::AtMost(k) if distance > k => None,
            _ => Some(distance),
        };
        assert_eq!(
            [case.mode.bitfront(q, t), e, wfa2.score(q, t)],
            [expected; 3],
            "[Bitfront, rust-edlib, WFA2] disagree with rust-edlib's unbounded distance {distance}: \
             {:?}, source {} bp, edit rate {}/10000, pair {i}",
            case.mode, case.length, case.rate_bp,
        );
        low = low.min(distance);
        high = high.max(distance);
        hits += usize::from(expected.is_some());
    }

    let [b, e, w] = medians(count, ordinal, min_time, |engine, passes| match engine {
        0 => batch(&pairs, passes, |q, t| case.mode.bitfront(q, t)),
        1 => batch(&pairs, passes, |q, t| case.mode.rust_edlib(q, t)),
        2 => batch(&pairs, passes, |q, t| wfa2.score(q, t)),
        _ => unreachable!(),
    });
    match case.mode {
        Mode::Global => print!("| global"),
        Mode::AtMost(k) => print!("| global ≤ {k} ({hits}/{count} hit)"),
        Mode::Prefix => print!("| prefix"),
        Mode::Suffix => print!("| suffix"),
    }
    print!(
        " | {} | {} | {:.2}% | {low}",
        case.length,
        pairs[0].1.len(),
        f64::from(case.rate_bp) / 100.
    );
    if low != high {
        print!("–{high}");
    }
    println!(
        " | {b:.0} | {e:.0} | {w:.0} | {:.2}× | {:.2}× |",
        e / b,
        w / b,
    );
}


fn main() {
    const USAGE: &str = "Usage: bitfront-bench [--smoke | --case MODE SOURCE_BP RATE_BP]\n  \
        MODE: global, prefix or suffix. RATE_BP: edits per 10,000 source bytes (0–10000).";
    let invalid = || -> ! {
        eprintln!("{USAGE}");
        std::process::exit(2);
    };
    let mut args = std::env::args().skip(1);
    let flag = args.next();
    let smoke = flag.as_deref() == Some("--smoke");
    let cases = match flag.as_deref() {
        None | Some("--smoke") => workload::cases(smoke),
        Some("--help" | "-h") if args.next().is_none() => {
            println!("{USAGE}");
            return;
        }
        Some("--case") => {
            let mode = match args.next().as_deref() {
                Some("global") => Mode::Global,
                Some("prefix") => Mode::Prefix,
                Some("suffix") => Mode::Suffix,
                _ => invalid(),
            };
            let length = args.next().and_then(|s| s.parse().ok())
                .filter(|&n| n > 0 && n <= i32::MAX as usize / 6);
            let rate = args.next().and_then(|s| s.parse().ok()).filter(|&r| r <= 10_000);
            let (Some(length), Some(rate_bp)) = (length, rate) else { invalid() };
            let flank = if mode == Mode::Global { 0 } else { 3 * length };
            vec![Case { mode, length, rate_bp, flank }]
        }
        _ => invalid(),
    };
    if args.next().is_some() {
        invalid();
    }
    let (pairs, min_time) = if smoke {
        (1, Duration::from_micros(100))
    } else {
        (4, Duration::from_millis(2))
    };
    println!("### Exact edit distance: Bitfront vs. rust-edlib and WFA2-lib\n");
    println!(
        "{pairs} seeded pair(s) per row. Median ns/pair over {SAMPLES} samples; ratios above 1 \
         favor Bitfront. WFA2: score only, MemoryHigh, no heuristics.\n"
    );
    println!(
        "| Search | Source bp | Target bp | Nominal edit rate | Exact distance | Bitfront ns/pair \
         | rust-edlib ns/pair | WFA2 ns/pair | rust-edlib / Bitfront | WFA2 / Bitfront |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|");
    for (i, case) in cases.into_iter().enumerate() {
        compare(case, i, pairs, min_time);
    }
}
