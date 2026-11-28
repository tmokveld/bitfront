use crate::{Mode, Task, workload::{Case, REGRESSIONS, Scenario}};

pub(super) const USAGE: &str = "Usage: run.sh [--compare REV] [--task score|path|both] [--samples N] \
    [--smoke | --case NAME | --case MODE TEMPLATE_BP RATE_BP [CUTOFF] [--flank BP] [--swap]]\n  \
    --task: defaults to both; --compare supports only score and defaults to it.\n  \
    NAME: mixed-below, mixed-at or sparse-loose (fixed regression inputs).\n  \
    MODE: global, prefix or suffix. TEMPLATE_BP: positive unedited length in bases.\n  \
    RATE_BP: edits per 10,000 template bases (0–10000).\n  \
    CUTOFF: inclusive maximum distance (0–2147483646).\n  \
    --flank BP: append BP bases to the generated target (custom global cases only; default 0).\n  \
    --swap: exchange the generated query and target (custom global cases only).\n  \
    --flank and --swap can appear before or after --case.\n  \
    --samples N (or --sample-size N): Criterion samples per engine and workload; minimum 10, default 30 (10 for --smoke).";

// The native engines use signed 32-bit lengths, and Bitfront requires the sum
// of the input lengths to be less than i32::MAX.
const MAX_PAIR_BP: usize = i32::MAX as usize - 1;

#[derive(Default)]
pub(super) struct Options {
    pub smoke: bool,
    pub case: Option<Case>,
    pub tasks: Option<&'static [Task]>,
    samples: Option<usize>,
}

impl Options {
    pub fn samples(&self) -> usize {
        self.samples.unwrap_or(if self.smoke { 10 } else { 30 })
    }

    pub fn parse(args: impl Iterator<Item = String>) -> Result<Option<Self>, &'static str> {
        let mut args = args.peekable();
        let mut options = Self::default();
        let mut flank = None;
        let mut swap = false;
        while let Some(flag) = args.next() {
            match flag.as_str() {
                "--bench" => (), // Cargo supplies this to harness-free benchmarks.
                "-h" | "--help" => return Ok(None),
                "--task" if options.tasks.is_none() => {
                    options.tasks = Some(match args.next().as_deref() {
                        Some("score") => &[Task::Score],
                        Some("path") => &[Task::Path],
                        Some("both") => &[Task::Score, Task::Path],
                        _ => return Err("expected --task score, path, or both"),
                    });
                }
                "--samples" | "--sample-size" if options.samples.is_none() => {
                    let n = args.next().and_then(|s| s.parse::<usize>().ok()).filter(|&n| n >= 10)
                        .ok_or("sample count must be an integer of at least 10")?;
                    options.samples = Some(n);
                }
                "--smoke" if !options.smoke && options.case.is_none() => options.smoke = true,
                "--flank" if flank.is_none() => {
                    flank = Some(args.next().and_then(|s| s.parse::<usize>().ok())
                        .ok_or("invalid flank length")?);
                }
                "--swap" if !swap => swap = true,
                "--case" if !options.smoke && options.case.is_none() => {
                    let selector = args.next().ok_or("expected a search mode or named case")?;
                    if let Some(case) = REGRESSIONS.into_iter().find(|c| c.scenario.name() == Some(selector.as_str())) {
                        options.case = Some(case);
                        continue;
                    }
                    let mode = match selector.as_str() {
                        "global" => Mode::Global,
                        "prefix" => Mode::Prefix,
                        "suffix" => Mode::Suffix,
                        _ => return Err("unknown search mode or named case"),
                    };
                    let length = args.next().and_then(|s| s.parse().ok())
                        .filter(|&n| n > 0 && n <= MAX_PAIR_BP / 6).ok_or("invalid template length")?;
                    let rate_bp = args.next().and_then(|s| s.parse().ok())
                        .filter(|&r| r <= 10_000).ok_or("invalid edit rate")?;
                    let cutoff = args.next_if(|s| !s.starts_with('-')).map(|s| {
                        s.parse::<usize>().ok().filter(|&k| k < i32::MAX as usize).ok_or("invalid cutoff")
                    }).transpose()?;
                    let flank = if mode == Mode::Global { 0 } else { 3 * length };
                    options.case = Some(Case::seeded(mode, length, rate_bp, flank, cutoff));
                }
                _ => return Err("unexpected or repeated benchmark option"),
            }
        }
        if flank.is_some() || swap {
            let Some(case) = options.case.as_mut().filter(|c| c.mode == Mode::Global && matches!(c.scenario, Scenario::Seeded(_))) else {
                return Err("--flank and --swap require a custom global --case");
            };
            if let Some(flank) = flank {
                // The generated query can be up to twice the template length.
                if flank > MAX_PAIR_BP - 3 * case.length {
                    return Err("flank makes the input pair too long");
                }
                case.flank = flank;
            }
            case.swap = swap;
        }
        Ok(Some(options))
    }
}

#[cfg(test)]
mod tests {
    use super::Options;
    use crate::{Mode, workload};

    fn parse(args: &str) -> Result<Option<Options>, &'static str> {
        Options::parse(args.split_ascii_whitespace().map(str::to_owned))
    }

    #[test]
    fn custom_flanked_case_preserves_the_standard_pairs_and_can_swap_them() {
        let standard = workload::cases(false)
            .into_iter()
            .find(|c| c.mode == Mode::Global && c.length == 256 && c.rate_bp == 500 && c.flank == 768)
            .unwrap();
        let custom = parse("--case global 256 500 --flank 768").unwrap().unwrap().case.unwrap();
        let swapped = parse("--swap --flank 768 --case global 256 500").unwrap().unwrap().case.unwrap();
        for i in 0..4 {
            let original = standard.pair(i);
            assert_eq!(custom.pair(i), original);
            let (q, t) = original;
            assert_eq!(swapped.pair(i), (t, q));
        }
        let bounded = parse("--flank 768 --case global 256 1250 31 --swap").unwrap().unwrap().case.unwrap();
        assert_eq!(bounded.scenario.cutoff(32), Some(31));
        assert!(bounded.swap);
    }

    #[test]
    fn rejects_invalid_flank_and_swap_options() {
        for args in [
            "--flank 5",
            "--smoke --swap",
            "--case prefix 256 500 --flank 5",
            "--swap --case suffix 256 500",
            "--case global 0 500 --flank 5",
            "--case global 256 500 --flank -1",
            "--case global 256 500 --flank 2147483647",
            "--case global 256 500 --flank 5 --flank 5",
            "--swap --case global 256 500 --swap",
            "--case mixed-at --swap",
            "--flank 1 --case sparse-loose",
        ] {
            assert!(parse(args).is_err(), "{args}");
        }
    }
}
