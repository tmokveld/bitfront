use crate::{Mode, workload::Case};

pub(super) const USAGE: &str = "Usage: run.sh [--compare REV] [--samples N] [--smoke | --case MODE SOURCE_BP RATE_BP [CUTOFF]]\n  \
    MODE: global, prefix or suffix. RATE_BP: edits per 10,000 source bytes (0–10000).\n  \
    CUTOFF: inclusive maximum distance (0–2147483646).\n  \
    --samples N (or --sample-size N): Criterion samples per engine and workload; minimum 10, default 30 (10 for --smoke).";

#[derive(Default)]
pub(super) struct Options {
    pub smoke: bool,
    pub case: Option<Case>,
    samples: Option<usize>,
}

impl Options {
    pub fn samples(&self) -> usize {
        self.samples.unwrap_or(if self.smoke { 10 } else { 30 })
    }

    pub fn parse(args: impl Iterator<Item = String>) -> Result<Option<Self>, &'static str> {
        let mut args = args.peekable();
        let mut options = Self::default();
        while let Some(flag) = args.next() {
            match flag.as_str() {
                "--bench" => (), // Cargo supplies this to harness-free benchmarks.
                "-h" | "--help" => return Ok(None),
                "--samples" | "--sample-size" if options.samples.is_none() => {
                    let n = args.next().and_then(|s| s.parse::<usize>().ok()).filter(|&n| n >= 10)
                        .ok_or("sample count must be an integer of at least 10")?;
                    options.samples = Some(n);
                }
                "--smoke" if !options.smoke && options.case.is_none() => options.smoke = true,
                "--case" if !options.smoke && options.case.is_none() => {
                    let mode = match args.next().as_deref() {
                        Some("global") => Mode::Global,
                        Some("prefix") => Mode::Prefix,
                        Some("suffix") => Mode::Suffix,
                        _ => return Err("expected a search mode"),
                    };
                    let length = args.next().and_then(|s| s.parse().ok())
                        .filter(|&n| n > 0 && n <= i32::MAX as usize / 6).ok_or("invalid source length")?;
                    let rate_bp = args.next().and_then(|s| s.parse().ok())
                        .filter(|&r| r <= 10_000).ok_or("invalid edit rate")?;
                    let cutoff = args.next_if(|s| !s.starts_with('-')).map(|s| {
                        s.parse::<usize>().ok().filter(|&k| k < i32::MAX as usize).ok_or("invalid cutoff")
                    }).transpose()?;
                    let flank = if mode == Mode::Global { 0 } else { 3 * length };
                    options.case = Some(Case { mode, length, rate_bp, flank, cutoff });
                }
                _ => return Err("unexpected or repeated benchmark option"),
            }
        }
        Ok(Some(options))
    }
}
