use crate::{CASE_ALIGNMENT, CASE_COLUMNS, Case, Mode, Pair, Task, Wfa2, heap::Probe, write_case};
use std::hint::black_box;
use std::io::{self, Write};

pub(super) fn report(
    out: &mut impl Write, cases: &[Case], count: usize, tasks: &[Task],
    previous: Option<(&str, &impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>)>,
) -> io::Result<()> {
    let probe = Probe::load()?;
    // Group by task in the report, but generate the input and reference distances only once per case.
    let mut tables = Vec::with_capacity(tasks.len());
    for &task in tasks {
        let mut table = Vec::new();
        header(&mut table, task, previous.map(|(rev, _)| rev))?;
        tables.push((task, table));
    }
    for case in cases {
        let pairs = crate::case_pairs(case, count);
        for (task, table) in &mut tables {
            row(table, &probe, case, &pairs, *task, previous)?;
        }
    }
    for (_, table) in tables {
        out.write_all(&table)?;
    }
    out.flush()
}

fn header(out: &mut impl Write, task: Task, previous: Option<&str>) -> io::Result<()> {
    let (title, extra_columns, extra_alignment) = match task {
        Task::Score => {
            if previous.is_some() {
                ("Score", " previous peak B | previous / Bitfront |", "---:|---:|")
            } else {
                ("Score", "", "")
            }
        }
        Task::Path => (
            "Path",
            " Bitfront returned B | rust-edlib returned B | WFA2 returned B |",
            "---:|---:|---:|",
        ),
    };
    writeln!(out, "\n### {title} memory\n")?;
    writeln!(out, "Peak live heap in allocator-usable bytes, maximum across the seeded pairs. \
        Includes owned output and WFA2's reusable native aligner. Excludes inputs, stack and allocator metadata; not RSS.")?;
    writeln!(out, "Ratios are engine peak / Bitfront peak: above 1× means Bitfront used less heap; \
        — means Bitfront used no heap.")?;
    if task == Task::Path {
        writeln!(out, "Returned bytes are the maximum heap freed by dropping one result, already included in peak.")?;
    } else if let Some(revision) = previous {
        writeln!(out, "Previous Bitfront: `{revision}`.")?;
    }
    writeln!(out)?;
    writeln!(out, "|{CASE_COLUMNS} Bitfront peak B | rust-edlib peak B | WFA2 peak B \
        | rust-edlib / Bitfront | WFA2 / Bitfront |{extra_columns}")?;
    writeln!(out, "|{CASE_ALIGNMENT}---:|---:|---:|---:|---:|{extra_alignment}")
}

fn write_peak_ratio(out: &mut impl Write, engine: usize, bitfront: usize) -> io::Result<()> {
    if bitfront == 0 {
        write!(out, " — |")
    } else {
        write!(out, " {:.2}× |", engine as f64 / bitfront as f64)
    }
}

struct Heap {
    peak: usize,
    returned: usize,
}

fn measure<S, T>(
    probe: &Probe, pairs: &[Pair], task: Task, setup: impl FnOnce() -> S,
    mut align: impl FnMut(&mut S, &[u8], &[u8], Option<usize>) -> T,
) -> Heap {
    probe.begin();
    let mut state = setup();
    let mut returned = 0;
    for pair in pairs {
        let result = black_box(align(
            &mut state,
            black_box(pair.query.as_slice()),
            black_box(pair.target.as_slice()),
            pair.cutoff,
        ));
        if task == Task::Path {
            let before = probe.live();
            drop(result);
            returned = returned.max(before - probe.live());
        } else {
            drop(result);
        }
    }
    drop(state);
    Heap { peak: probe.end(), returned }
}

fn row(
    out: &mut impl Write, probe: &Probe, case: &Case, pairs: &[Pair], task: Task,
    previous: Option<(&str, &impl Fn(Mode, &[u8], &[u8], Option<usize>) -> Option<usize>)>,
) -> io::Result<()> {
    let mode = case.mode;
    let new_wfa2 = || {
        let aligner = Wfa2::new(mode, pairs[0].cutoff, task);
        assert!(probe.live() > 0, "heap sampler did not observe the native WFA2 aligner");
        aligner
    };
    let (bitfront, edlib, wfa2) = match task {
        Task::Score => (
            measure(probe, pairs, task, || (), |_, q, t, k| mode.bitfront(q, t, k)),
            measure(probe, pairs, task, || (), |_, q, t, k| mode.rust_edlib(q, t, k)),
            measure(probe, pairs, task, new_wfa2, Wfa2::score),
        ),
        Task::Path => (
            measure(probe, pairs, task, || (), |_, q, t, k| mode.bitfront_path(q, t, k)),
            measure(probe, pairs, task, || (), |_, q, t, k| mode.edlib_path(q, t, k)),
            measure(probe, pairs, task, new_wfa2, Wfa2::path),
        ),
    };
    let previous = if task == Task::Score {
        previous.map(|(_, score)| measure(probe, pairs, task, || (), |_, q, t, k| score(mode, q, t, k)).peak)
    } else {
        None
    };

    write_case(out, case, pairs)?;
    write!(out, " {} | {} | {} |", bitfront.peak, edlib.peak, wfa2.peak)?;
    write_peak_ratio(out, edlib.peak, bitfront.peak)?;
    write_peak_ratio(out, wfa2.peak, bitfront.peak)?;
    if task == Task::Path {
        write!(out, " {} | {} | {} |", bitfront.returned, edlib.returned, wfa2.returned)?;
    }
    if let Some(peak) = previous {
        write!(out, " {peak} |")?;
        write_peak_ratio(out, peak, bitfront.peak)?;
    }
    writeln!(out)
}
