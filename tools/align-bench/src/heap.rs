use std::{env, ffi::OsString, fs::File, io::{self, Write}, process::Command};

const WORKER: &str = "BITFRONT_BENCH_HEAP_WORKER";

pub(super) fn in_worker() -> bool {
    env::var_os(WORKER).is_some()
}

/// Re-exec the same benchmark (including --compare) with allocation hooks only in the child.
/// Its stdout shares the flushed timing report's file offset; Criterion never loads the probe.
pub(super) fn append_tables(mut output: File) -> io::Result<()> {
    let loader = if cfg!(target_os = "macos") {
        "DYLD_INSERT_LIBRARIES"
    } else if cfg!(all(target_os = "linux", target_env = "gnu")) {
        "LD_PRELOAD"
    } else {
        return writeln!(output, "\nHeap memory tables require macOS or Linux with glibc.");
    };
    // A relative loader path works even when Cargo's OUT_DIR contains spaces (LD_PRELOAD
    // splits on spaces as well as colons). The worker uses generated inputs and inherited stdout.
    let mut libraries = OsString::from("./bitfront-heap.so");
    if let Some(existing) = env::var_os(loader).filter(|s| !s.is_empty()) {
        libraries.push(":");
        libraries.push(existing);
    }
    let status = Command::new(env::current_exe()?).args(env::args_os().skip(1))
        .current_dir(env!("OUT_DIR")).env(loader, libraries).env(WORKER, "1")
        .stdout(output).status()?;
    if status.success() { Ok(()) } else { Err(io::Error::other(format!("heap measurement failed: {status}"))) }
}

pub(super) struct Probe {
    begin: extern "C" fn(),
    live: extern "C" fn() -> usize,
    end: extern "C" fn() -> usize,
}

impl Probe {
    #[cfg(unix)]
    pub(super) fn load() -> io::Result<Self> {
        let process = libloading::os::unix::Library::this();
        // SAFETY: these signatures are the heap.c ABI. The preloaded image remains resident
        // for the life of the process, independently of this symbol lookup handle.
        let probe = unsafe { Self {
            begin: *process.get(b"bf_heap_begin\0").map_err(io::Error::other)?,
            live: *process.get(b"bf_heap_live\0").map_err(io::Error::other)?,
            end: *process.get(b"bf_heap_end\0").map_err(io::Error::other)?,
        } };
        // An ignored preload or a non-system Rust allocator must not produce convincing zeros.
        probe.begin();
        let allocation = std::hint::black_box(Box::new([0u8; 64]));
        let live = probe.live();
        drop(allocation);
        if probe.end() < 64 || live < 64 {
            return Err(io::Error::other("heap sampler did not observe a Rust allocation"));
        }
        Ok(probe)
    }

    #[cfg(not(unix))]
    pub(super) fn load() -> io::Result<Self> {
        Err(io::Error::other("heap measurement requires macOS or Linux with glibc"))
    }

    pub(super) fn begin(&self) { (self.begin)(); }
    pub(super) fn live(&self) -> usize { (self.live)() }
    pub(super) fn end(&self) -> usize { (self.end)() }
}
