import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


BENCH = Path(__file__).resolve().parent
ROOT = BENCH.parents[1]
USAGE = "Usage: run.sh --compare REV [--samples N] [--smoke | --case MODE SOURCE_BP RATE_BP [CUTOFF]]"


def main(args):
    if not args:
        print(USAGE, file=sys.stderr)
        return 2
    if args[0] in ("--help", "-h"):
        print(USAGE)
        return 0
    ref = subprocess.run(
        ["git", "-C", ROOT, "rev-parse", "--verify", "--end-of-options", args[0] + "^{commit}"],
        stdout=subprocess.PIPE, text=True,
    )
    if ref.returncode:
        return ref.returncode
    commit = ref.stdout.strip()

    with tempfile.TemporaryDirectory(prefix="bitfront-compare-") as temp:
        project = Path(temp)
        # Distinct Cargo sources let the working tree and the Git revision keep
        # the same package name and version. Only the temporary lockfile changes.
        manifest = project / "Cargo.toml"
        manifest.write_text(f'''[package]
name = "bitfront-revision-bench"
version = "0.0.0"
edition = "2024"
publish = false

[workspace]

[[bench]]
name = "bitfront-revision-bench"
path = {json.dumps(str(BENCH / "revision.rs"), ensure_ascii=False)}
harness = false

[dependencies]
bitfront-bench = {{ path = {json.dumps(str(BENCH), ensure_ascii=False)} }}
previous = {{ package = "bitfront", git = "{ROOT.as_uri()}", rev = "{commit}" }}

[profile.bench]
lto = "thin"
codegen-units = 1
''', encoding="utf-8")
        shutil.copyfile(BENCH / "Cargo.lock", project / "Cargo.lock")
        env = dict(os.environ, BITFRONT_PREVIOUS_COMMIT=commit)
        env.setdefault("CARGO_TARGET_DIR", str(BENCH / "target"))
        return subprocess.run(
            ["cargo", "bench", "--manifest-path", manifest, "--bench", "bitfront-revision-bench", "--", *args[1:]], env=env,
        ).returncode


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
