# Bounded local builds and tests

Run trusted POSIX build/test commands through `tools/run_guarded.py --root <scratch> -- <command>`.
Use the same private scratch root for all runs; a nonblocking lock prevents overlap.
Build with `cmake --build <build> --parallel 1` and test with `ctest --test-dir <build> -j1`.

Defaults stop the command group at 1 GiB observed RSS, 120 seconds, less than 20% host
memory headroom, less than 5 GiB free disk, 768 MiB of scratch files, 256 MiB of private
temporary files, or 8 MiB of captured output. Cores are disabled and each child-written
file is capped at 128 MiB. The guard sets two OpenMP workers and one BLAS/build/test worker.
It reports elapsed time, peak sampled RSS, exit status and a unique bounded log.

Memory and directory sizes are sampled: they are not a hard kernel memory guarantee.
Commands must keep children in their process group. Do not use this to sandbox hostile
programs or commands that daemonize. An unavailable host/process probe prevents launch
or stops the run. Only the guard's unique temporary directory is cleaned; source files,
build outputs and logs are preserved. Keep the build beneath the scratch root and check
available space before creating additional source worktrees.

Run the small guard regressions with `python3 tests/test_resource_guard.py`. These include
intentional timeout, memory, output, disk-reserve and temporary-file-limit violations.
