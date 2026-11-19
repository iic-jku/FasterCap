#!/usr/bin/env python3
"""Run one trusted build/test at a time with conservative POSIX resource guards.

RSS, host memory and directory sizes are sampled, not hard kernel limits. Children
must remain in the process group (ordinary cmake/ctest/FasterCap commands do).
"""
import argparse
from contextlib import ExitStack
import fcntl
import json
import math
import os
from pathlib import Path
import re
import resource
import selectors
import shutil
import signal
import subprocess
import sys
import tempfile
import time

MIB = 1024 * 1024


def tree_bytes(root):
    total = 0
    for directory, _, files in os.walk(root, followlinks=False):
        for name in files:
            try:
                path = Path(directory, name)
                if not path.is_symlink():
                    total += path.stat().st_size
            except FileNotFoundError:
                pass
    return total


def group_rss(pgid):
    result = subprocess.run(['/bin/ps', '-axo', 'pgid=,rss='],
                            check=True, capture_output=True, text=True, timeout=2)
    return sum(int(row[1]) * 1024 for line in result.stdout.splitlines()
               if len(row := line.split()) == 2 and int(row[0]) == pgid)


def host_free_percent():
    if sys.platform == 'darwin':
        result = subprocess.run(['/usr/bin/memory_pressure', '-Q'], check=True,
                                capture_output=True, text=True, timeout=2)
        match = re.search(r'System-wide memory free percentage:\s*(\d+)%', result.stdout)
        if not match:
            raise RuntimeError('cannot read host memory headroom')
        return float(match[1])
    values = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        key, value = line.split(':', 1)
        values[key] = int(value.split()[0])
    return 100.0 * values['MemAvailable'] / values['MemTotal']


def terminate_group(process):
    # Kill any descendants still holding a pipe even if the group leader exited.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def positive(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be finite and positive')
    return number


def main():
    def interrupted(signum, frame):
        raise KeyboardInterrupt('received signal ' + str(signum))

    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGHUP, interrupted)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True,
                        help='private scratch directory shared by all guarded runs')
    parser.add_argument('--label', default='run')
    parser.add_argument('--seconds', type=positive, default=120)
    parser.add_argument('--rss-mib', type=positive, default=1024)
    parser.add_argument('--reserve-mib', type=positive, default=5120)
    parser.add_argument('--scratch-mib', type=positive, default=768)
    parser.add_argument('--temp-mib', type=positive, default=256)
    parser.add_argument('--log-mib', type=positive, default=8)
    parser.add_argument('--host-free-percent', type=positive, default=20)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.command[:1] == ['--']:
        args.command.pop(0)
    if not args.command or not re.fullmatch(r'[A-Za-z0-9_.-]+', args.label):
        parser.error('command and a simple log label are required')
    root = args.root.resolve()
    root.mkdir(parents=True, exist_ok=True, mode=0o700)
    with (root / 'run.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            print('guard: another run owns this scratch directory', file=sys.stderr)
            return 125
        return run(args, root)


def run(args, root):
    process = None
    peak = 0
    written = 0
    reason = None
    started = time.monotonic()
    log_path = root / (args.label + '-' + str(time.time_ns()) + '.log')
    try:
        # Fail closed if ps or the host memory probe is unavailable.
        group_rss(os.getpgrp())
        if host_free_percent() < args.host_free_percent:
            raise RuntimeError('host memory headroom below reserve')
        if shutil.disk_usage(root).free < args.reserve_mib * MIB:
            raise RuntimeError('free disk below reserve')
        if tree_bytes(root) > args.scratch_mib * MIB:
            raise RuntimeError('scratch directory already exceeds budget')
        with tempfile.TemporaryDirectory(prefix='run-', dir=root) as temporary, \
                log_path.open('xb') as log, selectors.DefaultSelector() as selector, \
                ExitStack() as children:
            env = os.environ.copy()
            env.update(TMPDIR=temporary, TMP=temporary, TEMP=temporary,
                       OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='1',
                       VECLIB_MAXIMUM_THREADS='1', MKL_NUM_THREADS='1',
                       CMAKE_BUILD_PARALLEL_LEVEL='1', CTEST_PARALLEL_LEVEL='1',
                       PYTHONDONTWRITEBYTECODE='1')

            def limits():
                resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
                # Inherited per-file bound also limits solver files outside TMPDIR.
                resource.setrlimit(resource.RLIMIT_FSIZE, (128 * MIB, 128 * MIB))

            process = subprocess.Popen(args.command, env=env, start_new_session=True,
                                       preexec_fn=limits, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT)
            children.callback(process.stdout.close)
            children.callback(terminate_group, process)
            selector.register(process.stdout, selectors.EVENT_READ)
            next_memory = next_disk = 0.0
            while selector.get_map() or process.poll() is None:
                now = time.monotonic()
                if now - started > args.seconds:
                    raise RuntimeError('elapsed time limit')
                if now >= next_memory:
                    peak = max(peak, group_rss(process.pid))
                    if peak > args.rss_mib * MIB:
                        raise RuntimeError('process-group RSS limit')
                    next_memory = now + 0.1
                if now >= next_disk:
                    if host_free_percent() < args.host_free_percent:
                        raise RuntimeError('host memory headroom below reserve')
                    if shutil.disk_usage(root).free < args.reserve_mib * MIB:
                        raise RuntimeError('free disk below reserve')
                    if tree_bytes(temporary) > args.temp_mib * MIB:
                        raise RuntimeError('temporary-file limit')
                    if tree_bytes(root) > args.scratch_mib * MIB:
                        raise RuntimeError('scratch directory limit')
                    next_disk = now + 1.0
                for key, _ in selector.select(timeout=0.1):
                    chunk = os.read(key.fileobj.fileno(), 65536)
                    if not chunk:
                        selector.unregister(key.fileobj)
                        continue
                    room = max(0, int(args.log_mib * MIB) - written)
                    log.write(chunk[:room])
                    written += min(room, len(chunk))
                    if len(chunk) > room:
                        raise RuntimeError('log output limit')
            code = process.wait()
    except (OSError, RuntimeError, subprocess.SubprocessError, KeyboardInterrupt) as error:
        reason = str(error) or 'interrupted'
        code = 125
    elapsed = time.monotonic() - started
    report = {'exit': code, 'guard': reason, 'seconds': round(elapsed, 2),
              'peak_rss_mib': round(peak / MIB, 2), 'log': str(log_path),
              'command': args.command, 'cwd': os.getcwd()}
    with log_path.with_suffix('.json').open('x') as summary:
        summary.write(json.dumps(report) + '\n')
    print(json.dumps(report))
    if log_path.exists():
        with log_path.open('rb') as log:
            log.seek(max(0, log_path.stat().st_size - (2000 if code == 0 else 12000)))
            print(log.read().decode(errors='replace'), end='')
    return code if code >= 0 else 128 - code


if __name__ == '__main__':
    sys.exit(main())
