"""Check opt-in timing diagnostics and bounded benchmark preflight."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import geometries
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
import benchmark_solver as benchmark


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", required=True)
    parser.add_argument("--counter", required=True)
    args = parser.parse_args()
    executable, counter = str(Path(args.fastercap).resolve()), str(Path(args.counter).resolve())
    env = dict(os.environ, OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
               OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(prefix="fc-timing-") as temporary:
        root = Path(temporary)
        model = geometries.GEOMETRIES["three_cubes"].write(str(root))
        for options in (["-m0.1", "-d0.1", "-t1e-9", "-pj"], ["-a0.01"]):
            reference = None
            for detailed in (False, True):
                report = root / "report.json"
                output, _ = benchmark.invoke([executable, "-b", model, *options, "-f0",
                    "--report=" + str(report), *(["-i"] if detailed else [])], root, env)
                records = json.loads(report.read_text())["iterations"]
                report.unlink()
                result = [(r["matrix"], r["panels"], r["links"], r["solve_status"]) for r in records]
                if reference is None:
                    reference = result
                    assert "Hierarchical products:" not in output
                else:
                    assert result == reference, "-i changes full-precision results or work"
                    sections = output.split("Solve statistics:")[1:]
                    assert len(sections) == len(records)
                    for section in sections:
                        times = benchmark.timings(section)
                        assert times["link_io"] == 0, "in-core multiplication recorded link loading"
                        assert times["charges"] + times["fields"] + times["links_including_io"] > 0
        # Exercise the whole harness, including cleanup and same-work checks.
        destination = root / "benchmark.json"
        subprocess.run([sys.executable, str(Path(benchmark.__file__)), "--fastercap", executable,
                        "--counter", counter, "--workdir", str(root), "--output", str(destination),
                        "--cases", "high", "--threads", "1", "2", "--repeats", "1"],
                       check=True, env=env, timeout=20)
        data = json.loads(destination.read_text())
        assert len(data["records"]) == 2 and len(data["summaries"]) == 2
        assert not list(root.glob("benchmark-*")), "generated benchmark geometry was retained"
        # Tight depth-3 grouping exceeds the cap; fail before link allocation.
        large = root / "count-only"
        large.mkdir()
        model = geometries.coated_spheres_geometry(100, depth=3).write(str(large))
        process = subprocess.run([counter, model, "1e-36"], capture_output=True, text=True,
                                 timeout=15, env=env)
        assert process.returncode == 3 and "BENCHMARK panels 3840" in process.stdout
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
