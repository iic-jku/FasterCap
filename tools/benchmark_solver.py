"""Bounded fixed-mesh benchmarks; invoke through tools/run_guarded.py on POSIX."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import sys
import tempfile
import time

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tests"))
import geometries


PHASES = {
    "charges": "Hierarchy charge/moments wall time",
    "links_including_io": "Hierarchy links wall time (includes I/O)",
    "fields": "Hierarchy downward field wall time",
    "link_io": "Hierarchy link loading wall time",
    "coefficient_generation": "Time for building potential matrix",
    "solve": "Time for gmres solving",
}


def timings(output):
    result = {}
    for name, label in PHASES.items():
        values = re.findall(re.escape(label) + r": ([^\s]+)s", output)
        if len(values) != 1:
            raise ValueError("expected one timing for " + label)
        value = float(values[0])
        if not math.isfinite(value) or value < 0:
            raise ValueError("invalid timing for " + label)
        result[name] = value
    if result["link_io"] > result["links_including_io"] + 1e-8:
        raise ValueError("I/O exceeds inclusive link time")
    result["link_compute"] = max(0.0, result["links_including_io"] - result["link_io"])
    return result


def fixture(name, depth):
    if name == "conductor":
        return geometries.layered_spheres_geometry([1.0], [], depth=depth)
    if name == "uneven":
        return geometries.coated_spheres_geometry(100, depth=2, split_coat=True, split_half=True)
    er = {"moderate": 10, "high": 100, "lossy": 10-2j}[name]
    return geometries.coated_spheres_geometry(er, depth=depth)


def invoke(command, directory, env, timeout=30):
    start = time.perf_counter()
    process = subprocess.run(command, cwd=directory, env=env, capture_output=True,
                             text=True, timeout=timeout)
    elapsed = time.perf_counter() - start
    output = process.stdout + process.stderr
    if process.returncode:
        raise RuntimeError("command failed (%d): %s\n%s" %
                           (process.returncode, command, output[-4000:]))
    return output, elapsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", type=Path, required=True)
    parser.add_argument("--counter", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--size", choices=("small", "medium", "large"), default="small")
    parser.add_argument("--cases", nargs="+", choices=("conductor", "moderate", "high", "uneven", "lossy"),
                        default=["conductor", "moderate", "high", "uneven", "lossy"])
    parser.add_argument("--threads", nargs="+", type=int, choices=(1, 2, 4), default=[1, 2, 4])
    parser.add_argument("--repeats", type=int, choices=range(1, 6), default=3)
    args = parser.parse_args()
    executable, counter = str(args.fastercap.resolve()), str(args.counter.resolve())
    args.workdir.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise ValueError("output already exists: " + str(args.output))
    if 1 not in args.threads:
        raise ValueError("include one thread as the numerical reference")
    depth = {"small": 2, "medium": 3, "large": 4}[args.size]
    environment = dict(os.environ, OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1",
                       VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1")
    records, summaries = [], []
    with tempfile.TemporaryDirectory(prefix="benchmark-", dir=args.workdir.resolve()) as temporary:
        root = Path(temporary)
        for name in dict.fromkeys(args.cases):
            directory = root / name
            directory.mkdir()
            geometry = fixture(name, depth)
            model = geometry.write(str(directory))
            # Fixed physical mesh and interaction setting, identical to the solver.
            output, _ = invoke([counter, model, "1e-34"], directory,
                               dict(environment, OMP_NUM_THREADS="1"))
            match = re.search(r"BENCHMARK panels (\d+) nodes (\d+) links (\d+)", output)
            if not match:
                raise ValueError("preflight returned no counts")
            panels, nodes, links = map(int, match.groups())
            if panels > 16000 or links > 8000000:
                raise ValueError("fixture exceeds panel/link budget")
            reference = None
            for repeat in range(args.repeats):
                # Interleave thread counts, rather than timing all serial runs first.
                for threads in sorted(set(args.threads)):
                    report = directory / "report.json"
                    command = [executable, "-b", model, "-m1e32", "-d1e-34", "-t1e-9",
                               "-pj", "-f0", "-i", "--report=" + str(report)]
                    output, elapsed = invoke(command, directory,
                                             dict(environment, OMP_NUM_THREADS=str(threads)))
                    data = json.loads(report.read_text())
                    report.unlink()
                    record = data["iterations"][-1]
                    if data["status"]["exit_code"] or record["solve_status"]["exit_code"] or not record["finite"]:
                        raise ValueError("incomplete or nonfinite solve")
                    if record["panels"] != panels or record["links"] != links:
                        raise ValueError("solve work differs from count-only preflight")
                    iterations = [len(re.findall(r"\d+", line)) for line in
                                  re.findall(r"GMRES Iteration:([^\n]*)", output)]
                    products = int(re.search(r"Hierarchical products: (\d+)", output)[1])
                    work = (record["matrix"], iterations, products)
                    if reference is None:
                        reference = work
                    elif work != reference:
                        raise ValueError("matrix or algebraic work changes across repeats/threads: " +
                                         repr({"reference": reference, "observed": work}))
                    row = dict(case=name, depth=2 if name == "uneven" else depth, threads=threads,
                               repeat=repeat, elapsed_seconds=elapsed, wall_seconds=timings(output),
                               panels=panels, nodes=nodes, links=links, gmres_steps=iterations,
                               hierarchical_products=products, matrix=record["matrix"], quality=record["quality"],
                               allocated_kib=int(re.search(r"Total allocated memory: (\d+) kilobytes", output)[1]),
                               potential_estimates=int(re.search(r"Number of potential estimates: (\d+)", output)[1]))
                    if geometry.analytic:
                        matrix = record["matrix"]
                        row["analytic_error"] = max(abs(complex(matrix["real"][i-1][j-1],
                            matrix["imaginary"][i-1][j-1])-value)/abs(value)
                            for (i,j),value in geometry.analytic.items())
                    records.append(row)
                    print("%s threads=%d repeat=%d elapsed=%.4fs" % (name, threads, repeat, elapsed), flush=True)
            serial = statistics.median(r["elapsed_seconds"] for r in records
                                       if r["case"] == name and r["threads"] == 1)
            for threads in sorted(set(args.threads)):
                rows = [r for r in records if r["case"] == name and r["threads"] == threads]
                elapsed = statistics.median(r["elapsed_seconds"] for r in rows)
                summaries.append(dict(case=name, threads=threads, median_seconds=elapsed,
                    speedup=serial/elapsed, wall_seconds={key:statistics.median(r["wall_seconds"][key] for r in rows)
                                                       for key in rows[0]["wall_seconds"]}))
    result = dict(build=data["build"], executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
                  size=args.size, records=records, summaries=summaries,
                  limits={"panels":16000,"links":8000000,"seconds_per_process":30},
                  memory_note="allocated_kib is solver accounting; sampled RSS and scratch bounds are in the outer guard receipt")
    with args.output.open("x") as output:
        json.dump(result, output, indent=2)
        output.write("\n")
    for row in summaries:
        print("%s %d threads: %.4fs, %.2fx" % (row["case"],row["threads"],row["median_seconds"],row["speedup"]))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, RuntimeError, KeyError, subprocess.TimeoutExpired) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
