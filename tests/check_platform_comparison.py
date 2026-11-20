"""Capture the original #10 plate and compare full-precision platform reports.

The budgets detect numerical portability regressions; they are not physical error
bounds. No console-precision fallback or stored platform matrix is accepted.
"""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import geometries

THREADS = dict(OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1", BLIS_NUM_THREADS="1")
FIXED_OPTIONS = ["-m1e32", "-d1e-50", "-t1e-10", "-pj", "-f0"]
AUTO_OPTIONS = ["-a0.05", "-d0.5", "-m0.5", "-ap", "-f0"]
# Tight solves have a 1e-10 algebraic target. A 1e-7 cross-build investigation
# budget allows amplification without accepting the historical percent-level drift.
FIXED_BUDGET = 1e-7
# One tenth of the original requested 5% automatic target: detects the observed
# 3% path divergence without requiring all compilers to round identically.
ADAPTIVE_BUDGET = 0.005
CASES = ("fixed_mesh", "adaptive", "manual_replay")
FIXED_PANELS = 1404
FIXED_LINK_BOUND = FIXED_PANELS * (FIXED_PANELS - 1)
MAX_PANELS = 16000
MAX_LINKS = 8000000
MAX_OUTPUT_BYTES = 2 * 1024 * 1024


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_json(path):
    def invalid(value):
        raise ValueError("nonstandard JSON number: " + value)
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=invalid)


def matrix_values(matrix):
    rows, columns = matrix["row_conductors"], matrix["column_conductors"]
    if len(rows) != 2 or rows != columns or len(set(rows)) != 2 \
            or matrix["units"] != "F" or matrix["reference_column"] is not None:
        raise ValueError("expected the plate's full 2 by 2 matrix in solver farads")
    arrays = [matrix["real"], matrix["imaginary"]]
    for array in arrays:
        if len(array) != 2 or any(len(row) != 2 for row in array):
            raise ValueError("invalid matrix dimensions")
        if any(isinstance(v, bool) or not isinstance(v, (float, int)) or not math.isfinite(v)
               for row in array for v in row):
            raise ValueError("invalid matrix value")
    return [complex(arrays[0][i][j], arrays[1][i][j]) for i in range(2) for j in range(2)]


def validate_report(report):
    records = report["iterations"]
    if report["schema_version"] != 1 or report["status"]["exit_code"] != 0 or not records \
            or report["last_completed_iteration"] != len(records) - 1:
        raise ValueError("missing successful full-precision result")
    for i, record in enumerate(records):
        if record["iteration"] != i or record["solve_status"]["exit_code"] != 0 or record["finite"] is not True:
            raise ValueError("invalid or unconverged iteration")
        if not 0 < record["panels"] <= MAX_PANELS or not 0 < record["links"] <= MAX_LINKS:
            raise ValueError("plate capture exceeded panel/link investigation caps")
        matrix_values(record["matrix"])
    return records[-1]


def work(record):
    # Matching these fields is evidence of matching work, not a topology hash.
    return {key: record[key] for key in ("panels", "links", "settings")}


def matrix_difference(a, b):
    if a["row_conductors"] != b["row_conductors"]:
        raise ValueError("conductor ordering differs")
    av, bv = matrix_values(a), matrix_values(b)
    differences = [abs(x - y) for x, y in zip(av, bv)]
    relative = [d / max(abs(x), abs(y)) if x or y else 0.0
                for d, x, y in zip(differences, av, bv)]
    scale = max(math.sqrt(sum(abs(x)**2 for x in av)), math.sqrt(sum(abs(x)**2 for x in bv)))
    return {"max_entry_relative": max(relative), "max_entry_absolute": max(differences),
            "frobenius_relative": math.sqrt(sum(d*d for d in differences)) / scale if scale else 0.0,
            "entry_relative": relative, "entry_absolute": differences}


def host_metadata():
    data = {"platform": platform.platform(), "machine": platform.machine(),
            "processor": platform.processor(), "thread_environment": THREADS}
    commands = {}
    if sys.platform == "darwin":
        commands = {"cpu": ["sysctl", "-n", "machdep.cpu.brand_string"]}
    elif sys.platform.startswith("linux"):
        commands = {"cpu": ["lscpu"], "packages": ["dpkg-query", "-W", "libopenblas*", "libwx*", "libgomp*"]}
    elif os.name == "nt":
        data["processor_identifier"] = os.environ.get("PROCESSOR_IDENTIFIER")
    for name, command in commands.items():
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            data[name] = {"command": command, "status": result.returncode, "output": result.stdout[:16000]}
        except (OSError, subprocess.SubprocessError) as error:
            data[name] = {"command": command, "error": str(error)}
    return data


def replay_options(record):
    settings = record["settings"]
    supported = {"automatic", "auto_tolerance", "gmres_tolerance", "mesh_refinement",
                 "interaction_ratio", "curvature", "scheme", "preconditioner"}
    if set(settings) != supported:
        raise ValueError("manual replay cannot represent automatic settings: " +
                         ", ".join(sorted(set(settings) - supported)) +
                         "; use the exact exported refined mesh for this algorithm")
    if settings["scheme"] != "collocation":
        raise ValueError("unexpected plate discretization scheme")
    return ["-m%.17g" % settings["mesh_refinement"], "-d%.17g" % settings["interaction_ratio"],
            "-mc%.17g" % settings["curvature"], "-t1e-10", "-pj", "-f0"]



def run_solver(command, destination, report_path, log_path, timeout):
    started = time.monotonic()
    paths = (report_path, log_path)
    oversized = False
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT,
                                   env=dict(os.environ, **THREADS), cwd=destination)
        try:
            while process.poll() is None:
                if time.monotonic() - started > timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
                if any(path.exists() and path.stat().st_size > MAX_OUTPUT_BYTES for path in paths):
                    raise ValueError("solver report/log exceeded 2 MiB; retained prefix only")
                try:
                    process.wait(timeout=0.02)
                except subprocess.TimeoutExpired:
                    pass
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            for path in paths:
                if path.exists() and path.stat().st_size > MAX_OUTPUT_BYTES:
                    with path.open("r+b") as stream:
                        stream.truncate(MAX_OUTPUT_BYTES)
                    oversized = True
    if oversized:
        raise ValueError("solver report/log exceeded 2 MiB; retained prefix only")
    return process.returncode, time.monotonic() - started


def capture(executable, destination, timeout, fixed_only=False):
    destination.mkdir(parents=True, exist_ok=False)
    models, hashes = {}, {}
    for name, triangles in (("fixed_triangles", True), ("original_quads", False)):
        model_dir = destination / name
        model_dir.mkdir()
        models[name] = Path(geometries.write_plate_over_substrate(str(model_dir), triangles)).resolve()
        hashes[name] = {p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                        for p in sorted(model_dir.iterdir())}
    artifact = {"schema_version": 1, "host": host_metadata(), "geometry_sha256": hashes,
                "physical_capacitance_scale": 1e-6, "cases": {}}
    try:
        for name in CASES[:1] if fixed_only else CASES:
            if name == "fixed_mesh":
                options = FIXED_OPTIONS
            elif name == "adaptive":
                options = AUTO_OPTIONS
            else:
                options = replay_options(artifact["cases"]["adaptive"]["report"]["iterations"][-1])
            geometry = "fixed_triangles" if name == "fixed_mesh" else "original_quads"
            model = models[geometry]
            report_path = destination / (name + ".json")
            command = [str(executable), "-b", str(model), *options, "--report=" + str(report_path.resolve())]
            # Stream logs so a timeout or outer guard leaves inspectable evidence.
            returncode, elapsed = run_solver(command, destination, report_path, destination / (name + ".log"), timeout)
            record = {"geometry": geometry, "options": options, "returncode": returncode, "elapsed_seconds": elapsed}
            artifact["cases"][name] = record
            if returncode != 0:
                raise ValueError(name + " solver returned " + str(returncode))
            report = read_json(report_path)
            final = validate_report(report)
            record["report"] = report
            if name == "fixed_mesh" and ((len(report["iterations"]), final["panels"]) != (1, FIXED_PANELS) or final["links"] > FIXED_LINK_BOUND):
                raise ValueError("fixed triangle panel count or link budget changed")
            if name == "manual_replay":
                adaptive = artifact["cases"]["adaptive"]["report"]["iterations"][-1]
                record["matches_adaptive_counts"] = all(final[key] == adaptive[key] for key in ("panels", "links"))
                record["difference_from_adaptive"] = matrix_difference(final["matrix"], adaptive["matrix"])
                if not record["matches_adaptive_counts"]:
                    raise ValueError("manual settings did not reproduce automatic mesh/link counts; exact mesh export required")
            write_json(destination / "capture.json", artifact)
            print(name, final["panels"], final["links"], "panels/links", "%.3fs" % elapsed)
    except (ValueError, KeyError, TypeError, OSError, subprocess.SubprocessError) as error:
        artifact["error"] = str(error)
        write_json(destination / "capture.json", artifact)
        print("FAILED:", error, file=sys.stderr)
        return 1
    return 0


def compare(artifacts, labels, allow_dirty=False):
    if len(artifacts) < 2 or len(labels) != len(artifacts):
        raise ValueError("at least two labeled platform captures are required")
    sources = None
    for artifact in artifacts:
        if artifact["schema_version"] != 1 or "error" in artifact or set(artifact["cases"]) != set(CASES):
            raise ValueError("incomplete or failed capture")
        if artifact["geometry_sha256"] != artifacts[0]["geometry_sha256"]:
            raise ValueError("input geometry hashes differ")
        for name, case in artifact["cases"].items():
            if case["returncode"] != 0:
                raise ValueError("failed solver capture")
            final = validate_report(case["report"])
            build = case["report"]["build"]
            for key in ("compiler", "configuration", "target_system", "target_processor", "wxwidgets", "blas_link_libraries"):
                if not isinstance(build.get(key), str) or not build[key]:
                    raise ValueError("missing build provenance: " + key)
            if build.get("omp_max_threads") != 1 or any(
                    build.get("thread_environment", {}).get(key) != "1"
                    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")):
                raise ValueError("capture did not use the controlled thread environment")
            revisions = build["source_revisions"]
            if set(revisions) != {"FasterCap", "LinAlgebra", "Geometry"} or not all(revisions.values()):
                raise ValueError("missing source revisions")
            if not allow_dirty and any(build["source_dirty"].get(k) is not False for k in revisions):
                raise ValueError("dirty or unknown build sources")
            if sources is None:
                sources = revisions
            elif sources != revisions:
                raise ValueError("build source revisions differ")
            expected = FIXED_OPTIONS if name == "fixed_mesh" else AUTO_OPTIONS if name == "adaptive" else \
                replay_options(artifact["cases"]["adaptive"]["report"]["iterations"][-1])
            if case.get("geometry") != ("fixed_triangles" if name == "fixed_mesh" else "original_quads"):
                raise ValueError("wrong plate discretization for capture case")
            if case["options"] != expected:
                raise ValueError("unexpected solver options")
            if name != "adaptive" and (final["settings"]["automatic"] or final["settings"]["gmres_tolerance"] != 1e-10):
                raise ValueError("fixed/manual capture did not use tight manual GMRES")
            if name == "fixed_mesh" and ((len(case["report"]["iterations"]), final["panels"]) != (1, FIXED_PANELS) or final["links"] > FIXED_LINK_BOUND):
                raise ValueError("fixed mesh panel count or link budget changed")
    comparisons = []
    for left, right in itertools.combinations(range(len(artifacts)), 2):
        for name in CASES:
            a, b = [artifacts[k]["cases"][name]["report"]["iterations"] for k in (left, right)]
            same_work = work(a[-1]) == work(b[-1])
            budget = FIXED_BUDGET if name == "fixed_mesh" or (name == "manual_replay" and same_work) else ADAPTIVE_BUDGET
            difference = matrix_difference(a[-1]["matrix"], b[-1]["matrix"])
            comparisons.append({"left": labels[left], "right": labels[right], "case": name,
                                "relative_budget": budget, "same_reported_final_work": same_work,
                                "same_reported_path": [work(v) for v in a] == [work(v) for v in b],
                                "iteration_counts": [len(a), len(b)], "difference": difference,
                                "passed": difference["max_entry_relative"] <= budget and (name != "fixed_mesh" or same_work)})
    return {"schema_version": 1, "source_revisions": sources,
            "passed": all(item["passed"] for item in comparisons), "comparisons": comparisons}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    make = sub.add_parser("capture")
    make.add_argument("--fastercap", required=True, type=Path)
    make.add_argument("--output-dir", required=True, type=Path)
    make.add_argument("--timeout", type=float, default=60)
    make.add_argument("--fixed-only", action="store_true", help="bounded local first probe; not a complete CI artifact")
    diff = sub.add_parser("compare")
    diff.add_argument("artifacts", nargs="+", type=Path)
    diff.add_argument("--output", required=True, type=Path)
    diff.add_argument("--allow-dirty", action="store_true", help="controlled local comparisons only")
    args = parser.parse_args()
    try:
        if args.command == "capture":
            if not math.isfinite(args.timeout) or not 0 < args.timeout <= 120:
                raise ValueError("timeout must be positive and at most 120 seconds per solve")
            return capture(args.fastercap.resolve(), args.output_dir.resolve(), args.timeout, args.fixed_only)
        result = compare([read_json(path) for path in args.artifacts], [str(path.parent) for path in args.artifacts], args.allow_dirty)
        write_json(args.output, result)
        for item in result["comparisons"]:
            print("%s / %s %s: %s, max entry difference %.9g, budget %.9g, matching work %s" %
                  (item["left"], item["right"], item["case"], "PASS" if item["passed"] else "FAIL",
                   item["difference"]["max_entry_relative"], item["relative_budget"], item["same_reported_final_work"]))
        return 0 if result["passed"] else 1
    except (ValueError, TypeError, KeyError, IndexError, OSError) as error:
        print("FAILED:", error, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
