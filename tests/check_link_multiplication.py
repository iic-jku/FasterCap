"""Preserve full-precision products across workers and genuinely loaded link blocks."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import geometries


def run(executable, directory, model, workers, ratio):
    report = directory / "report.json"
    process = subprocess.run([executable, "-b", model, "-m1e32", "-d1e-36", "-t1e-9",
                              "-pj", "-f" + str(ratio), "--report=" + str(report)],
                             cwd=directory, capture_output=True, text=True, timeout=15,
                             env=dict(os.environ, OMP_NUM_THREADS=str(workers), OMP_DYNAMIC="FALSE",
                                      OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1"))
    output = process.stdout + process.stderr
    if ratio == 1e12:
        # Deliberately insufficient storage: read the allocation estimate before
        # any link array is allocated, then choose an 18 MiB one-chunk budget.
        free = re.search(r"Available free memory left is: (\d+) Mbytes", output)
        required = re.search(r"links is: (\d+) Mbytes", output)
        assert free and required and int(required[1]) > 16, output[-4000:]
        if report.exists():
            report.unlink()
        return int(free[1]) / 18.0
    assert process.returncode == 0, output[-4000:]
    data = json.loads(report.read_text())
    report.unlink()
    record = data["iterations"][-1]
    assert data["status"]["exit_code"] == 0 and record["solve_status"]["exit_code"] == 0
    assert record["finite"] and record["panels"] == 1280
    assert 2**20 < record["links"] <= 1280 * 1279
    assert ("Going out-of-core" in output) == (ratio != 0), output[-4000:]
    steps = [len(re.findall(r"\d+", line)) for line in re.findall(r"GMRES Iteration:([^\n]*)", output)]
    return record["matrix"], record["links"], steps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", required=True)
    args = parser.parse_args()
    executable = str(Path(args.fastercap).resolve())
    with tempfile.TemporaryDirectory(prefix="fc-link-products-") as temporary:
        root = Path(temporary)
        geometry = geometries.coated_spheres_geometry(100, depth=2, split_coat=True,
                                                       offset=(.4, -.2, .15))
        model = geometry.write(str(root))
        # The fixed triangular mesh bounds all possible directed links before
        # invoking the solver: 1280*1279 pairs, below 26 MiB of link payload.
        triangles = sum(line.startswith("T ") for path in root.glob("*.txt")
                        for line in path.read_text().splitlines())
        assert triangles == 1280
        reference = run(executable, root, model, 1, 0)
        for workers in (2, 4):
            assert run(executable, root, model, workers, 0) == reference
        for workers in (1, 2, 4):
            ratio = run(executable, root, model, workers, 1e12)
            assert run(executable, root, model, workers, ratio) == reference
    print("Identical matrices and work at 1/2/4 workers, in-core and across loaded blocks")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
