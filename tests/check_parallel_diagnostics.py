"""Require stable diagnostic counts and full-precision matrices across workers."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import check_capacitance
import geometries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", required=True)
    parser.add_argument("--geometry", choices=("three_cubes", "coated_coax", "near_triangles"), required=True)
    args = parser.parse_args()
    executable = str(Path(args.fastercap).resolve())
    geometry = geometries.GEOMETRIES[args.geometry]
    options = ["-a0.01"] if args.geometry != "near_triangles" else ["-m1e32", "-pn"]
    with tempfile.TemporaryDirectory(prefix="fc-parallel-") as temporary:
        root = Path(temporary)
        model = geometry.write(str(root))
        reference = None
        for threads in (1, 2, 4):
            report = root / "report.json"
            command = [executable, "-b", model, *options, "-f0", "--report=" + str(report)]
            process = subprocess.run(command, capture_output=True, text=True, timeout=15,
                env=dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
                         OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1"))
            output = process.stdout + process.stderr
            if process.returncode:
                raise AssertionError(output[-4000:])
            data = json.loads(report.read_text()); report.unlink()
            assert data["status"]["exit_code"] == 0
            counts = [int(n) for n in re.findall(r"Number of potential estimates: (\d+)", output)]
            matrices = [r["matrix"] for r in data["iterations"]]
            assert counts and all(n > 0 for n in counts)
            for record in data["iterations"]:
                assert record["finite"] and record["solve_status"]["exit_code"] == 0
            matrix = matrices[-1]
            rows = [[complex(a,b) for a,b in zip(real,imag)]
                    for real,imag in zip(matrix["real"],matrix["imaginary"])]
            # The 2D report retains the reference conductor column; the physical
            # checker expects the reduced square matrix used by console output.
            if matrix["reference_column"] is not None:
                rows = [row[:len(rows)] for row in rows]
            errors = check_capacitance.check_matrix(rows, geometry, 0.01)
            assert not errors, "\n".join(errors)
            if args.geometry == "near_triangles":
                assert counts == [2], counts
                assert output.count("Warning: panel distance too small found during potential calculation") == 1
            result = (matrices, counts)
            print("%d threads: potential estimates %s" % (threads, counts))
            if reference is None:
                reference = result
            else:
                assert counts == reference[1], "potential counts differ: %s vs %s" % (counts, reference[1])
                assert matrices == reference[0], "full-precision matrices differ"
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
