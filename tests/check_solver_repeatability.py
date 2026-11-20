"""Check repeat/process/thread independence of a fixed dielectric solve."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import geometries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", required=True)
    args = parser.parse_args()
    executable = str(Path(args.fastercap).resolve())
    with tempfile.TemporaryDirectory(prefix="fc-repeat-") as temporary:
        root = Path(temporary)
        for permittivity in (100, 10-2j):
            model = geometries.coated_spheres_geometry(permittivity, depth=2).write(str(root))
            reference = None
            for threads in (1, 2, 4, 1, 2, 4):
                report = root / "report.json"
                process = subprocess.run([executable, "-b", model, "-m1e32", "-d1e-34",
                    "-t1e-9", "-pj", "-f0", "--report=" + str(report)],
                    capture_output=True, text=True, timeout=15,
                    env=dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_DYNAMIC="FALSE",
                             OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1"))
                assert process.returncode == 0, process.stdout[-4000:] + process.stderr
                data = json.loads(report.read_text()); report.unlink()
                record = data["iterations"][-1]
                assert data["status"]["exit_code"] == 0 and record["solve_status"]["exit_code"] == 0
                assert record["finite"] and record["panels"] == 960
                result = (record["matrix"], record["links"])
                if reference is None:
                    reference = result
                else:
                    assert result == reference, "matrix/work changes across repetitions or thread counts"
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
