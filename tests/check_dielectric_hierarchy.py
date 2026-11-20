"""Compare compressed dielectric solves with the same fixed mesh using direct links."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import check_capacitance
import geometries


def solve(executable, directory, geometry, distance):
    directory.mkdir(parents=True, exist_ok=True)
    model = geometry.write(str(directory))
    report = directory / "report.json"
    command = [executable, "-b", model, "-m1e32", "-d" + distance,
               "-t1e-9", "-pj", "-f0", "--report=" + str(report)]
    process = subprocess.run(command, capture_output=True, text=True, timeout=20,
                             env=dict(os.environ, OMP_NUM_THREADS="1"))
    (directory / "solver.log").write_text(process.stdout + process.stderr)
    if process.returncode:
        raise AssertionError(process.stdout[-4000:] + process.stderr)
    data = json.loads(report.read_text())
    iteration = data["iterations"][-1]
    assert data["status"]["exit_code"] == 0 and iteration["solve_status"]["exit_code"] == 0
    assert iteration["finite"] and iteration["panels"] == 960
    matrix = iteration["matrix"]
    rows = [[complex(a, b) for a, b in zip(real, imaginary)]
            for real, imaginary in zip(matrix["real"], matrix["imaginary"])]
    errors = check_capacitance.check_matrix(rows, geometry, 0.02, reciprocity_tolerance=0.01)
    if errors:
        raise AssertionError("\n".join(errors))
    return rows, iteration["links"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastercap", required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--high-contrast", action="store_true")
    args = parser.parse_args()
    root = args.workdir.resolve()
    executable = str(Path(args.fastercap).resolve())
    root.mkdir(parents=True, exist_ok=True)
    # Numerical reports intentionally refuse to overwrite an existing file.
    with tempfile.TemporaryDirectory(prefix="run-", dir=root) as temporary:
        root = Path(temporary)
        cases = [("concentric", 10.0, (0, 0, 0), 2.0),
                 ("offset", 10.0, (0.4, -0.2, 0.15), 2.0),
                 ("lossy", 10.0-2.0j, (0.4, -0.2, 0.15), 2.0)]
        if args.high_contrast:
            cases = [("concentric100", 100.0, (0, 0, 0), 2.0),
                     ("offset100", 100.0, (0.4, -0.2, 0.15), 2.0),
                     ("thin100", 100.0, (0, 0, 0), 1.125)]
        matrix_tolerance = 0.015 if args.high_contrast else 0.005
        for name, permittivity, offset, radius in cases:
            geometry = geometries.coated_spheres_geometry(permittivity, depth=2,
                                                          offset=offset, coat_radius=radius)
            reference, full_links = solve(executable, root / name / "direct", geometry, "1e-40")
            assert full_links == 960 * 959
            for distance in ("1e-35", "1e-36"):
                matrix, links = solve(executable, root / name / distance, geometry, distance)
                assert 0 < links < full_links, "test must exercise grouped interactions"
                error = max(abs(a-b)/abs(b) for row, ref in zip(matrix, reference) for a, b in zip(row, ref))
                print("%s -d%s: %d links, %.6g%% difference from direct" %
                      (name, distance, links, 100*error))
                assert error <= matrix_tolerance, "hierarchy exceeds its matrix-error budget"
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, OSError, ValueError, KeyError, subprocess.TimeoutExpired) as error:
        print("FAILED:", error, file=sys.stderr)
        sys.exit(1)
