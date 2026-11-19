"""Validate report serialization and its CLI contract with Python's JSON parser."""

import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import check_capacitance
import geometries


def read_report(path):
    def invalid_constant(value):
        raise AssertionError("non-JSON numerical constant: " + value)
    data = json.loads(path.read_text(), parse_constant=invalid_constant)
    assert data["schema_version"] == 1
    assert data["last_completed_iteration"] == (len(data["iterations"]) - 1 if data["iterations"] else None)
    build = data["build"]
    assert build["compiler"] and build["configuration"]
    assert set(build["source_revisions"]) == {"FasterCap", "LinAlgebra", "Geometry"}
    assert all(value is None or len(value) == 40 for value in build["source_revisions"].values())
    assert all(value in (True, False, None) for value in build["source_dirty"].values())
    return data


def fixture_checks(writer, directory):
    subprocess.run([writer, str(directory)], check=True)
    real = read_report(directory / "real.json")
    matrix = real["iterations"][0]["matrix"]
    expected = [[math.nextafter(1.0, 2.0), -0.12345678901234567], [-0.125, sys.float_info.min]]
    assert matrix["real"] == expected, "finite doubles did not round trip"
    assert matrix["imaginary"] == [[0.0, 0.0], [0.0, 0.0]]
    assert math.copysign(1.0, matrix["imaginary"][1][1]) == -1.0
    assert real["input_file"] == 'input\\with"quotes\tand\ncontrols.lst'
    assert matrix["row_conductors"] == ['first\\"\t\n\r\b\f\x01', "second"]
    quality = real["iterations"][0]["quality"]
    assert quality["negative_row_sum_count"] == 1
    assert quality["nonpositive_diagonal_count"] == 0
    assert quality["nonnegative_offdiagonal_count"] == 0
    assert quality["worst_reciprocal_pair"]["row"] == 1
    assert quality["worst_reciprocal_pair"]["column"] == 2

    lossy = read_report(directory / "complex.json")
    assert lossy["status"]["exit_code"] == 75
    imaginary = [[-math.nextafter(0.25, 1.0), 0.0625], [0.03125, -0.25]]
    assert lossy["iterations"][0]["matrix"]["imaginary"] == imaginary
    matrix = [[complex(a, b) for a, b in zip(row_a, row_b)] for row_a, row_b in zip(expected, imaginary)]
    defect = math.sqrt(sum(abs(matrix[i][j] - matrix[j][i]) ** 2 for i in range(2) for j in range(2)))
    norm = math.sqrt(sum(abs(value) ** 2 for row in matrix for value in row))
    assert math.isclose(lossy["iterations"][0]["quality"]["reciprocity_frobenius"], defect / norm, rel_tol=1e-14)

    two_d = read_report(directory / "2d.json")["iterations"][0]
    scale = 5.563250279878e-11  # TWO_PI_TIMES_E0, physical 2D output scale
    assert two_d["matrix"]["real"] == [[2 * scale, -2 * scale]]
    assert two_d["matrix"]["imaginary"] == [[-0.5 * scale, 0.5 * scale]]
    assert two_d["matrix"]["reference_column"] == 2
    assert two_d["matrix"]["units"] == "F/m"
    assert two_d["quality"]["worst_reciprocal_pair"] is None
    assert two_d["quality"]["row_sums_real"] == [2 * scale]

    nonfinite = read_report(directory / "nonfinite.json")
    assert nonfinite["status"]["exit_code"] == 65
    assert not nonfinite["iterations"][0]["finite"]
    assert nonfinite["iterations"][0]["quality"] is None
    assert nonfinite["iterations"][0]["matrix"]["real"] == [[None, None]]
    assert nonfinite["iterations"][0]["matrix"]["imaginary"][0][0] is None
    empty = read_report(directory / "empty.json")
    assert empty["iterations"] == [] and empty["status"]["exit_code"] == 66


def run(executable, directory, options, expected=0):
    result = subprocess.run([executable, "-b"] + options, cwd=directory, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, timeout=40)
    assert result.returncode == expected, "exit %d, expected %d:\n%s" % (result.returncode, expected, result.stdout)
    return result.stdout


def compare_printed(record, output):
    names, printed = check_capacitance.parse_last_matrix(output)
    matrix = record["matrix"]
    n = len(matrix["row_conductors"])
    assert names == matrix["row_conductors"]
    for i in range(n):
        for j in range(n):
            raw = complex(matrix["real"][i][j], matrix["imaginary"][i][j])
            assert check_capacitance.relative_difference(raw, printed[i][j]) <= 5e-6


def solver_checks(executable, directory):
    geometry = geometries.GEOMETRIES["three_cubes"]
    path = Path(geometry.write(str(directory)))
    input_bytes = path.read_bytes()
    baseline = run(executable, directory, [path.name, "-e"])
    csv_path = path.with_suffix(".csv")
    csv_bytes = csv_path.read_bytes()
    report_path = directory / "report with spaces.json"
    output = run(executable, directory, [path.name, "-e", "--report=" + report_path.name])
    report = read_report(report_path)
    assert report["status"]["exit_code"] == 0
    assert len(report["iterations"]) == 1
    assert check_capacitance.parse_last_matrix(baseline) == check_capacitance.parse_last_matrix(output)
    assert csv_path.read_bytes() == csv_bytes, "report changed legacy CSV"
    compare_printed(report["iterations"][-1], output)

    output = run(executable, directory, [path.name, "-a0.01", "--report=auto.json"])
    report = read_report(directory / "auto.json")
    assert len(report["iterations"]) == output.count("Capacitance matrix is:") > 1
    for number, record in enumerate(report["iterations"]):
        assert record["iteration"] == number and record["panels"] > 0 and record["links"] > 0
    compare_printed(report["iterations"][-1], output)

    output = run(executable, directory, [path.name, "-t1e-30", "--report=unconverged.json"], 75)
    report = read_report(directory / "unconverged.json")
    assert report["status"]["exit_code"] == 75
    assert report["iterations"][-1]["solve_status"]["exit_code"] == 75
    compare_printed(report["iterations"][-1], output)

    for name in ("coated_sphere_lossy", "coated_coax"):
        case = directory / name
        case.mkdir(exist_ok=True)
        input_path = Path(geometries.GEOMETRIES[name].write(str(case)))
        output = run(executable, case, [input_path.name, "-a0.01", "--report=result.json"])
        report = read_report(case / "result.json")
        compare_printed(report["iterations"][-1], output)
        if name.endswith("lossy"):
            assert any(v != 0 for row in report["iterations"][-1]["matrix"]["imaginary"] for v in row)

    run(executable, directory, ["missing.lst", "--report=missing.json"], 66)
    report = read_report(directory / "missing.json")
    assert report["iterations"] == [] and report["status"]["exit_code"] == 66
    invalid_path = Path(geometries.GEOMETRIES["degenerate_triangle"].write(str(directory)))
    run(executable, directory, [invalid_path.name, "--report=invalid.json"], 65)
    report = read_report(directory / "invalid.json")
    assert report["status"]["exit_code"] == 65
    if report["iterations"]:
        assert not report["iterations"][-1]["finite"]

    for options in (["--report="], ["--report"], ["--report=a.json", "--report=b.json"],
                    ['--report="unterminated']):
        run(executable, directory, [path.name] + options, 64)
    run(executable, directory, [path.name, "--report=" + "x" * 5000], 74)
    # The former 256-byte temporary overflowed on a valid report path of this length.
    long_directory = directory / ("a" * 130) / ("b" * 130)
    long_directory.mkdir(parents=True, exist_ok=True)
    long_report = long_directory / "result.json"
    run(executable, directory, [path.name, "--report=" + str(long_report.relative_to(directory))])
    assert read_report(long_report)["status"]["exit_code"] == 0
    run(executable, directory, [path.name, "--report=nonexistent/output.json"], 74)
    run(executable, directory, [path.name, "--report=" + path.name], 74)
    assert path.read_bytes() == input_bytes, "report overwrote primary input"
    existing = directory / "existing.json"
    existing.write_bytes(b"existing report must survive")
    panel = Path(check_capacitance.panel_files(str(path))[0])
    panel_bytes = panel.read_bytes()
    hardlink = directory / "panel-hardlink.json"
    os.link(panel, hardlink)
    targets = [(existing, existing.read_bytes()), (panel, panel_bytes), (hardlink, panel_bytes)]
    symlink = directory / "panel-symlink.json"
    try:
        symlink.symlink_to(panel)
    except OSError as error:
        if os.name != "nt" or getattr(error, "winerror", None) != 1314:
            raise
        print("symbolic-link check unavailable: Windows requires link creation privilege")
    else:
        targets.append((symlink, panel_bytes))
    for target, original in targets:
        run(executable, directory, [path.name, "--report=" + str(target)], 74)
        assert target.read_bytes() == original, "report modified existing file or alias"
        assert panel.read_bytes() == panel_bytes, "report overwrote referenced panel input"



def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--writer")
    mode.add_argument("--fastercap")
    parser.add_argument("--workdir", required=True, type=Path)
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="run-", dir=args.workdir) as scratch:
        if args.writer:
            fixture_checks(os.path.abspath(args.writer), Path(scratch))
        else:
            solver_checks(os.path.abspath(args.fastercap), Path(scratch))
    print("numerical report checks passed")


if __name__ == "__main__":
    main()
