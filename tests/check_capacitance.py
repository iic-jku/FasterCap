"""Run FasterCap on a test geometry and check the capacitance matrix it prints.

Usage: check_capacitance.py --fastercap <exe> --geometry <name> --tolerance <rel>
                            --workdir <dir> -- <FasterCap options>

Checks, for the last matrix printed:
- the exit code (0 unless --exit-code), and one row per conductor of the geometry
- finite entries, positive diagonal, negative off-diagonal terms
- symmetry, and row sums not below zero (capacitance to infinity), within the tolerance
- the analytic values and the entries that are equal by symmetry of the geometry

Complex permittivities give complex entries: their real parts must have the signs above, and
the differences are relative to the magnitudes.

With --out-of-core, FasterCap runs again out-of-core with blocks of one link chunk, and the
two matrices must agree.

With --no-matrix, FasterCap must stop without printing a matrix.

With --reference, FasterCap runs again on another geometry made of the same panels, and the
two matrices must agree.

With --names, the conductors must have the given names.

With --same-as, FasterCap runs again with other options, which must give the same matrix.

With --dump-input, FasterCap writes the input geometry back (-oi), and the dump must give the
same matrix.

With --output-files, the run also writes the refined geometry (-o), the matrix (-e) and the
charge densities (-c): the refined geometry must hold all refined panels, the matrix file the
printed matrix, and the charges of each excitation must add up to its row of the matrix.

With --interrupt, an out-of-core run gets SIGINT once it has written link chunks: it must stop
with FC_USER_BREAK and delete its temporary files.
"""

import argparse
import cmath
import math
import os
import re
import signal
import subprocess
import sys
import time

# keep the source tree free of __pycache__
sys.dont_write_bytecode = True

import geometries  # noqa: E402

OUTPUT_TAIL_LINES = 60

# a link chunk holds 2^20 links (AUTOREFINE_LINK_CHUNK_SIZE) of 16 bytes; blocks of 24 Mbytes
# hold one chunk each, also if the free memory changes by a third between the runs
CHUNK_MBYTES = 16
BLOCK_MBYTES = 24
# in-core and out-of-core runs add the same terms in the same order: the printed digits agree
OUT_OF_CORE_TOLERANCE = 1e-5
# the same panels, given differently
REFERENCE_TOLERANCE = 1e-5
# files written next to the console output: the same digits, charges with 7 of them
OUTPUT_TOLERANCE = 1e-5
# below the CTest timeout, so that a hanging run still reports its output
RUN_TIMEOUT = 240
# FC_USER_BREAK
USER_BREAK_EXIT_CODE = 125
TEMPORARY_FILE = re.compile(r"frc[lp][0-9A-F]+\.tmp$")


def parse_last_matrix(output):
    blocks = output.split("Capacitance matrix is:")
    if len(blocks) < 2:
        return None, []
    lines = blocks[-1].strip().splitlines()
    m = re.match(r"Dimension (\d+) x (\d+)", lines[0]) if lines else None
    if not m:
        return None, []
    n = int(m.group(1))
    names, rows = [], []
    for line in lines[1:1 + n]:
        fields = line.split()
        names.append(fields[0])
        rows.append([number(v) for v in fields[1:]])
    return names, rows


def number(text):
    """a matrix entry: real, or complex as '1.8e-10-9.9e-12j' for complex permittivities"""
    value = complex(text)
    return value.real if value.imag == 0 else value


def relative_difference(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale else 0.0


def check_matrix(rows, geometry, tolerance):
    errors = []
    n = len(rows)
    if n != geometry.conductors or any(len(r) != n for r in rows):
        return ["expected a %d x %d matrix, got %d rows" % (geometry.conductors, geometry.conductors, n)]
    for i in range(n):
        for j in range(n):
            if not cmath.isfinite(rows[i][j]):
                errors.append("C%d%d is %r" % (i + 1, j + 1, rows[i][j]))
    if errors:
        return errors
    for i in range(n):
        if rows[i][i].real <= 0:
            errors.append("C%d%d = %g is not positive" % (i + 1, i + 1, rows[i][i]))
        if sum(rows[i]).real < -tolerance * rows[i][i].real:
            errors.append("row %d sums to %s, below zero" % (i + 1, sum(rows[i])))
        for j in range(n):
            if i == j:
                continue
            if rows[i][j].real >= 0:
                errors.append("C%d%d = %g is not negative" % (i + 1, j + 1, rows[i][j]))
            if j > i and relative_difference(rows[i][j], rows[j][i]) > tolerance:
                errors.append("C%d%d = %g and C%d%d = %g differ by more than %g"
                              % (i + 1, j + 1, rows[i][j], j + 1, i + 1, rows[j][i], tolerance))
    for (i, j), expected in geometry.analytic.items():
        obtained = rows[i - 1][j - 1]
        if relative_difference(obtained, expected) > tolerance:
            errors.append("C%d%d = %s, analytic %s (%.2f %% off)"
                          % (i, j, obtained, expected, 100.0 * relative_difference(obtained, expected)))
    for (i, j), (k, l) in geometry.equal_entries:
        if relative_difference(rows[i - 1][j - 1], rows[k - 1][l - 1]) > tolerance:
            errors.append("C%d%d = %g and C%d%d = %g should be equal"
                          % (i, j, rows[i - 1][j - 1], k, l, rows[k - 1][l - 1]))
    return errors


def run(fastercap, list_file, options):
    command = [os.path.abspath(fastercap), os.path.basename(list_file), "-b"] + options
    try:
        result = subprocess.run(command, cwd=os.path.dirname(list_file), stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, universal_newlines=True,
                                timeout=RUN_TIMEOUT)
    except subprocess.TimeoutExpired as exc:
        output = exc.output or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return command, "timeout after %d s" % RUN_TIMEOUT, output
    return command, result.returncode, result.stdout


def check_run(command, returncode, output, geometry, tolerance, exit_code=0, no_matrix=False):
    errors = []
    if isinstance(returncode, str):
        errors.append(returncode)
    elif returncode != exit_code:
        errors.append("exit code %d, expected %d" % (returncode, exit_code))
    names, rows = parse_last_matrix(output)
    if no_matrix:
        if names is not None:
            errors.append("a capacitance matrix in the output")
    elif names is None:
        errors.append("no capacitance matrix in the output")
    else:
        errors += check_matrix(rows, geometry, tolerance)
    print(" ".join(command))
    for name, row in zip(names or [], rows):
        print("%-12s %s" % (name, " ".join("%.6g" % v if isinstance(v, float) else str(v) for v in row)))
    return rows, errors


def panel_files(list_file):
    """the panel files a list file names in its 'C' lines, in their order"""
    directory = os.path.dirname(list_file)
    with open(list_file) as f:
        return [os.path.join(directory, line.split()[1]) for line in f if line.startswith("C ")]


def read_panels(path):
    """(vertices, value after them) of the T and Q panels of a panel file"""
    panels = []
    with open(path) as f:
        for line in f:
            fields = line.split()
            if fields and fields[0] in ("T", "Q"):
                corners = 3 if fields[0] == "T" else 4
                values = [float(v) for v in fields[2:]]
                vertices = [values[3 * k:3 * k + 3] for k in range(corners)]
                rest = values[3 * corners:]
                panels.append((vertices, rest[0] if len(rest) == 1 else None))
    return panels


def area(vertices):
    total = 0.0
    for k in range(1, len(vertices) - 1):
        u = [vertices[k][m] - vertices[0][m] for m in range(3)]
        w = [vertices[k + 1][m] - vertices[0][m] for m in range(3)]
        total += 0.5 * math.sqrt((u[1] * w[2] - u[2] * w[1]) ** 2 + (u[2] * w[0] - u[0] * w[2]) ** 2
                                 + (u[0] * w[1] - u[1] * w[0]) ** 2)
    return total


def check_output_files(list_file, output, names, rows):
    """errors in the files of -o (refined geometry), -e (matrix) and -c (charge densities)"""
    errors = []
    base = os.path.splitext(list_file)[0]
    n = len(rows)

    refined = re.findall(r"Number of panels after refinement: (\d+)", output)
    refined_file = base + "_ref.lst"
    if not os.path.exists(refined_file):
        errors.append("-o: no %s" % os.path.basename(refined_file))
    else:
        count = sum(len(read_panels(f)) for f in panel_files(refined_file))
        if not refined or count != int(refined[-1]):
            errors.append("-o: %d panels in %s, %s after refinement"
                          % (count, os.path.basename(refined_file), refined[-1] if refined else "?"))

    matrix_file = base + ".csv"
    if not os.path.exists(matrix_file):
        errors.append("-e: no %s" % os.path.basename(matrix_file))
    else:
        with open(matrix_file) as f:
            written = [[float(v) for v in line.split(",")] for line in f if line.strip()]
        if len(written) != n or any(len(r) != n for r in written):
            errors.append("-e: %s is not a %d x %d matrix" % (os.path.basename(matrix_file), n, n))
        else:
            for i in range(n):
                for j in range(n):
                    if relative_difference(written[i][j], rows[i][j]) > OUTPUT_TOLERANCE:
                        errors.append("-e: C%d%d = %g in the file, %g printed"
                                      % (i + 1, j + 1, written[i][j], rows[i][j]))

    # FasterCap prints one excitation per row: with conductor i at 1 V, the charge on conductor
    # j is Cij; <base>_<conductor i>.lst names the panel files of all conductors in their order
    for i, name in enumerate(names):
        charge_file = "%s_%s.lst" % (base, name)
        if not os.path.exists(charge_file):
            errors.append("-c: no %s" % os.path.basename(charge_file))
            continue
        files = panel_files(charge_file)
        if len(files) != n:
            errors.append("-c: %s names %d conductors" % (os.path.basename(charge_file), len(files)))
            continue
        for j, path in enumerate(files):
            charge = sum(area(vertices) * density for vertices, density in read_panels(path))
            if abs(charge - rows[i][j]) > OUTPUT_TOLERANCE * rows[i][i]:
                errors.append("-c: charge %g on conductor %d with %d at 1 V, C%d%d = %g"
                              % (charge, j + 1, i + 1, i + 1, j + 1, rows[i][j]))
    return errors


def out_of_core_ratio(fastercap, list_file, options):
    """-f value for blocks of BLOCK_MBYTES, from the memory a probe run reports"""
    command, _, output = run(fastercap, list_file, options + ["-f1e12"])
    links = re.search(r"links is: (\d+) Mbytes", output)
    free = re.search(r"Available free memory left is: (\d+) Mbytes", output)
    if not links or not free:
        return None, "%s reported no link or free memory" % " ".join(command)
    if int(links.group(1)) <= CHUNK_MBYTES:
        return None, "the links need %s Mbytes, not more than one chunk" % links.group(1)
    return int(free.group(1)) / float(BLOCK_MBYTES), None


def interrupt_out_of_core(fastercap, list_file, options, ratio, workdir):
    """errors and output of an out-of-core run interrupted after writing two link chunks"""
    tmpdir = os.path.join(workdir, "tmp")
    os.makedirs(tmpdir, exist_ok=True)

    def temporary_files():
        return [f for f in os.listdir(tmpdir) if TEMPORARY_FILE.match(f)]

    for f in temporary_files():
        os.remove(os.path.join(tmpdir, f))
    command = [os.path.abspath(fastercap), os.path.basename(list_file), "-b"] + options + ["-f%g" % ratio]
    log = os.path.join(workdir, "interrupted.log")
    with open(log, "w") as out:
        process = subprocess.Popen(command, cwd=os.path.dirname(list_file), stdout=out,
                                   stderr=subprocess.STDOUT, env=dict(os.environ, TMPDIR=tmpdir))
        deadline = time.time() + RUN_TIMEOUT
        while process.poll() is None and time.time() < deadline \
                and sum(f.startswith("frcl") for f in temporary_files()) < 2:
            time.sleep(0.01)
        errors = []
        if process.poll() is not None:
            errors.append("the run ended before writing two link chunks")
        else:
            process.send_signal(signal.SIGINT)
        try:
            returncode = process.wait(timeout=RUN_TIMEOUT)
        except subprocess.TimeoutExpired:
            process.kill()
            returncode = process.wait()
            errors.append("no exit %d s after SIGINT" % RUN_TIMEOUT)
    print(" ".join(command))
    if returncode != USER_BREAK_EXIT_CODE:
        errors.append("exit code %d, expected %d" % (returncode, USER_BREAK_EXIT_CODE))
    left = temporary_files()
    if left:
        errors.append("%d temporary files left, e.g. %s" % (len(left), left[0]))
    with open(log) as f:
        return errors, f.read()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--fastercap", required=True, help="FasterCap executable")
    parser.add_argument("--geometry", required=True, choices=sorted(geometries.GEOMETRIES))
    parser.add_argument("--tolerance", required=True, type=float, help="relative tolerance")
    parser.add_argument("--workdir", required=True, help="directory for the input files")
    parser.add_argument("--out-of-core", action="store_true",
                        help="run again out-of-core and compare the matrices")
    parser.add_argument("--interrupt", action="store_true",
                        help="interrupt an out-of-core run instead (not on Windows)")
    parser.add_argument("--exit-code", type=int, default=0, help="expected exit code")
    parser.add_argument("--no-matrix", action="store_true",
                        help="FasterCap must stop without a capacitance matrix")
    parser.add_argument("--reference", choices=sorted(geometries.GEOMETRIES),
                        help="geometry of the same panels, which must give the same matrix")
    parser.add_argument("--names", help="the conductor names, comma-separated")
    parser.add_argument("--same-as", help="other FasterCap options that must give the same matrix")
    parser.add_argument("--dump-input", action="store_true",
                        help="dump the input geometry (-oi), which must give the same matrix")
    parser.add_argument("--output-files", action="store_true",
                        help="write refined geometry, matrix and charges (-o -e -c) and check them")
    parser.add_argument("fastercap_args", nargs=argparse.REMAINDER,
                        help="FasterCap options after --, e.g. -- -a0.01")
    args = parser.parse_args()
    fastercap_args = args.fastercap_args[1:] if args.fastercap_args[:1] == ["--"] else args.fastercap_args

    geometry = geometries.GEOMETRIES[args.geometry]
    os.makedirs(args.workdir, exist_ok=True)
    list_file = os.path.abspath(geometry.write(args.workdir))

    if args.interrupt:
        ratio, error = out_of_core_ratio(args.fastercap, list_file, fastercap_args)
        if error:
            errors, output = [error], ""
        else:
            errors, output = interrupt_out_of_core(args.fastercap, list_file, fastercap_args, ratio,
                                                   args.workdir)
        return report(errors, output)

    output_options = ["-o", "-e", "-c"] if args.output_files else []
    command, returncode, output = run(args.fastercap, list_file, fastercap_args + output_options)
    rows, errors = check_run(command, returncode, output, geometry, args.tolerance,
                             args.exit_code, args.no_matrix)

    if args.names and not errors:
        names, _ = parse_last_matrix(output)
        if names != args.names.split(","):
            errors.append("conductors %s, expected %s" % (",".join(names), args.names))

    if args.same_as is not None and not errors:
        command, returncode, same_output = run(args.fastercap, list_file, args.same_as.split())
        same_rows, same_errors = check_run(command, returncode, same_output, geometry, args.tolerance)
        errors += ["%s: %s" % (args.same_as, e) for e in same_errors]
        if not same_errors:
            for i, (row, same_row) in enumerate(zip(rows, same_rows)):
                for j, (value, same_value) in enumerate(zip(row, same_row)):
                    if relative_difference(value, same_value) > REFERENCE_TOLERANCE:
                        errors.append("C%d%d = %s, %s gives %s" % (i + 1, j + 1, value, args.same_as, same_value))

    if args.output_files and not errors:
        names, _ = parse_last_matrix(output)
        errors += check_output_files(list_file, output, names, rows)

    if args.dump_input and not errors:
        command, returncode, dump_output = run(args.fastercap, list_file, fastercap_args + ["-oi"])
        dump_file = os.path.splitext(list_file)[0] + "_dump.lst"
        if returncode != 0 or parse_last_matrix(dump_output)[0] is not None:
            errors.append("-oi: exit code %s, expected 0 without a matrix" % returncode)
        elif not os.path.exists(dump_file):
            errors.append("-oi: no %s" % os.path.basename(dump_file))
        else:
            command, returncode, dump_output = run(args.fastercap, dump_file, fastercap_args)
            dump_rows, dump_errors = check_run(command, returncode, dump_output, geometry,
                                               args.tolerance)
            errors += ["dump: " + e for e in dump_errors]
            if not dump_errors:
                for i, (row, dump_row) in enumerate(zip(rows, dump_rows)):
                    for j, (value, dump_value) in enumerate(zip(row, dump_row)):
                        if relative_difference(value, dump_value) > REFERENCE_TOLERANCE:
                            errors.append("dump: C%d%d = %g, input gives %g"
                                          % (i + 1, j + 1, dump_value, value))

    if args.out_of_core and not errors:
        ratio, error = out_of_core_ratio(args.fastercap, list_file, fastercap_args)
        if error:
            errors.append(error)
        else:
            command, returncode, output = run(args.fastercap, list_file, fastercap_args + ["-f%g" % ratio])
            ooc_rows, ooc_errors = check_run(command, returncode, output, geometry, args.tolerance)
            errors += ["out-of-core: " + e for e in ooc_errors]
            percent = re.search(r"more than (\S+)% of the available free memory", output)
            if "Going out-of-core" not in output:
                errors.append("out-of-core: FasterCap stayed in-core")
            elif not percent or relative_difference(float(percent.group(1)), 100.0 / ratio) > 0.01:
                errors.append("out-of-core: logged links memory %s%% of the free memory, expected %.3g%%"
                              % (percent.group(1) if percent else "?", 100.0 / ratio))
            elif not ooc_errors:
                for i, (row, ooc_row) in enumerate(zip(rows, ooc_rows)):
                    for j, (value, ooc_value) in enumerate(zip(row, ooc_row)):
                        if relative_difference(value, ooc_value) > OUT_OF_CORE_TOLERANCE:
                            errors.append("out-of-core C%d%d = %g, in-core %g" % (i + 1, j + 1, ooc_value, value))

    if args.reference and not errors:
        reference = geometries.GEOMETRIES[args.reference]
        reference_file = os.path.abspath(reference.write(args.workdir))
        command, returncode, output = run(args.fastercap, reference_file, fastercap_args)
        reference_rows, reference_errors = check_run(command, returncode, output, reference,
                                                     args.tolerance)
        errors += ["reference: " + e for e in reference_errors]
        if not reference_errors:
            for i, (row, reference_row) in enumerate(zip(rows, reference_rows)):
                for j, (value, reference_value) in enumerate(zip(row, reference_row)):
                    if relative_difference(value, reference_value) > REFERENCE_TOLERANCE:
                        errors.append("C%d%d = %g, %s gives %g"
                                      % (i + 1, j + 1, value, args.reference, reference_value))

    return report(errors, output)


def report(errors, output):
    if errors:
        tail = output.splitlines()[-OUTPUT_TAIL_LINES:]
        print("FAILED:\n  " + "\n  ".join(errors))
        print("FasterCap output (last %d lines):\n%s" % (len(tail), "\n".join(tail)))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
