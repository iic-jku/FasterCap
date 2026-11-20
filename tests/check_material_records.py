"""Check C/D record validation and parsed geometry using tiny input dumps only."""
import argparse
import fnmatch
from pathlib import Path
import re
import subprocess
import sys
import tempfile

sys.dont_write_bytecode = True
import geometries


def run(executable, root, dimensions, records, expected=0, field=None, line=2, child=None,
        filename="input.lst"):
    root.mkdir()
    geometries.write_material_record_case(str(root), dimensions, records)
    if child is not None:
        (root / "child.lst").write_text("* child\n" + "\n".join(child) + "\n")
    result = subprocess.run([executable, "-b", "-oi", "input.lst"], cwd=root,
                            capture_output=True, text=True, timeout=5)
    output = result.stdout + result.stderr
    if result.returncode != expected:
        raise AssertionError("exit %d, expected %d\n%s" % (result.returncode, expected, output))
    if "Capacitance matrix is:" in output or "assert" in output.lower():
        raise AssertionError("unexpected solve or assertion\n" + output)
    dumps = list(root.glob("input*dump*"))
    if expected:
        location = re.search(r'of input file "([^"]+)"', output)
        if not ("missing or invalid " + field in output and "line %d " % line in output
                and location and Path(location.group(1)).name == filename):
            raise AssertionError("missing field/file/line diagnostic\n" + output)
        if dumps or "cannot open file" in output.lower():
            raise AssertionError("invalid record reached its child or output stage\n" + output)
        return None
    if not (root / "input_dump.lst").exists():
        raise AssertionError("missing parsed geometry dump\n" + output)
    # All runs use the same basename. Compare material records, names and geometry;
    # discard comments identifying how a file was generated.
    return {path.name: [entry for entry in path.read_text().splitlines()
                        if entry and not entry.startswith(("*", "0 "))]
            for path in dumps}


def failures():
    for dimensions in (3, 2):
        zero = " ".join(["0"] * dimensions)
        for kind in ("C", "D"):
            prefix = kind + " panels.txt " + ("1" if kind == "C" else "1 2")
            for count in range(dimensions):
                # The first case is the original missing-offset acceptance defect.
                record = prefix + " " + " ".join(["0"] * count)
                yield "missing-%dd-%s-%d" % (dimensions, kind.lower(), count), dimensions, [record], "translation coordinates", 2
            for field, before in [("translation coordinates", prefix)] + (
                    [("dielectric reference coordinates", prefix + " " + zero)] if kind == "D" else []):
                if field.startswith("dielectric"):
                    for count in range(dimensions):
                        record = before + " " + " ".join(["0"] * count)
                        yield "reference-%dd-%d" % (dimensions, count), dimensions, [record], field, 2
                for value in ("text", "nan", "inf", "-inf", "1e100", "0junk"):
                    record = before + " " + " ".join(["0"] * (dimensions - 1) + [value])
                    if field == "translation coordinates" and kind == "D":
                        record += " " + zero
                    yield "%dd-%s-%s-%s" % (dimensions, kind, field.split()[0], value), dimensions, [record], field, 2
            # A preceding valid record must not provide the missing record's values.
            previous = "C panels.txt 1 " + " ".join(["2"] * dimensions)
            yield "stale-%dd-%s" % (dimensions, kind), dimensions, [previous, prefix], "translation coordinates", 3
            for prefix, field in [(kind, "file name"),
                                  (kind + " absent.txt", "outer permittivity"),
                                  (kind + " absent.txt bad", "outer permittivity"),
                                  (kind + " absent.txt 1-j", "outer permittivity"),
                                  (kind + " absent.txt 1-jinf", "outer permittivity")]:
                yield "header-%dd-%s-%s" % (dimensions, kind, prefix), dimensions, [prefix], field, 2
        for value in ("", "bad", "nan", "2+j", "2-jinf"):
            yield "inner-%dd-%s" % (dimensions, value), dimensions, ["D absent.txt 1 " + value], "inner permittivity", 2
        for value in ("nan", "inf", "1e100"):
            for kind in ("C", "D"):
                record = kind + " absent.txt " + value + (" 2 " if kind == "D" else " ") + zero
                if kind == "D":
                    record += " " + zero
                yield "nonfinite-%dd-%s-%s" % (dimensions, kind, value), dimensions, [record], "outer permittivity", 2
            record = "D absent.txt 1 " + value + " " + zero + " " + zero
            yield "nonfinite-inner-%dd-%s" % (dimensions, value), dimensions, [record], "inner permittivity", 2
        yield "dummy-reference-%dd" % dimensions, dimensions, ["D absent.txt 1 -2 " + zero], "dielectric reference coordinates", 2
        # The last field must be validated before attempting to open the child.
        yield "before-open-%dd" % dimensions, dimensions, ["C absent.txt 1 " + zero[:-1]], "translation coordinates", 2


def check_valid(executable, root, dimensions):
    zero = " ".join(["0"] * dimensions)
    shift = "-1 2 +3" if dimensions == 3 else "-1 +2"
    reference = "4 5 6" if dimensions == 3 else "4 5"
    records = ["C panels.txt 3-j0.02 " + shift + " + * trailing comment",
               "C panels.txt 3-j0.02 " + zero,
               "D panels.txt 4+j0.5 2-j0.25 " + zero + " " + reference + " -"]
    equivalent = ["C\tpanels.txt\t3 - j 0.02\t" + shift + " +\r",
                  "C panels.txt 3 - j 0.02 " + zero + " * trailing comment",
                  "D panels.txt 2 - j 0.25 4 + j 0.5 " + zero + " " + reference + "\t\r"]
    first = run(executable, root / "compact", dimensions, records)
    second = run(executable, root / "spaced", dimensions, equivalent)
    if first != second:
        raise AssertionError("complex spelling or D '-' changes parsed materials/geometry")
    attached = [records[0].replace(" + *", "+ *"), records[1], records[2][:-2] + "-"]
    if first != run(executable, root / "attached", dimensions, attached):
        raise AssertionError("attached optional marker changes parsed materials/geometry")
    materials = first["input_dump.lst"]
    if sum(entry.startswith("C ") for entry in materials) != 1:
        raise AssertionError("C '+' failed to merge the two conductor records")
    if sum(entry.startswith("D ") for entry in materials) != 1:
        raise AssertionError("valid dielectric record was skipped")

    # Real permittivity followed by a signed coordinate must remain real; signed
    # inner permittivity retains the existing dummy-interface skip behavior.
    conductor = "C panels.txt 1 " + shift
    plain = run(executable, root / "real", dimensions, [conductor])
    dummy = run(executable, root / "dummy", dimensions,
                [conductor, "D panels.txt 1 -2 " + zero + " " + reference])
    if plain != dummy:
        raise AssertionError("negative inner permittivity changed dummy-interface behavior")
    exponent = shift.replace("-1", "-1e0").replace("2", "2e0").replace("+3", "+3e0")
    if plain != run(executable, root / "exponents", dimensions, ["C panels.txt 1 " + exponent]):
        raise AssertionError("signed/exponent coordinates changed geometry")
    zero_dump = run(executable, root / "zero", dimensions, ["C panels.txt 1 " + zero])
    underflow = " ".join(["1e-50"] * dimensions)
    if zero_dump != run(executable, root / "underflow", dimensions, ["C panels.txt 1 " + underflow]):
        raise AssertionError("finite float underflow changed coordinate handling")

    # A parent translates the reference; the D record's local translation does not.
    child = ["C panels.txt 1 " + zero,
             "D panels.txt 1 2 " + " ".join(["1"] * dimensions) + " " + reference]
    external = run(executable, root / "external", dimensions, ["C child.lst 1 " + shift], child=child)
    embedded = run(executable, root / "embedded", dimensions,
                   ["C child.lst 1 " + shift, "F child.lst"] + child + ["End"])
    if external != embedded:
        raise AssertionError("embedded and external child records differ")
    dielectric = next(entry.split() for entry in external["input_dump.lst"] if entry.startswith("D "))
    expected_reference = [float(a) + float(b) for a, b in zip(shift.split(), reference.split())]
    if [float(value) for value in dielectric[7:7 + dimensions]] != expected_reference:
        raise AssertionError("dielectric reference was not translated by the parent only")


def check_nested_invalid(executable, root, dimensions):
    zero = " ".join(["0"] * dimensions)
    child = ["C panels.txt 1 " + zero,
             "D absent.txt 1 2 " + zero + " " + zero[:-1]]
    run(executable, root / "external", dimensions, ["C child.lst 1 " + zero],
        74, "dielectric reference coordinates", 3, child, "child.lst")
    run(executable, root / "embedded", dimensions,
        ["C child.lst 1 " + zero, "F child.lst"] + child + ["End"],
        74, "dielectric reference coordinates", 5, filename="child.lst")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable")
    parser.add_argument("--case", default="*", help="select case names, e.g. missing-3d-c-0")
    args = parser.parse_args()
    executable = str(Path(args.executable).resolve())
    count = 0
    with tempfile.TemporaryDirectory(prefix="fc-material-records-") as temporary:
        root = Path(temporary)
        for name, dimensions, records, field, line in failures():
            if fnmatch.fnmatchcase(name, args.case):
                try:
                    run(executable, root / str(count), dimensions, records, 74, field, line)
                except AssertionError as error:
                    raise AssertionError(name + ": " + str(error)) from error
                count += 1
        for dimensions in (2, 3):
            for prefix, check in (("valid", check_valid), ("nested", check_nested_invalid)):
                name = "%s-%dd" % (prefix, dimensions)
                if fnmatch.fnmatchcase(name, args.case):
                    directory = root / name
                    directory.mkdir()
                    check(executable, directory, dimensions)
                    count += 1
    if not count:
        raise AssertionError("case selector matched no tests")
    print("%d material-record validation cases passed" % count)


if __name__ == "__main__":
    main()
