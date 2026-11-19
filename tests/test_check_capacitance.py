"""Regression tests of matrix validation, including the paths that report failures."""

import contextlib
import io
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.dont_write_bytecode = True

import check_capacitance as checker
import geometries


class MatrixChecks(unittest.TestCase):
    def geometry(self, conductors=2, analytic=None, equal_entries=None):
        return geometries.Geometry(None, conductors, analytic or {}, equal_entries or [])

    def test_complex_reciprocity_uses_transpose(self):
        rows = [[2-0.2j, -1+0.1j], [-1+0.1j, 2-0.2j]]
        self.assertEqual(checker.check_matrix(rows, self.geometry(), 0.01), [])

    def test_complex_failures_are_reported(self):
        cases = [
            ([[-2+0.2j]], self.geometry(1), "not positive"),
            ([[2-0.2j, 1-0.1j], [1-0.1j, 2-0.2j]], self.geometry(), "not negative"),
            ([[2-0.2j, -1+0.1j], [-0.9+0.09j, 2-0.2j]], self.geometry(), "differ"),
            ([[2-0.2j, -1+0.1j], [-1+0.1j, 3-0.3j]],
             self.geometry(equal_entries=[((1, 1), (2, 2))]), "should be equal"),
        ]
        for rows, geometry, diagnostic in cases:
            with self.subTest(diagnostic=diagnostic):
                errors = checker.check_matrix(rows, geometry, 0.01)
                self.assertTrue(any(diagnostic in error for error in errors), errors)

    def test_nonfinite_entries_fail_before_physical_checks(self):
        for value in (math.nan, math.inf, -math.inf, complex(2, math.nan), complex(2, math.inf)):
            with self.subTest(value=value):
                errors = checker.check_matrix([[value]], self.geometry(1), 0.01)
                self.assertEqual(len(errors), 1)
                self.assertIn("C11 is", errors[0])

    def test_nonfinite_comparisons_fail(self):
        for value in (math.nan, math.inf, complex(2, math.nan), complex(2, math.inf)):
            with self.subTest(value=value):
                self.assertEqual(checker.relative_difference(value, 2), math.inf)
                self.assertEqual(checker.relative_difference(2, value), math.inf)
        self.assertEqual(checker.relative_difference(0, 0), 0)

    def test_analytic_tolerance_does_not_relax_reciprocity(self):
        rows = [[2.1, -1], [-0.98, 2.1]]
        geometry = self.geometry(analytic={(1, 1): 2.0})
        errors = checker.check_matrix(rows, geometry, 0.01, analytic_tolerance=0.1)
        self.assertEqual(len(errors), 1)
        self.assertIn("C12", errors[0])
        self.assertIn("differ", errors[0])

    def test_reciprocity_tolerance_does_not_relax_analytic(self):
        rows = [[2.1, -1], [-0.98, 2.1]]
        geometry = self.geometry(analytic={(1, 1): 2.0})
        errors = checker.check_matrix(rows, geometry, 0.01, reciprocity_tolerance=0.1)
        self.assertEqual(len(errors), 1)
        self.assertIn("analytic", errors[0])

    def test_row_sum_tolerance_is_independent(self):
        rows = [[1, -1.005], [-1.005, 2]]
        self.assertEqual(checker.check_matrix(rows, self.geometry(), 0.01), [])
        errors = checker.check_matrix(rows, self.geometry(), 0.01, row_sum_tolerance=0.001)
        self.assertEqual(len(errors), 1)
        self.assertIn("row 1 sums", errors[0])

    def test_equality_tolerance_is_independent(self):
        rows = [[2, -1], [-1, 2.01]]
        geometry = self.geometry(equal_entries=[((1, 1), (2, 2))])
        self.assertEqual(checker.check_matrix(rows, geometry, 0.01), [])
        errors = checker.check_matrix(rows, geometry, 0.01, equality_tolerance=0.001)
        self.assertEqual(len(errors), 1)
        self.assertIn("should be equal", errors[0])

    def test_zero_tolerance_is_not_replaced_by_default(self):
        rows = [[2, -1], [-1.000001, 2]]
        self.assertEqual(checker.check_matrix(rows, self.geometry(), 0.01), [])
        self.assertTrue(checker.check_matrix(rows, self.geometry(), 0.01, reciprocity_tolerance=0))

    def test_tolerances_must_be_finite_and_nonnegative(self):
        for text in ("nan", "inf", "-inf", "-0.1"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                checker.tolerance_value(text)
        self.assertEqual(checker.tolerance_value("0"), 0)
        self.assertEqual(checker.tolerance_value("1e-3"), 0.001)


class MatrixOutput(unittest.TestCase):
    def check_run(self, output, **kwargs):
        geometry = geometries.Geometry(None, 2, {}, [])
        with contextlib.redirect_stdout(io.StringIO()):
            return checker.check_run(["FasterCap"], 0, output, geometry, 0.01, **kwargs)

    def test_last_complex_matrix_is_parsed(self):
        output = ("Capacitance matrix is:\nDimension 1 x 1\nold 1\n"
                  "Capacitance matrix is:\nDimension 2 x 2\n"
                  "a 2-0.2j -1+0.1j\nb -1+0.1j 2-0.2j\nDone\n")
        names, rows = checker.parse_last_matrix(output)
        self.assertEqual(names, ["a", "b"])
        self.assertEqual(rows, [[2-0.2j, -1+0.1j], [-1+0.1j, 2-0.2j]])
        self.assertEqual(self.check_run(output)[1], [])

    def test_malformed_matrices_fail_without_traceback(self):
        for body in ("", "Dimension 2 x 3\na 2 -1\nb -1 2\n",
                     "Dimension 2 x 2\na 2 -1\n", "Dimension 2 x 2\na 2 -1\n\n",
                     "Dimension 2 x 2\na 2 -1\nb -1 broken\n",
                     "Dimension 2 x 2\na 2\nb -1 2\n", "Dimension 0 x 0\n"):
            with self.subTest(body=body):
                rows, errors = self.check_run("Capacitance matrix is:\n" + body)
                self.assertTrue(errors)
                self.assertIn("matrix", errors[0])
                self.assertEqual(rows, [])

    def test_malformed_matrix_does_not_pass_no_matrix(self):
        _, errors = self.check_run("Capacitance matrix is:\ninvalid\n", no_matrix=True)
        self.assertTrue(errors)

    def test_csv_accepts_complex_and_rejects_invalid_values(self):
        rows = [[2-0.2j, -1+0.1j], [-1+0.1j, 2-0.2j]]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "matrix.csv"
            path.write_text("2-0.2j, -1+0.1j\n-1+0.1j, 2-0.2j\n")
            self.assertEqual(checker.check_matrix_file(path, rows), [])
            for entry in ("nan", "inf", "2+nanj", "2-0.3j", "bad"):
                with self.subTest(entry=entry):
                    path.write_text(entry + ", -1+0.1j\n-1+0.1j, 2-0.2j\n")
                    errors = checker.check_matrix_file(path, rows)
                    self.assertTrue(errors)
                    self.assertTrue(all(error.startswith("-e:") for error in errors))
            path.write_text("2-0.2j\n")
            self.assertTrue(checker.check_matrix_file(path, rows))

    def test_nonfinite_charge_file_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / "model"
            panel = "Q c 0 0 0 1 0 0 1 1 0 0 1 0"
            (base.with_suffix(".csv")).write_text("1\n")
            (Path(directory) / "model_ref.lst").write_text("C panel.txt 1 0 0 0\n")
            (Path(directory) / "panel.txt").write_text(panel + "\n")
            (Path(directory) / "model_a.lst").write_text("C charge.txt 1 0 0 0\n")
            charge_path = Path(directory) / "charge.txt"
            output = "Number of panels after refinement: 1\n"
            charge_path.write_text(panel + " 1\n")
            self.assertEqual(checker.check_output_files(str(base) + ".lst", output, ["a"], [[1.0]]), [])
            for density in ("nan", "inf", "-inf"):
                with self.subTest(density=density):
                    charge_path.write_text(panel + " " + density + "\n")
                    errors = checker.check_output_files(str(base) + ".lst", output, ["a"], [[1.0]])
                    self.assertEqual(len(errors), 1)
                    self.assertIn("-c: charge", errors[0])

    def test_cli_overrides_reach_reference_runs(self):
        output = ("Capacitance matrix is:\nDimension 3 x 3\n"
                  "a 3 -1 -1\nb -1 3 -1\nc -1 -1 3\n")
        overrides = dict(analytic_tolerance=0.002, reciprocity_tolerance=0.003,
                         row_sum_tolerance=0.004, equality_tolerance=0.005)
        with tempfile.TemporaryDirectory() as directory:
            argv = ["check_capacitance.py", "--fastercap", "unused", "--geometry", "three_cubes",
                    "--tolerance", "0.01", "--workdir", directory, "--same-as=-a0.01",
                    "--reference", "three_cubes"]
            for name, value in overrides.items():
                argv += ["--" + name.replace("_", "-"), str(value)]
            with mock.patch.object(sys, "argv", argv), \
                    mock.patch.object(checker, "run", return_value=(["FasterCap"], 0, output)), \
                    mock.patch.object(checker, "check_matrix", wraps=checker.check_matrix) as check, \
                    contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(checker.main(), 0)
            self.assertEqual(check.call_count, 3)
            for call in check.call_args_list:
                self.assertEqual(call.kwargs, overrides)


if __name__ == "__main__":
    unittest.main()
