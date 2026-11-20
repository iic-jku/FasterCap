"""Reject misleading platform comparisons without invoking the solver."""
import copy
import hashlib
import math
from pathlib import Path
import sys
import subprocess
import tempfile
import unittest

sys.dont_write_bytecode = True
import check_platform_comparison as comparison


def fixture():
    matrix = {"units": "F", "reference_column": None, "row_conductors": ["substrate", "plate"],
              "column_conductors": ["substrate", "plate"], "real": [[4.3e-7, -4.1e-7], [-4.1e-7, 4.2e-7]],
              "imaginary": [[0.0, 0.0], [0.0, 0.0]]}
    settings = {"automatic": False, "auto_tolerance": 0.05, "gmres_tolerance": 1e-10,
                "mesh_refinement": 1e32, "interaction_ratio": 1e-50, "curvature": 3.0,
                "scheme": "collocation", "preconditioner": 1}
    record = {"iteration": 0, "solve_status": {"exit_code": 0}, "finite": True,
              "panels": comparison.FIXED_PANELS, "links": comparison.FIXED_LINK_BOUND, "matrix": matrix, "settings": settings}
    report = {"schema_version": 1, "status": {"exit_code": 0}, "last_completed_iteration": 0,
              "iterations": [record], "build": {
                  "compiler": "test compiler", "configuration": "RelWithDebInfo", "target_system": "test",
                  "target_processor": "test", "wxwidgets": "test wx", "blas_link_libraries": "test blas",
                  "omp_max_threads": 1, "thread_environment": dict(comparison.THREADS),
                  "source_revisions": {key: str(i)*40 for i, key in enumerate(("FasterCap", "LinAlgebra", "Geometry"), 1)},
                  "source_dirty": {key: False for key in ("FasterCap", "LinAlgebra", "Geometry")}}}
    artifact = {"schema_version": 1, "geometry_sha256": {"fixture": "a"*64}, "cases": {}}
    for name in comparison.CASES:
        artifact["cases"][name] = {"returncode": 0, "report": copy.deepcopy(report),
                                  "geometry": "fixed_triangles" if name == "fixed_mesh" else "original_quads",
                                  "options": comparison.FIXED_OPTIONS if name == "fixed_mesh" else comparison.AUTO_OPTIONS}
    adaptive = artifact["cases"]["adaptive"]["report"]["iterations"][-1]
    adaptive["settings"].update(automatic=True, mesh_refinement=0.01, interaction_ratio=0.5, gmres_tolerance=0.025)
    replay = artifact["cases"]["manual_replay"]
    replay["options"] = comparison.replay_options(adaptive)
    replay["report"]["iterations"][-1]["settings"].update(mesh_refinement=0.01, interaction_ratio=0.5)
    return artifact


class PlatformComparison(unittest.TestCase):
    def compare(self, a, b):
        return comparison.compare([a, b], ["first", "second"])

    def test_roundtrip_differences_do_not_require_bitwise_identity(self):
        a, b = fixture(), fixture()
        b["cases"]["fixed_mesh"]["report"]["iterations"][0]["matrix"]["real"][0][0] += 1e-16
        result = self.compare(a, b)
        self.assertTrue(result["passed"])
        self.assertGreater(result["comparisons"][0]["difference"]["max_entry_relative"], 0)

    def test_original_percent_divergence_fails_even_with_different_paths(self):
        a, b = fixture(), fixture()
        record = b["cases"]["adaptive"]["report"]["iterations"][0]
        record["matrix"]["real"][0][1] *= 1.03
        record["panels"] += 1
        result = self.compare(a, b)
        self.assertFalse(result["passed"])
        item = next(v for v in result["comparisons"] if v["case"] == "adaptive")
        self.assertFalse(item["same_reported_path"])
        self.assertGreater(item["difference"]["max_entry_relative"], 0.02)

    def test_replay_budget_depends_on_reported_work(self):
        a, b = fixture(), fixture()
        record = b["cases"]["manual_replay"]["report"]["iterations"][0]
        record["matrix"]["real"][0][0] *= 1.0001
        self.assertFalse(self.compare(a, b)["passed"])
        record["panels"] += 1
        self.assertTrue(self.compare(a, b)["passed"])

    def test_replay_rejects_unrepresented_automatic_refinement(self):
        record = fixture()["cases"]["adaptive"]["report"]["iterations"][-1]
        record["settings"]["uniform_enrichment"] = 1
        with self.assertRaisesRegex(ValueError, "cannot represent"):
            comparison.replay_options(record)

    def test_fixed_mesh_requires_matching_work(self):
        b = fixture()
        b["cases"]["fixed_mesh"]["report"]["iterations"][0]["links"] -= 2
        self.assertFalse(self.compare(fixture(), b)["passed"])

    def test_timeout_preserves_both_output_streams(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            program = "import sys,time; print('stdout',flush=True); print('stderr',file=sys.stderr,flush=True); time.sleep(5)"
            with self.assertRaises(subprocess.TimeoutExpired):
                comparison.run_solver([sys.executable, "-c", program], root, root/"report.json", root/"solver.log", 1.0)
            output = (root/"solver.log").read_text()
            self.assertIn("stdout", output)
            self.assertIn("stderr", output)

    def test_oversized_log_is_rejected_and_retention_is_capped(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            program = "import sys; sys.stdout.write('x' * %d)" % (comparison.MAX_OUTPUT_BYTES + 1024)
            with self.assertRaisesRegex(ValueError, "exceeded 2 MiB"):
                comparison.run_solver([sys.executable, "-c", program], root, root/"report.json", root/"solver.log", 2)
            self.assertEqual((root/"solver.log").stat().st_size, comparison.MAX_OUTPUT_BYTES)

    def test_every_pair_is_compared(self):
        result = comparison.compare([fixture(), fixture(), fixture()], ["a", "b", "c"])
        self.assertEqual(len(result["comparisons"]), 9)

    def test_incomplete_sources_geometry_options_and_status_are_rejected(self):
        for mutation in (
            lambda b: b["cases"].pop("manual_replay"),
            lambda b: b.update(error="timeout"),
            lambda b: b["geometry_sha256"].update(fixture="different"),
            lambda b: b["cases"]["fixed_mesh"].update(options=[]),
            lambda b: b["cases"]["fixed_mesh"].update(geometry="original_quads"),
            lambda b: b["cases"]["fixed_mesh"]["report"]["build"]["source_revisions"].update(Geometry="other"),
            lambda b: b["cases"]["fixed_mesh"]["report"]["build"]["source_dirty"].update(FasterCap=True),
            lambda b: b["cases"]["fixed_mesh"]["report"]["build"].update(compiler=None),
            lambda b: b["cases"]["fixed_mesh"]["report"]["build"].update(omp_max_threads=4),
            lambda b: b["cases"]["fixed_mesh"]["report"]["iterations"][0]["settings"].update(gmres_tolerance=0.025),
            lambda b: b["cases"]["fixed_mesh"]["report"]["status"].update(exit_code=75),
            lambda b: b["cases"]["fixed_mesh"]["report"]["iterations"][0].update(finite=False),
            lambda b: b["cases"]["fixed_mesh"]["report"]["iterations"][0].update(links=comparison.FIXED_LINK_BOUND+1),
        ):
            with self.subTest(mutation=mutation):
                b = fixture()
                mutation(b)
                with self.assertRaises(ValueError):
                    self.compare(fixture(), b)

    def test_nonfinite_imaginary_and_swapped_names_are_rejected(self):
        for value in (float("nan"), float("inf"), None, True):
            b = fixture()
            b["cases"]["fixed_mesh"]["report"]["iterations"][0]["matrix"]["imaginary"][0][1] = value
            with self.assertRaises(ValueError):
                self.compare(fixture(), b)
        b = fixture()
        b["cases"]["fixed_mesh"]["report"]["iterations"][0]["matrix"]["column_conductors"].reverse()
        with self.assertRaises(ValueError):
            self.compare(fixture(), b)

    def test_complex_difference_is_not_ignored(self):
        b = fixture()
        b["cases"]["fixed_mesh"]["report"]["iterations"][0]["matrix"]["imaginary"][0][1] = 1e-9
        self.assertFalse(self.compare(fixture(), b)["passed"])

    def test_explicit_fixed_triangles_have_bounded_input_size(self):
        with tempfile.TemporaryDirectory() as temporary:
            comparison.geometries.write_plate_over_substrate(temporary, triangulated=True)
            panels = [line.split() for p in Path(temporary).glob("*.txt")
                      for line in p.read_text().splitlines() if not line.startswith("0 ")]
        self.assertEqual(len(panels), comparison.FIXED_PANELS)
        self.assertTrue(all(row[0] == "T" and len(row) == 11 for row in panels))

    def test_fixed_triangles_cover_the_original_oriented_faces(self):
        def read(directory):
            result = {}
            for path in directory.glob("*.txt"):
                for line in path.read_text().splitlines()[1:]:
                    fields = line.split()
                    values = list(map(float, fields[2:]))
                    points = [values[i:i+3] for i in range(0, len(values), 3)]
                    plane = next((axis, points[0][axis]) for axis in range(3)
                                 if all(p[axis] == points[0][axis] for p in points))
                    key = (fields[1], plane)
                    a = [points[1][i] - points[0][i] for i in range(3)]
                    b = [points[-1][i] - points[0][i] for i in range(3)]
                    normal = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
                    factor = 0.5 if fields[0] == "T" else 1.0
                    result.setdefault(key, []).append([factor*x for x in normal])
                    if fields[0] == "T":
                        edges = [math.dist(points[i], points[(i+1)%3]) for i in range(3)]
                        self.assertGreater(min(edges), 0)
                        self.assertLessEqual(max(edges)/min(edges), math.sqrt(101)*(1+1e-8))
            return {key: [math.fsum(v[axis] for v in values) for axis in range(3)]
                    for key, values in result.items()}
        with tempfile.TemporaryDirectory() as temporary:
            quad, triangle = Path(temporary)/"quad", Path(temporary)/"triangle"
            quad.mkdir(); triangle.mkdir()
            comparison.geometries.write_plate_over_substrate(str(quad))
            comparison.geometries.write_plate_over_substrate(str(triangle), triangulated=True)
            original, fixed = read(quad), read(triangle)
        self.assertEqual(set(original), set(fixed))
        self.assertEqual(len(original), 18)
        for face in original:
            for expected, actual in zip(original[face], fixed[face]):
                self.assertTrue(math.isclose(expected, actual, rel_tol=1e-9, abs_tol=1e-12))

    def test_fixture_preserves_original_panel_records(self):
        # Digests from the generator in 9d7d5b8 (#23); protects historical face
        # ordering, oxide side/reference records, coordinates and units.
        expected = {
            "plate_over_substrate.lst": "40dfa9268add5d0435cb266913ea114d9f2b390026de996f8e30570234e12f6d",
            "substrate.txt": "a110b1ff6669337376abfef70b33d18b1e48ae409bd9407327dbaebbea943329",
            "plate.txt": "e98ecfbbdc37b08875c125783815ed52bc148f259bd8c75210235e33d103c223",
            "oxide.txt": "1298126a3acbf1e1a64232d4376e098758f849c5bd3286ae57cb47ea0581c4ce",
        }
        with tempfile.TemporaryDirectory() as temporary:
            comparison.geometries.write_plate_over_substrate(temporary)
            actual = {p.name: hashlib.sha256(p.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                      for p in Path(temporary).iterdir()}
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
