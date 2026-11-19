"""Physical properties of generated meshes, independent of the solver."""

import sys
import unittest

sys.dont_write_bytecode = True

import geometries


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


class CubeMesh(unittest.TestCase):
    def test_faces_point_outward_and_enclose_volume(self):
        for side in (1e-6, 1.0, 1e6):
            with self.subTest(side=side):
                volume = 0.0
                faces = geometries.cube_panels("cube", side)
                self.assertEqual(len(faces), 6)
                for _, _, vertices in faces:
                    a, b, c, d = vertices
                    normal = cross([b[k] - a[k] for k in range(3)], [d[k] - a[k] for k in range(3)])
                    outward = [sum(v[k] for v in vertices) / 4 - side / 2 for k in range(3)]
                    self.assertGreater(dot(normal, outward), 0, vertices)
                    volume += (dot(a, cross(b, c)) + dot(a, cross(c, d))) / 6
                self.assertAlmostEqual(volume / side ** 3, 1.0)


if __name__ == "__main__":
    unittest.main()
