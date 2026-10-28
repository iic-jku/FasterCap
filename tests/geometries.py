"""Test geometries in the FasterCap 3D input format, with what is known about their matrices.

Each geometry writes its panel files and list file into a directory. Lengths are in
meters, so the matrices are in farads.
"""

import collections
import math
import os

EPS0 = 8.8541878128e-12  # F/m
FOUR_PI_EPS0 = 4.0 * math.pi * EPS0

# Capacitance of a cube of side a in free space is CUBE_FACTOR * 4*pi*eps0 * a
# (Hwang and Mascagni, J. Appl. Phys. 95, 3798 (2004): 0.66067813(1))
CUBE_FACTOR = 0.66067813

# write(directory) -> list file path; analytic: {(i, j): farads}, 1-based;
# equal_entries: pairs of (i, j) that are equal by the symmetry of the geometry
Geometry = collections.namedtuple("Geometry", "write conductors analytic equal_entries")


def _fmt(values):
    return " ".join("%.12g" % v for v in values)


def _write_panels(path, title, panels):
    with open(path, "w") as f:
        f.write("0 %s\n" % title)
        for kind, name, vertices in panels:
            f.write("%s %s %s\n" % (kind, name, _fmt(c for v in vertices for c in v)))


def _write_list(path, lines):
    with open(path, "w") as f:
        f.write("".join(line + "\n" for line in lines))


def cube_panels(name, side):
    def v(x, y, z):
        return (x * side, y * side, z * side)

    faces = [
        [v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0)],
        [v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)],
        [v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)],
        [v(0, 1, 0), v(1, 1, 0), v(1, 1, 1), v(0, 1, 1)],
        [v(0, 0, 0), v(0, 1, 0), v(0, 1, 1), v(0, 0, 1)],
        [v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)],
    ]
    return [("Q", name, face) for face in faces]


def icosphere_panels(name, radius, depth):
    """Triangles of an icosahedron subdivided `depth` times, vertices on the sphere."""
    t = (1.0 + math.sqrt(5.0)) / 2.0
    vertices = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
                (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
                (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]

    def project(p):
        n = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
        return (radius * p[0] / n, radius * p[1] / n, radius * p[2] / n)

    triangles = [tuple(project(vertices[i]) for i in face) for face in faces]
    for _ in range(depth):
        finer = []
        for a, b, c in triangles:
            ab = project([(a[k] + b[k]) / 2 for k in range(3)])
            bc = project([(b[k] + c[k]) / 2 for k in range(3)])
            ca = project([(c[k] + a[k]) / 2 for k in range(3)])
            finer += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        triangles = finer
    return [("T", name, tri) for tri in triangles]


SIDE = 1.0
GAP = 0.5
RADIUS = 1.0
COAT_RADIUS = 2.0
COAT_PERMITTIVITY = 4.0


def write_cube(directory):
    _write_panels(os.path.join(directory, "cube.txt"), "cube", cube_panels("c", SIDE))
    path = os.path.join(directory, "cube.lst")
    _write_list(path, ["C cube.txt 1.0 0 0 0"])
    return path


def write_three_cubes(directory):
    """Cube 1 at the origin, cubes 2 and 3 beside it along x and along y."""
    _write_panels(os.path.join(directory, "cube.txt"), "cube", cube_panels("c", SIDE))
    pitch = SIDE + GAP
    path = os.path.join(directory, "three_cubes.lst")
    _write_list(path, ["C cube.txt 1.0 0 0 0",
                       "C cube.txt 1.0 %g 0 0" % pitch,
                       "C cube.txt 1.0 0 %g 0" % pitch])
    return path


def write_sphere(directory):
    """1280 triangles inscribed in the sphere: about 0.3 % below its capacitance."""
    _write_panels(os.path.join(directory, "sphere.txt"), "sphere",
                  icosphere_panels("s", RADIUS, 3))
    path = os.path.join(directory, "sphere.lst")
    _write_list(path, ["C sphere.txt 1.0 0 0 0"])
    return path


def write_coated_sphere(directory):
    """The sphere in a concentric dielectric shell, both of 1280 triangles."""
    _write_panels(os.path.join(directory, "sphere.txt"), "sphere",
                  icosphere_panels("s", RADIUS, 3))
    _write_panels(os.path.join(directory, "coat.txt"), "coat",
                  icosphere_panels("d", COAT_RADIUS, 3))
    path = os.path.join(directory, "coated_sphere.lst")
    # vacuum outside the shell; the reference point, the center, is on the inner side ('-')
    _write_list(path, ["C sphere.txt %g 0 0 0" % COAT_PERMITTIVITY,
                       "D coat.txt 1.0 %g 0 0 0 0 0 0 -" % COAT_PERMITTIVITY])
    return path


# the shell, 4*pi*eps0*er / (1/a - 1/b), in series with the vacuum outside it, 4*pi*eps0*b
COATED_SPHERE_C = FOUR_PI_EPS0 / ((1.0 / RADIUS - 1.0 / COAT_RADIUS) / COAT_PERMITTIVITY
                                  + 1.0 / COAT_RADIUS)


def write_degenerate_triangle(directory):
    """The cube and, beside it, a triangle whose vertices are on a line: no self-potential."""
    panels = cube_panels("c", SIDE) + [("T", "c", [(2, 0, 0), (3, 0, 0), (4, 0, 0)])]
    _write_panels(os.path.join(directory, "degenerate.txt"), "cube and degenerate triangle", panels)
    path = os.path.join(directory, "degenerate_triangle.lst")
    _write_list(path, ["C degenerate.txt 1.0 0 0 0"])
    return path


def write_patran_dielectric(directory):
    """The cube in a dielectric given as a PATRAN neutral file, which FasterCap can't read."""
    _write_panels(os.path.join(directory, "cube.txt"), "cube", cube_panels("c", SIDE))
    # title and end packets, as in FastCap2's examples/*.neu
    with open(os.path.join(directory, "coat.neu"), "w") as f:
        f.write("25       0       0       1       0       0       0       0       0\n"
                "PATRAN DIELECTRIC\n"
                "99       0       0       1       0       0       0       0       0\n")
    path = os.path.join(directory, "patran_dielectric.lst")
    _write_list(path, ["C cube.txt 1.0 0 0 0",
                       "D coat.neu 1.0 2.0 0 0 0 0.5 0.5 0.5"])
    return path



# a quadrilateral whose corner (1, 1, 1) is off the plane of the other three
SKEWED_QUAD = [(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)]


def write_skewed_quad(directory):
    """A non-planar Q panel, split by FasterCap along its smaller diagonal."""
    _write_panels(os.path.join(directory, "skewed.txt"), "skewed quad", [("Q", "c", SKEWED_QUAD)])
    path = os.path.join(directory, "skewed_quad.lst")
    _write_list(path, ["C skewed.txt 1.0 0 0 0"])
    return path


def write_skewed_quad_triangles(directory):
    """The two triangles of the skewed quad, split along the diagonal from corner 1 to 3."""
    a, b, c, d = SKEWED_QUAD
    _write_panels(os.path.join(directory, "skewed_triangles.txt"), "skewed quad as triangles",
                  [("T", "c", [a, b, d]), ("T", "c", [b, c, d])])
    path = os.path.join(directory, "skewed_quad_triangles.lst")
    _write_list(path, ["C skewed_triangles.txt 1.0 0 0 0"])
    return path


GEOMETRIES = {
    "cube": Geometry(write_cube, 1, {(1, 1): CUBE_FACTOR * FOUR_PI_EPS0 * SIDE}, []),
    "three_cubes": Geometry(write_three_cubes, 3, {}, [((1, 2), (1, 3)), ((2, 2), (3, 3))]),
    "sphere": Geometry(write_sphere, 1, {(1, 1): FOUR_PI_EPS0 * RADIUS}, []),
    "coated_sphere": Geometry(write_coated_sphere, 1, {(1, 1): COATED_SPHERE_C}, []),
    "degenerate_triangle": Geometry(write_degenerate_triangle, 1, {}, []),
    "patran_dielectric": Geometry(write_patran_dielectric, 1, {}, []),
    "skewed_quad": Geometry(write_skewed_quad, 1, {}, []),
    "skewed_quad_triangles": Geometry(write_skewed_quad_triangles, 1, {}, []),
}
