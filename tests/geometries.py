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
# concentric spheres of #7 (ediloren/FasterCap#8)
INNER_RADIUS = 0.05
OUTER_RADIUS = 0.1


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


def write_concentric_spheres(directory):
    """Inner and outer sphere of 1280 triangles each, in one file.

    The triangles are 'Q' panels with a repeated vertex: only 'Q' panels go through
    C2DTriangulate, whose assert #7 reports, and -a reads the input at each iteration.
    """
    triangles = (icosphere_panels("inner", INNER_RADIUS, 3)
                 + icosphere_panels("outer", OUTER_RADIUS, 3))
    panels = [("Q", name, vertices + (vertices[0],)) for _, name, vertices in triangles]
    _write_panels(os.path.join(directory, "spheres.txt"), "concentric spheres", panels)
    path = os.path.join(directory, "concentric_spheres.lst")
    _write_list(path, ["C spheres.txt 1.0 0 0 0"])
    return path


def concentric_spheres_analytic(a, b):
    """Sphere of radius a inside a thin shell of radius b, in free space."""
    # all field lines of the inner sphere end on the shell: C12 = -C11
    c11 = FOUR_PI_EPS0 * a * b / (b - a)
    return {(1, 1): c11, (1, 2): -c11, (2, 1): -c11, (2, 2): c11 + FOUR_PI_EPS0 * b}


GEOMETRIES = {
    "cube": Geometry(write_cube, 1, {(1, 1): CUBE_FACTOR * FOUR_PI_EPS0 * SIDE}, []),
    "three_cubes": Geometry(write_three_cubes, 3, {}, [((1, 2), (1, 3)), ((2, 2), (3, 3))]),
    "sphere": Geometry(write_sphere, 1, {(1, 1): FOUR_PI_EPS0 * RADIUS}, []),
    "concentric_spheres": Geometry(write_concentric_spheres, 2,
                                   concentric_spheres_analytic(INNER_RADIUS, OUTER_RADIUS), []),
}
