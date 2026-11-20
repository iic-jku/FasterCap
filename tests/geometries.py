"""Test geometries in the FasterCap input formats, with what is known about their matrices.

Each geometry writes its panel files and list file into a directory. Lengths are in
meters, so the matrices are in farads, and in farads per meter for the 2D ones.
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


def write_material_record_case(directory, dimensions, records):
    """One triangle or segment for C/D record parsing, without a numerical solve."""
    panel = ("T", "plate", [(0, 0, 0), (1, 0, 0), (0, 1, 0)]) if dimensions == 3 else (
        "S", "plate", [(0, 0), (1, 0)])
    _write_panels(os.path.join(directory, "panels.txt"), "material record panel", [panel])
    path = os.path.join(directory, "input.lst")
    _write_list(path, ["* material records" if dimensions == 3 else "* 2D material records"] + records)
    return path


def cube_panels(name, side):
    """Six faces with outward normals, needed for the solver's averaged panel normals."""
    def v(x, y, z):
        return (x * side, y * side, z * side)

    faces = [
        [v(0, 0, 0), v(0, 1, 0), v(1, 1, 0), v(1, 0, 0)],
        [v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)],
        [v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)],
        [v(0, 1, 0), v(0, 1, 1), v(1, 1, 1), v(1, 1, 0)],
        [v(0, 0, 0), v(0, 0, 1), v(0, 1, 1), v(0, 1, 0)],
        [v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)],
    ]
    return [("Q", name, face) for face in faces]


def write_plate_over_substrate(directory, triangulated=False):
    """Original #10/#23 plate fixture, including its historical Q face ordering.

    Coordinates are the original micrometer numbers, interpreted by FasterCap as
    meters. Divide reported capacitances by 1e6 for the physical farad values.
    This reproducer has no analytic reference; do not substitute the cube helper's
    corrected outward face orientation when comparing historical platform results.
    With triangulated=True, split each planar face into strips of aspect at most
    ten, then explicit triangles. This preserves the physical surfaces but uses
    a deliberately different mesh, independent of the solver's Q preprocessing.
    """
    def box(name, lower, upper):
        def v(x, y, z):
            return tuple(lo + (hi - lo) * t for lo, hi, t in zip(lower, upper, (x, y, z)))
        faces = [
            [v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0)],
            [v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)],
            [v(0, 0, 0), v(1, 0, 0), v(1, 0, 1), v(0, 0, 1)],
            [v(0, 1, 0), v(1, 1, 0), v(1, 1, 1), v(0, 1, 1)],
            [v(0, 0, 0), v(0, 1, 0), v(0, 1, 1), v(0, 0, 1)],
            [v(1, 0, 0), v(1, 1, 0), v(1, 1, 1), v(1, 0, 1)],
        ]
        return [("Q", name, face) for face in faces]

    for filename, title, name, lower, upper in (
        ("substrate.txt", "substrate", "substrate", (-9, -9, -0.43), (109, 109, -0.1)),
        ("plate.txt", "plate", "plate", (0, 0, 0.936), (100, 100, 1.036)),
        ("oxide.txt", "field oxide", "oxide", (-9, -9, 0), (109, 109, 0.326)),
    ):
        panels = box(name, lower, upper)
        if triangulated:
            triangles = []
            for _, _, vertices in panels:
                def length(a, b):
                    return math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))
                if length(vertices[0], vertices[1]) < length(vertices[1], vertices[2]):
                    vertices = vertices[1:] + vertices[:1]
                count = max(1, math.ceil(length(vertices[0], vertices[1]) /
                                        (10 * length(vertices[1], vertices[2]))))
                def point(a, b, index):
                    return tuple(x + (y-x) * index / count for x, y in zip(a, b))
                for i in range(count):
                    a = point(vertices[0], vertices[1], i)
                    b = point(vertices[0], vertices[1], i+1)
                    c = point(vertices[3], vertices[2], i+1)
                    d = point(vertices[3], vertices[2], i)
                    triangles.extend([("T", name, [a, b, c]), ("T", name, [a, c, d])])
            panels = triangles
        _write_panels(os.path.join(directory, filename), title, panels)
    path = os.path.join(directory, "plate_over_substrate.lst")
    _write_list(path, ["C substrate.txt 3.9 0 0 0", "C plate.txt 3.9 0 0 0",
                       "D oxide.txt 3.9 4.632 0 0 0 50 50 0.163 -"])
    return path


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
COAT_RADIUS = 2.0
COAT_PERMITTIVITY = 4.0
# 2D coaxial cable: inner conductor of RADIUS, coat, outer conductor
COAX_COAT_RADIUS = 1.5
COAX_OUTER_RADIUS = 2.0
COAX_PERMITTIVITY = 10.0
COAX_SEGMENTS = 256


def circle_segments(name, radius, count):
    """Sides of the polygon inscribed in a circle around the origin, as 2D 'S' segments."""
    points = [(radius * math.cos(2.0 * math.pi * i / count), radius * math.sin(2.0 * math.pi * i / count))
              for i in range(count + 1)]
    return [("S", name, [points[i], points[i + 1]]) for i in range(count)]


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


def write_two_cubes(directory):
    """Two cubes beside each other along x, one conductor in one panel file."""
    pitch = SIDE + GAP
    second = [(kind, name, [(x + pitch, y, z) for x, y, z in vertices])
              for kind, name, vertices in cube_panels("c", SIDE)]
    _write_panels(os.path.join(directory, "two_cubes.txt"), "two cubes",
                  cube_panels("c", SIDE) + second)
    path = os.path.join(directory, "two_cubes.lst")
    _write_list(path, ["C two_cubes.txt 1.0 0 0 0"])
    return path


def write_two_cubes_merged(directory):
    """The two cubes as two conductors, merged into one by '+'."""
    _write_panels(os.path.join(directory, "cube.txt"), "cube", cube_panels("c", SIDE))
    path = os.path.join(directory, "two_cubes_merged.lst")
    _write_list(path, ["C cube.txt 1.0 0 0 0 +",
                       "C cube.txt 1.0 %g 0 0" % (SIDE + GAP)])
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


def write_coated_sphere_local_reference(directory):
    """The coated sphere with a reference point outside each coat panel, as klayout-pex writes them."""
    _write_panels(os.path.join(directory, "sphere.txt"), "sphere", icosphere_panels("s", RADIUS, 3))
    with open(os.path.join(directory, "coat.txt"), "w") as f:
        f.write("0 coat, reference points outside\n")
        for _, name, vertices in icosphere_panels("d", COAT_RADIUS, 3):
            outside = [1.1 * sum(v[k] for v in vertices) / 3.0 for k in range(3)]
            f.write("T %s %s\n" % (name, _fmt([c for v in vertices for c in v] + outside)))
    path = os.path.join(directory, "coated_sphere_local_reference.lst")
    # the panels' reference points are on the side of the first permittivity, the vacuum
    _write_list(path, ["C sphere.txt %g 0 0 0" % COAT_PERMITTIVITY,
                       "D coat.txt 1.0 %g 0 0 0 0 0 0" % COAT_PERMITTIVITY])
    return path


def write_coated_sphere_lossy(directory):
    """The coated sphere with a lossy coat, LOSSY_PERMITTIVITY written as FasterCap reads it."""
    write_coated_sphere(directory)
    permittivity = "%g-j%g" % (LOSSY_PERMITTIVITY.real, -LOSSY_PERMITTIVITY.imag)
    path = os.path.join(directory, "coated_sphere_lossy.lst")
    _write_list(path, ["C sphere.txt %s 0 0 0" % permittivity,
                       "D coat.txt 1.0 %s 0 0 0 0 0 0 -" % permittivity])
    return path


LOSSY_PERMITTIVITY = 4.0 - 1.0j
COATED_SPHERE_LOSSY_C = FOUR_PI_EPS0 / ((1.0 / RADIUS - 1.0 / COAT_RADIUS) / LOSSY_PERMITTIVITY
                                        + 1.0 / COAT_RADIUS)


def permittivity_text(value):
    """FasterCap writes the imaginary unit before the magnitude, e.g. '4-j1'."""
    value = complex(value)
    if value.imag == 0:
        return "%g" % value.real
    return "%g%sj%g" % (value.real, "-" if value.imag < 0 else "+", abs(value.imag))


def nested_dielectric_cubes_geometry(permittivity):
    """Planar nested cubes: the closed outer conductor shields the inner one."""
    def write(directory):
        lines = []
        er = permittivity_text(permittivity)
        for name, side, offset in [("inner", 1.0, (.15,-.1,.05)),
                                   ("coat", 2.0, (0,0,0)), ("outer", 3.0, (0,0,0))]:
            panels = [(kind, label, [tuple(v[k]-side/2+offset[k] for k in range(3))
                                    for v in vertices])
                      for kind, label, vertices in cube_panels(name, side)]
            _write_panels(os.path.join(directory, name+".txt"), name, panels)
            lines.append("D coat.txt 1 %s 0 0 0 0 0 0 -" % er if name == "coat"
                         else "C %s.txt %s 0 0 0" % (name, er if name == "inner" else "1"))
        path = os.path.join(directory, "cubes.lst")
        _write_list(path, lines)
        return path
    return Geometry(write, 2, {}, [])


def layered_spheres_geometry(permittivities, interface_radii, depth=3, split_coat=False,
                             offset=(0.0, 0.0, 0.0), split_half=False):
    """Nested conductors with concentric dielectric layers between radii 1 and 3.

    Offset moves the inner conductor and all interfaces. Splitting changes only
    the first interface's triangulation, leaving its physical surface unchanged.
    """
    radii = [1.0, *interface_radii, 3.0]
    assert len(permittivities) == len(radii) - 1
    assert all(a < b for a, b in zip(radii, radii[1:]))
    names = ["inner", *["coat" if i == 0 else "coat%d" % (i+1)
                        for i in range(len(interface_radii))], "outer"]

    def write(directory):
        for name, radius in zip(names, radii):
            panels = icosphere_panels(name, radius, depth)
            if split_coat and name == "coat":
                # Bisect facets without moving their support plane. Splitting
                # only one hemisphere also exercises unequal source-panel areas.
                split = []
                for kind, label, (a, b, c) in panels:
                    if split_half and a[0] + b[0] + c[0] <= 0:
                        split.append((kind, label, (a, b, c)))
                        continue
                    midpoint = tuple((x+y)/2 for x, y in zip(a, b))
                    split.extend([(kind, label, (a, midpoint, c)),
                                  (kind, label, (midpoint, b, c))])
                panels = split
            if name != "outer":
                panels = [(kind, label, [tuple(x + shift for x, shift in zip(v, offset))
                                         for v in vertices])
                          for kind, label, vertices in panels]
            _write_panels(os.path.join(directory, name + ".txt"), name, panels)
        path = os.path.join(directory, "coated_spheres.lst")
        lines = ["C inner.txt %s 0 0 0" % permittivity_text(permittivities[0])]
        for i, name in enumerate(names[1:-1]):
            # The shared center is on the inner ('-') side of each interface.
            lines.append("D %s.txt %s %s 0 0 0 %s -" %
                         (name, permittivity_text(permittivities[i+1]),
                          permittivity_text(permittivities[i]), _fmt(offset)))
        lines.append("C outer.txt %s 0 0 0" % permittivity_text(permittivities[-1]))
        _write_list(path, lines)
        return path

    # Integrate E = Q/(4*pi*eps*r*r) through each annulus in series. The outer
    # conductor adds capacitance to infinity in the outermost medium.
    coupling = FOUR_PI_EPS0 / sum((1.0/a - 1.0/b) / er
                                  for a, b, er in zip(radii, radii[1:], permittivities))
    analytic = {(1, 1): coupling, (1, 2): -coupling, (2, 1): -coupling,
                (2, 2): coupling + FOUR_PI_EPS0 * permittivities[-1] * radii[-1]}
    if any(offset):
        # Reciprocity and shielding still hold; the spherical formula does not.
        analytic = {}
    return Geometry(write, 2, analytic, [])


def coated_spheres_geometry(permittivity, coat_radius=2.0, depth=3, split_coat=False,
                            offset=(0.0, 0.0, 0.0), split_half=False):
    """Nested conductors with one coat; offset moves the inner conductor and coat."""
    return layered_spheres_geometry([permittivity, 1.0], [coat_radius], depth,
                                    split_coat, offset, split_half)


def write_three_cubes_single_file(directory):
    """The three cubes in one input file, the cube as an internal 'File' section."""
    path = os.path.join(directory, "three_cubes_single_file.lst")
    pitch = SIDE + GAP
    with open(path, "w") as f:
        f.write("* three cubes, one file\n")
        f.write("C cube 1.0 0 0 0\nC cube 1.0 %g 0 0\nC cube 1.0 0 %g 0\nEnd\n" % (pitch, pitch))
        f.write("File cube\n0 cube\n")
        for kind, name, vertices in cube_panels("c", SIDE):
            f.write("%s %s %s\n" % (kind, name, _fmt(c for v in vertices for c in v)))
        f.write("End\n")
    return path


def write_three_cubes_renamed(directory):
    """The three cubes, their panels' conductor renamed by 'N' from c to box."""
    _write_panels(os.path.join(directory, "box.txt"), "cube", cube_panels("c", SIDE))
    with open(os.path.join(directory, "box.txt"), "a") as f:
        f.write("N c box\n")
    pitch = SIDE + GAP
    path = os.path.join(directory, "three_cubes_renamed.lst")
    _write_list(path, ["C box.txt 1.0 0 0 0",
                       "C box.txt 1.0 %g 0 0" % pitch,
                       "C box.txt 1.0 0 %g 0" % pitch])
    return path


def write_two_cubes_named(directory):
    """The two cubes as conductors a and b in one file, merged by 'N b a'."""
    pitch = SIDE + GAP
    second = [(kind, "b", [(x + pitch, y, z) for x, y, z in vertices])
              for kind, _, vertices in cube_panels("b", SIDE)]
    _write_panels(os.path.join(directory, "named_cubes.txt"), "cubes a and b",
                  cube_panels("a", SIDE) + second)
    with open(os.path.join(directory, "named_cubes.txt"), "a") as f:
        f.write("N b a\n")
    path = os.path.join(directory, "two_cubes_named.lst")
    _write_list(path, ["C named_cubes.txt 1.0 0 0 0"])
    return path


def write_near_triangles(directory):
    """Two finite, unequal panels closer than AUTOPANEL_EPS in parallel planes."""
    panels = [("T", "c", [(-1, -1, 0), (1, -1, 0), (0, 2, 0)]),
              ("T", "c", [(-2, -1, 1e-13), (2, -1, 1e-13), (0, 2, 1e-13)])]
    _write_panels(os.path.join(directory, "near.txt"), "near parallel triangles", panels)
    path = os.path.join(directory, "near_triangles.lst")
    _write_list(path, ["C near.txt 1.0 0 0 0"])
    return path


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


def write_bowtie_quad(directory):
    """A Q panel whose sides 1-2 and 3-4 cross: two triangles touching at a point."""
    _write_panels(os.path.join(directory, "bowtie.txt"), "bow-tie",
                  [("Q", "c", [(0, 0, 0), (2, 2, 0), (2, 0, 0), (0, 1, 0)])])
    path = os.path.join(directory, "bowtie_quad.lst")
    _write_list(path, ["C bowtie.txt 1.0 0 0 0"])
    return path


def write_coated_coax(directory, segments=COAX_SEGMENTS):
    """2D coaxial cable whose inner conductor is coated with a dielectric."""
    _write_panels(os.path.join(directory, "inner.txt"), "inner conductor",
                  circle_segments("inner", RADIUS, segments))
    _write_panels(os.path.join(directory, "coat.txt"), "coat",
                  circle_segments("coat", COAX_COAT_RADIUS, segments))
    _write_panels(os.path.join(directory, "outer.txt"), "outer conductor",
                  circle_segments("outer", COAX_OUTER_RADIUS, segments))
    path = os.path.join(directory, "coated_coax.lst")
    # "2D" in the first line selects the 2D solver; the reference point, the center, is on the
    # side of the first permittivity, the coat
    _write_list(path, ["* 2D coated coax",
                       "C inner.txt %g 0 0" % COAX_PERMITTIVITY,
                       "D coat.txt %g 1.0 0 0 0 0" % COAX_PERMITTIVITY,
                       "C outer.txt 1.0 0 0"])
    return path


# per meter: the coat, 2*pi*eps0*er / ln(b/a), in series with the vacuum, 2*pi*eps0 / ln(c/b)
COATED_COAX_C = 2.0 * math.pi * EPS0 / (math.log(COAX_COAT_RADIUS / RADIUS) / COAX_PERMITTIVITY
                                       + math.log(COAX_OUTER_RADIUS / COAX_COAT_RADIUS))


GEOMETRIES = {
    "near_triangles": Geometry(write_near_triangles, 1, {}, []),
    "layered_spheres_flux": layered_spheres_geometry([10.0, 3.0, 1.0], [1.6, 2.3], depth=2),
    "coated_spheres_flux_uneven": coated_spheres_geometry(100.0, depth=2,
                                                          split_coat=True, split_half=True),
    "coated_spheres_flux": coated_spheres_geometry(10.0, depth=2),
    "coated_spheres_flux_lossy": coated_spheres_geometry(10.0-2.0j, depth=2),
    "coated_spheres_flux_high_contrast": coated_spheres_geometry(100.0, depth=2),
    "coated_spheres_flux_fine": coated_spheres_geometry(10.0, depth=3),
    "coated_spheres_flux_offset": coated_spheres_geometry(100.0, depth=2,
                                                          offset=(0.4, -0.2, 0.15)),
    "coated_spheres_flux_thin": coated_spheres_geometry(100.0, coat_radius=1.125, depth=2),
    "cube": Geometry(write_cube, 1, {(1, 1): CUBE_FACTOR * FOUR_PI_EPS0 * SIDE}, []),
    "three_cubes": Geometry(write_three_cubes, 3, {}, [((1, 2), (1, 3)), ((2, 2), (3, 3))]),
    "two_cubes": Geometry(write_two_cubes, 1, {}, []),
    "two_cubes_merged": Geometry(write_two_cubes_merged, 1, {}, []),
    "sphere": Geometry(write_sphere, 1, {(1, 1): FOUR_PI_EPS0 * RADIUS}, []),
    "concentric_spheres": Geometry(write_concentric_spheres, 2,
                                   concentric_spheres_analytic(INNER_RADIUS, OUTER_RADIUS), []),
    "coated_sphere": Geometry(write_coated_sphere, 1, {(1, 1): COATED_SPHERE_C}, []),
    "coated_sphere_local_reference": Geometry(write_coated_sphere_local_reference, 1, {(1, 1): COATED_SPHERE_C}, []),
    "coated_sphere_lossy": Geometry(write_coated_sphere_lossy, 1, {(1, 1): COATED_SPHERE_LOSSY_C}, []),
    "three_cubes_single_file": Geometry(write_three_cubes_single_file, 3, {}, [((1, 2), (1, 3)), ((2, 2), (3, 3))]),
    "three_cubes_renamed": Geometry(write_three_cubes_renamed, 3, {}, [((1, 2), (1, 3)), ((2, 2), (3, 3))]),
    "two_cubes_named": Geometry(write_two_cubes_named, 1, {}, []),
    "degenerate_triangle": Geometry(write_degenerate_triangle, 1, {}, []),
    "patran_dielectric": Geometry(write_patran_dielectric, 1, {}, []),
    "skewed_quad": Geometry(write_skewed_quad, 1, {}, []),
    "skewed_quad_triangles": Geometry(write_skewed_quad_triangles, 1, {}, []),
    "bowtie_quad": Geometry(write_bowtie_quad, 1, {}, []),
    # the 2D solver prints no row for the last conductor, the reference
    "coated_coax": Geometry(write_coated_coax, 1, {(1, 1): COATED_COAX_C}, []),
}
