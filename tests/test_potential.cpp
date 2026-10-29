// Unit tests of CPotential: the analytic formulas against closed forms, the numerical
// quadratures against the analytic formulas or a quadrature of their own, and the field
// routines against finite differences of the potential.
//
// Coefficients of potential (Auto, QAuto, AutoNumerical, Mutual*) include 1/(4*pi*eps0);
// the potentials and fields of a single panel (Potential*, InsidePotential, CornerPotential,
// EnField*) don't. Quadrature rules 14, 15, 18 and 19 have no points: the routines return 0
// for them, so the tests use rule 17, which the refinement uses.

#include <cmath>
#include <cstdio>

#include "Potential.h"

// FasterCap's own 4*pi*eps0, so that the closed forms check the formulas, not the constant
static const double COULOMB = 1.0 / FOUR_PI_TIMES_E0;
static const int RULE = 17;

static int failures = 0;

static void check(const char *what, double obtained, double expected, double tolerance)
{
    double difference = std::fabs(obtained - expected) / std::fabs(expected);
    if (!(difference <= tolerance)) {
        std::printf("FAILED %s: %.10g, expected %.10g (%.2g, tolerance %.2g)\n",
                    what, obtained, expected, difference, tolerance);
        failures++;
    }
    else {
        std::printf("ok     %s (%.2g)\n", what, difference);
    }
}

static void triangle(C3DVector_float t[3], double x0, double y0, double z0, double x1, double y1, double z1,
                     double x2, double y2, double z2)
{
    t[0].pos(x0, y0, z0);
    t[1].pos(x1, y1, z1);
    t[2].pos(x2, y2, z2);
}

static double area(const C3DVector_float t[3])
{
    double ux = t[1].x - t[0].x, uy = t[1].y - t[0].y, uz = t[1].z - t[0].z;
    double wx = t[2].x - t[0].x, wy = t[2].y - t[0].y, wz = t[2].z - t[0].z;
    double cx = uy * wz - uz * wy, cy = uz * wx - ux * wz, cz = ux * wy - uy * wx;
    return 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
}

// mean of f over a triangle: centroids of its n*n similar sub-triangles, equal weights
template <class F> static double meanOverTriangle(const C3DVector_float t[3], int n, F f)
{
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; i + j < n; j++) {
            // upright sub-triangle (i,j), (i+1,j), (i,j+1) and, if any, the inverted one
            for (int inverted = 0; inverted < 2; inverted++) {
                if (inverted && i + j + 2 > n) {
                    break;
                }
                double a = (i + (inverted ? 2.0 : 1.0) / 3.0) / n;
                double b = (j + (inverted ? 2.0 : 1.0) / 3.0) / n;
                C3DVector p(t[0].x + a * (t[1].x - t[0].x) + b * (t[2].x - t[0].x),
                            t[0].y + a * (t[1].y - t[0].y) + b * (t[2].y - t[0].y),
                            t[0].z + a * (t[1].z - t[0].z) + b * (t[2].z - t[0].z));
                sum += f(p);
            }
        }
    }
    return sum / (n * n);
}

static C3DVector centroid(const C3DVector_float t[3])
{
    return C3DVector((t[0].x + t[1].x + t[2].x) / 3.0, (t[0].y + t[1].y + t[2].y) / 3.0,
                     (t[0].z + t[1].z + t[2].z) / 3.0);
}

static double distance(C3DVector a, C3DVector b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

static void shifted(C3DVector_float out[3], const C3DVector_float t[3], double dz)
{
    for (int k = 0; k < 3; k++) {
        out[k].pos(t[k].x, t[k].y, t[k].z + dz);
    }
}

static void testTriangleClosedForms(CPotential &pot)
{
    // equilateral triangle of side 1: the potential at the centroid is 4 ln(2 + sqrt(3)), at a
    // vertex 2 ln(3), and the mean over all pairs of points of 1/r is 4 ln(3), as sums over the
    // sides of d ln((l+ + R+) / (l- + R-)) give. The vertices are floats: 1e-7
    C3DVector_float eq[3];
    triangle(eq, 0, 0, 0, 1, 0, 0, 0.5, std::sqrt(3.0) / 2, 0);
    C3DVector vertex(0, 0, 0);
    check("Auto(), equilateral triangle", pot.Auto(eq), COULOMB * 4.0 * std::log(3.0), 1e-6);
    check("PotentialOpt() at the centroid of an equilateral triangle", pot.PotentialOpt(centroid(eq), eq),
          4.0 * std::log(2.0 + std::sqrt(3.0)), 1e-6);
    check("PotentialOpt() at a vertex of an equilateral triangle", pot.PotentialOpt(vertex, eq),
          2.0 * std::log(3.0), 1e-6);
    check("CornerPotential(), equilateral triangle", pot.CornerPotential(eq), 2.0 * std::log(3.0), 1e-6);
    check("PotentialOpt() without dividing by the area", pot.PotentialOpt(centroid(eq), eq, false),
          std::sqrt(3.0) * std::log(2.0 + std::sqrt(3.0)), 1e-6);

    // unit square: the mean of 1/r over all pairs of points is 4 (ln(1 + sqrt(2)) - (sqrt(2) - 1) / 3);
    // QAuto() integrates the coupling of its two triangles numerically: 1.3e-4
    C3DVector square[4];
    square[0].pos(0, 0, 0);
    square[1].pos(1, 0, 0);
    square[2].pos(1, 1, 0);
    square[3].pos(0, 1, 0);
    check("QAuto(), unit square", pot.QAuto(square),
          COULOMB * 4.0 * (std::log(1.0 + std::sqrt(2.0)) - (std::sqrt(2.0) - 1.0) / 3.0), 5e-4);
}

static void testTrianglePotential(CPotential &pot)
{
    C3DVector_float tri[3];
    triangle(tri, 0, 0, 0, 1.3, 0.1, 0, 0.4, 0.9, 0);
    C3DVector inside(0.5, 0.3, 0), above(0.5, 0.3, 0.7), far(0.5, 0.3, 1000.0);

    // 1/r is smooth at 0.7 above the triangle: 200*200 sub-triangles give 1e-5
    double mean = meanOverTriangle(tri, 200, [&](C3DVector p) { return 1.0 / distance(above, p); });
    check("PotentialOpt() above a triangle", pot.PotentialOpt(above, tri), mean, 1e-4);
    // the same integral, rounded differently in single precision: 3e-8
    check("Potential() equals PotentialOpt()", pot.Potential(above, tri), pot.PotentialOpt(above, tri), 1e-6);
    check("InsidePotential() equals PotentialOpt()", pot.InsidePotential(inside, tri),
          pot.PotentialOpt(inside, tri), 1e-6);
    // the quadrupole term is (size / distance)^2 = 1e-6 of the monopole
    check("PotentialOpt() far from a triangle", pot.PotentialOpt(far, tri), 1.0 / distance(far, centroid(tri)), 1e-5);
    // rule 17 has precision 9: 3e-6 at this distance
    check("PotentialNumerical() above a triangle", pot.PotentialNumerical(above, tri, RULE),
          pot.PotentialOpt(above, tri), 1e-5);

    // the field along the normal is minus the derivative of the potential: 2.4e-5
    const double delta = 1e-5;
    C3DVector up(above.x, above.y, above.z + delta), down(above.x, above.y, above.z - delta);
    C3DVector normal(0, 0, 1);
    check("EnFieldNumerical() above a triangle", pot.EnFieldNumerical(above, tri, normal, RULE),
          -(pot.PotentialOpt(up, tri) - pot.PotentialOpt(down, tri)) / (2.0 * delta), 1e-4);

    // the mean of the potential over the triangle itself: the potential has a logarithmic
    // derivative at the sides, 64*64 sub-triangles give 1e-4; the quadrature of AutoNumerical()
    // misses the singularity by 8e-4
    double self = COULOMB * meanOverTriangle(tri, 64, [&](C3DVector p) { return pot.PotentialOpt(p, tri); });
    check("Auto(), triangle", pot.Auto(tri), self, 1e-3);
    check("AutoNumerical() with InsidePotential()", pot.AutoNumerical(tri, RULE, 0, true), pot.Auto(tri), 2e-3);
    check("AutoNumerical() with Potential()", pot.AutoNumerical(tri, RULE, 1, true), pot.Auto(tri), 2e-3);
}

static void testTriangleMutual(CPotential &pot)
{
    // a triangle and a shifted copy 1 above it: the potential of the second is smooth on the
    // first, so 200*200 sub-triangles give 1e-6
    C3DVector_float tri1[3], tri2[3];
    triangle(tri1, 0, 0, 0, 1.3, 0.1, 0, 0.4, 0.9, 0);
    triangle(tri2, 0.2, -0.1, 1, 1.5, 0, 1, 0.6, 0.8, 1);
    double mutual = COULOMB * meanOverTriangle(tri1, 200, [&](C3DVector p) { return pot.PotentialOpt(p, tri2); });
    check("Mutual_9thOrd_HalfNum()", pot.Mutual_9thOrd_HalfNum(tri1, tri2), mutual, 1e-5);
    // the areas from float vertices: 4e-8
    check("Mutual_9thOrd_HalfNum() without dividing by the areas", pot.Mutual_9thOrd_HalfNum(tri1, tri2, false),
          pot.Mutual_9thOrd_HalfNum(tri1, tri2) * area(tri1) * area(tri2), 1e-6);
    check("MutualHalfNumerical()", pot.MutualHalfNumerical(tri1, tri2, RULE), mutual, 1e-5);
    check("MutualFullNumerical()", pot.MutualFullNumerical(tri1, tri2, RULE), mutual, 1e-5);
    check("Mutual_4thOrd_FullNum()", pot.Mutual_4thOrd_FullNum(tri1, tri2), mutual, 1e-4);
    // low-order rules: 0.19 % and 0.33 %
    check("Mutual_2thOrd_FullNum()", pot.Mutual_2thOrd_FullNum(tri1, tri2), mutual, 1e-2);
    check("Mutual_3thOrd_FullNum()", pot.Mutual_3thOrd_FullNum(tri1, tri2), mutual, 1e-2);

    // the field on the first triangle along its normal is minus the derivative of the mutual
    // coefficient when the first triangle moves along it
    const double delta = 1e-4;
    C3DVector_float up[3], down[3];
    shifted(up, tri1, delta);
    shifted(down, tri1, -delta);
    double field = -(pot.Mutual_9thOrd_HalfNum(up, tri2) - pot.Mutual_9thOrd_HalfNum(down, tri2)) / (2.0 * delta);
    C3DVector normal(0, 0, 1);
    check("MutualDFullNumerical()", pot.MutualDFullNumerical(tri1, tri2, normal, RULE), field, 1e-4);
    // divided differences of step rdist / 20, as the refinement uses them
    double h = distance(centroid(tri1), centroid(tri2)) / 20.0;
    check("MutualDHalfNumerical()", pot.MutualDHalfNumerical(tri1, tri2, normal, h, RULE), field, 5e-3);
    // 2nd order: 0.74 %
    check("MutualD_2thOrd_FullNum()", pot.MutualD_2thOrd_FullNum(tri1, tri2, normal), field, 2e-2);

    // two unit squares 1 apart: the sum over the pairs of their triangles
    C3DVector sq1[4], sq2[4];
    sq1[0].pos(0, 0, 0);
    sq1[1].pos(1, 0, 0);
    sq1[2].pos(1, 1, 0);
    sq1[3].pos(0, 1, 0);
    for (int k = 0; k < 4; k++) {
        sq2[k].pos(sq1[k].x + 0.3, sq1[k].y, 1.0);
    }
    C3DVector_float t[2][2][3];
    for (int s = 0; s < 2; s++) {
        C3DVector *q = s ? sq2 : sq1;
        triangle(t[s][0], q[0].x, q[0].y, q[0].z, q[1].x, q[1].y, q[1].z, q[2].x, q[2].y, q[2].z);
        triangle(t[s][1], q[2].x, q[2].y, q[2].z, q[0].x, q[0].y, q[0].z, q[3].x, q[3].y, q[3].z);
    }
    double pairs = 0.0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            pairs += pot.Mutual_9thOrd_HalfNum(t[0][i], t[1][j], false);
        }
    }
    check("QMutual_9thOrd_HalfNum()", pot.QMutual_9thOrd_HalfNum(sq1, sq2), pairs, 1e-6);
}

static void testSegment(CPotential &pot)
{
    // segment from -1/2 to 1/2 on the x axis: minus the mean of ln(r) over it, from
    // the integral u ln(sqrt(u^2 + y^2)) - u + y atan(u / y) along u = x' - x
    C2DVector_float segment[2];
    segment[0].pos(-0.5, 0);
    segment[1].pos(0.5, 0);
    const double x = 0.2;
    for (double y : {0.1, 1.0, 10.0}) {
        auto primitive = [&](double u) { return u * std::log(std::sqrt(u * u + y * y)) - u + y * std::atan(u / y); };
        C2DVector r(x, y);
        char what[80];
        std::snprintf(what, sizeof(what), "PotentialOpt() of a segment at y = %g", y);
        check(what, pot.PotentialOpt(r, segment), -(primitive(0.5 - x) - primitive(-0.5 - x)), 1e-9);
    }

    // the field along a normal is minus the derivative of the potential
    const double delta = 1e-5;
    C2DVector r(0.2, 0.3), normal(0.6, 0.8);
    C2DVector plus(r.x + delta * normal.x, r.y + delta * normal.y), minus(r.x - delta * normal.x, r.y - delta * normal.y);
    check("EnField() of a segment", pot.EnField(r, segment, normal),
          -(pot.PotentialOpt(plus, segment) - pot.PotentialOpt(minus, segment)) / (2.0 * delta), 1e-6);
}

int main()
{
    CPotential pot;

    testTriangleClosedForms(pot);
    testTrianglePotential(pot);
    testTriangleMutual(pot);
    testSegment(pot);

    if (failures) {
        std::printf("%d checks failed\n", failures);
        return 1;
    }
    return 0;
}
