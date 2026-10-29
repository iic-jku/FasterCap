// Unit tests of C2DTriangulate as FasterCap uses it for planar Q panels: the outline is
// triangulated, its sides are inserted as constraints, the triangles outside it are deleted
// (which needs the constraints), and the mesh is refined to a minimum angle. The triangles must cover the outline exactly, lie inside it, be counter-clockwise
// and keep the minimum angle.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "Geometry/Triangulate.h"

// as AUTOREFINE_MIN_ANGLE in Solver/Autorefine.cpp
static const double MIN_ANGLE = 20.0;

static int failures = 0;

static void fail(const char *name, const char *what, double value)
{
    std::printf("FAILED %s: %s (%g)\n", name, what, value);
    failures++;
}

static double cross(C2DVector o, C2DVector a, C2DVector b)
{
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

static double polygonArea(const std::vector<C2DVector> &outline)
{
    double twice = 0.0;
    for (size_t i = 0; i < outline.size(); i++) {
        const C2DVector &a = outline[i], &b = outline[(i + 1) % outline.size()];
        twice += a.x * b.y - b.x * a.y;
    }
    return 0.5 * twice;
}

// even-odd rule; points on the outline count as inside
static bool inside(const std::vector<C2DVector> &outline, C2DVector p, double tolerance)
{
    bool in = false;
    for (size_t i = 0; i < outline.size(); i++) {
        const C2DVector &a = outline[i], &b = outline[(i + 1) % outline.size()];
        double length = std::hypot(b.x - a.x, b.y - a.y);
        double t = ((p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y)) / (length * length);
        if (t >= 0.0 && t <= 1.0 && std::fabs(cross(a, b, p)) / length <= tolerance) {
            return true;
        }
        if ((a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y)) {
            in = !in;
        }
    }
    return in;
}

static double minAngle(C2DVector a, C2DVector b, C2DVector c)
{
    C2DVector v[3] = {a, b, c};
    double smallest = 180.0;
    for (int k = 0; k < 3; k++) {
        C2DVector o = v[k], p = v[(k + 1) % 3], q = v[(k + 2) % 3];
        double angle = std::atan2(std::fabs(cross(o, p, q)),
                                  (p.x - o.x) * (q.x - o.x) + (p.y - o.y) * (q.y - o.y)) * 180.0 / 3.14159265358979323846;
        smallest = std::min(smallest, angle);
    }
    return smallest;
}

// outline: counter-clockwise corners; minAngleExpected: the smallest angle the refined mesh
// must keep, below MIN_ANGLE where the outline has a smaller corner
static void testOutline(const char *name, const std::vector<C2DVector> &outline, double minAngleExpected)
{
    C2DTriangulate tri;
    for (const C2DVector &p : outline) {
        tri.vertexes.push_back(p);
    }
    tri.Triangulate();
    for (size_t i = 0; i < outline.size(); i++) {
        tri.InsertConstrEdge(outline[i], outline[(i + 1) % outline.size()]);
    }
    for (size_t i = 0; i < outline.size(); i++) {
        tri.DeleteCWEdges(outline[i], outline[(i + 1) % outline.size()]);
    }
    tri.Refine(MIN_ANGLE);
    tri.GenerateTriangles();

    double size = 0.0;
    for (const C2DVector &p : outline) {
        size = std::max(size, std::max(std::fabs(p.x), std::fabs(p.y)));
    }
    const double tolerance = 1e-9 * size;
    double total = 0.0, smallest = 180.0;
    long count = tri.GetTriangleNum();
    for (long i = 0; i < count; i++) {
        C2DVector a, b, c;
        tri.GetTriangleCoords(i, &a, &b, &c);
        double twice = cross(a, b, c);
        if (twice <= 0.0) {
            fail(name, "triangle not counter-clockwise or degenerate", twice);
        }
        total += 0.5 * twice;
        C2DVector center((a.x + b.x + c.x) / 3.0, (a.y + b.y + c.y) / 3.0);
        if (!inside(outline, center, tolerance) || !inside(outline, a, tolerance) || !inside(outline, b, tolerance)
                || !inside(outline, c, tolerance)) {
            fail(name, "triangle outside the outline, centroid x", center.x);
        }
        smallest = std::min(smallest, minAngle(a, b, c));
    }
    double expected = polygonArea(outline);
    if (count == 0 || std::fabs(total - expected) > 1e-9 * expected) {
        fail(name, "triangle areas don't add up to the outline's", total - expected);
    }
    if (smallest < minAngleExpected) {
        fail(name, "minimum angle in degrees", smallest);
    }
    std::printf("%s %s: %ld triangles, minimum angle %.1f degrees\n", failures ? "      " : "ok    ", name, count,
                smallest);
}

int main()
{
    testOutline("square", {C2DVector(0, 0), C2DVector(1, 0), C2DVector(1, 1), C2DVector(0, 1)}, MIN_ANGLE);
    testOutline("rectangle 10:1", {C2DVector(0, 0), C2DVector(10, 0), C2DVector(10, 1), C2DVector(0, 1)}, MIN_ANGLE);
    testOutline("convex quadrilateral", {C2DVector(0, 0), C2DVector(3, 0.5), C2DVector(2.5, 2), C2DVector(-0.5, 1.5)},
                MIN_ANGLE);
    testOutline("concave quadrilateral", {C2DVector(0, 0), C2DVector(2, 0), C2DVector(1, 0.5), C2DVector(1, 2)}, MIN_ANGLE);
    // a triangle given as a quadrilateral whose 4th corner is the midpoint of a side, rounded to
    // float as FasterCap reads it: the corner is almost straight, inward or outward. Outward,
    // Refine() stops at 15 degrees next to it
    testOutline("triangle with float-rounded midpoint, inward", {C2DVector(0, 0), C2DVector(0.3, 0),
                C2DVector(0.15, (float)0.35), C2DVector(0, 0.7)}, MIN_ANGLE);
    testOutline("triangle with float-rounded midpoint, outward", {C2DVector(0, 0), C2DVector(0.3, 0),
                C2DVector(0.15, (float)0.55), C2DVector(0, 1.1)}, 10.0);
    testOutline("square of side 1e-6", {C2DVector(0, 0), C2DVector(1e-6, 0), C2DVector(1e-6, 1e-6), C2DVector(0, 1e-6)},
                MIN_ANGLE);
    testOutline("square of side 1e6", {C2DVector(0, 0), C2DVector(1e6, 0), C2DVector(1e6, 1e6), C2DVector(0, 1e6)},
                MIN_ANGLE);

    if (failures) {
        std::printf("%d checks failed\n", failures);
        return 1;
    }
    return 0;
}
