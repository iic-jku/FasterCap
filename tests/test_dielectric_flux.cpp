// Closed-surface Gauss identities on exactly representable polyhedra. No mesh or solve.
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include "Potential.h"

typedef decltype(CAutoPanel::m_clsVertex) Triangle;
static int failures = 0;

static double dimension(Triangle t, C3DVector &normal)
{
    C3DVector a(t[0].x, t[0].y, t[0].z), b(t[1].x, t[1].y, t[1].z);
    C3DVector c(t[2].x, t[2].y, t[2].z);
    normal = CrossProd(b-a, c-a);
    double twiceArea = Mod(normal);
    normal = normal / twiceArea;
    return twiceArea / 2.0;
}

static void check(const char *name, double got, double want, double tolerance)
{
    if (!std::isfinite(got) || std::fabs(got-want) > tolerance) {
        std::printf("FAILED %s: %.17g, expected %.17g, absolute tolerance %.3g\n",
                    name, got, want, tolerance);
        ++failures;
    }
}

static void closedCube(CPotential &pot, double scale, double offset, int rule, bool inward, bool skew = false)
{
    // Every face is counterclockwise from outside; split into 12 triangles.
    const int xyz[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                          {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    const int faces[6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},
                            {3,7,6,2},{0,4,7,3},{1,2,6,5}};
    Triangle targets[12];
    int nt = 0;
    for (int face=0; face<6; ++face) {
        for (int part=0; part<2; ++part) {
            int indices[3] = {faces[face][0], faces[face][part+1], faces[face][part+2]};
            for (int k=0; k<3; ++k) {
                int index=indices[inward ? 2-k : k];
                targets[nt][k].pos(offset+scale*xyz[index][0], offset+scale*xyz[index][1],
                                   offset+scale*xyz[index][2]);
            }
            ++nt;
        }
    }
    if (skew) {
        for (int i=0;i<nt;++i) for (int k=0;k<3;++k) {
            double x=targets[i][k].x, y=targets[i][k].y, z=targets[i][k].z;
            targets[i][k].pos(x+z/4, y+x/2, z+y/8);
        }
    }
    // Triangles inside, very near either side of the top face, and on its interior.
    const double heights[] = {0.0, 1.0-1.0/1024, 1.0+1.0/1024, 1.0};
    for (double height : heights) {
        Triangle source;
        source[0].pos(offset-scale/2, offset-scale/4, offset+scale*height);
        source[1].pos(offset+scale/4, offset-scale/4, offset+scale*height);
        source[2].pos(offset-scale/2, offset+scale/4, offset+scale*height);
        if (skew) {
            for (int k=0;k<3;++k) {
                double x=source[k].x, y=source[k].y, z=source[k].z;
                source[k].pos(x+z/4, y+x/2, z+y/8);
            }
        }
        double flux=0;
        for (int i=0; i<nt; ++i) {
            C3DVector normal;
            double area=dimension(targets[i],normal);
            double coefficient=pot.MutualDHalfNumerical(targets[i],source,normal,scale/20,rule);
            flux += area*coefficient*FOUR_PI_TIMES_E0/(4.0*PI);
        }
        double expected=height<1 ? 1.0 : height>1 ? 0.0 : 0.5;
        if (inward) expected=-expected;
        check("closed cube flux",flux,expected,2e-12);
    }
}

static void arbitraryDirection(CPotential &pot)
{
    Triangle target, source;
    target[0].pos(0,0,0); target[1].pos(1,0,0); target[2].pos(0,1,0);
    source[0].pos(0.25,0.125,2); source[1].pos(0.5,0.125,2); source[2].pos(0.25,0.5,2);
    const C3DVector directions[] = {C3DVector(1,0,0),C3DVector(0,1,0),C3DVector(0,0,1),C3DVector(1,-2,3)};
    for (C3DVector direction : directions) {
        double expected=pot.MutualDFullNumerical(target,source,direction,17)*FOUR_PI_TIMES_E0;
        double result=pot.MutualDHalfNumerical(target,source,direction,0.05,17)*FOUR_PI_TIMES_E0;
        check("arbitrary field direction",result,expected,2e-7);
        C3DVector n;
        double areas=dimension(target,n)*dimension(source,n);
        check("area normalization",pot.MutualDHalfNumerical(target,source,direction,0.05,17,false)*FOUR_PI_TIMES_E0,
              result*areas,2e-14);
        CAutoPanel a,b;
        for (int i=0;i<3;++i) { a.m_clsVertex[i]=target[i]; b.m_clsVertex[i]=source[i]; }
        a.CalcPanelGeomPar(); b.CalcPanelGeomPar();
        check("panel overload",pot.MutualDHalfNumerical(a,b,direction,0.05,17)*FOUR_PI_TIMES_E0,result,2e-14);
    }
}

static void coplanarPanel(CPotential &pot)
{
    // Identical source/target triangles are exactly coplanar even when a rounded
    // scalar triple product fails to cancel. Their normal principal value is zero.
    Triangle target, source;
    target[0].pos(-0.8269853591918945, -0.5172462463378906, -0.2020111083984375);
    target[1].pos(0.5741386413574219, 0.3886747360229492, 0.15745925903320312);
    target[2].pos(-0.40696144104003906, 0.04431343078613281, -0.5693292617797852);
    C3DVector normal;
    double area = dimension(target, normal);
    for (int shift = 0; shift < 3; ++shift) {
        for (int i = 0; i < 3; ++i) source[i] = target[(i+shift)%3];
        check("coplanar principal value",
              pot.MutualDHalfNumerical(target, source, normal, 0.05, 17)
                  * area * FOUR_PI_TIMES_E0, 0.0, 2e-12);
    }
    // A tiny but representable gap is not coplanar. It retains the one-sided
    // field jump even though its distance is much smaller than machine epsilon.
    target[0].pos(0,0,0); target[1].pos(1,0,0); target[2].pos(0,1,0);
    area = dimension(target, normal);
    for (double side : {-1.0, 1.0}) {
        for (int i = 0; i < 3; ++i) {
            source[i] = target[i];
            source[i].z = side * std::ldexp(1.0, -60);
        }
        check("separated panel field jump",
              pot.MutualDHalfNumerical(target, source, normal, 0.05, 17)
                  * area * FOUR_PI_TIMES_E0, -side * 2.0 * PI, 2e-12);
    }
}

int main()
{
    CPotential pot;
    for (double scale : {1.0/8,1.0,8.0})
        for (int rule : {1,5,17})
            for (bool inward : {false,true})
                closedCube(pot,scale,0,rule,inward);
    closedCube(pot,1.0,1024.0,17,false);
    closedCube(pot,1.0,0.0,17,false,true);
    arbitraryDirection(pot);
    coplanarPanel(pot);
    if (failures) {
        std::printf("%d dielectric flux checks failed\n",failures);
        return 1;
    }
    std::puts("Closed-surface flux, orientation, scale, translation and field-direction checks passed");
    return 0;
}
