// Investigation of the actual 3D default interaction dispatch, without a solve.
#include "Autorefine.h"
#include <cmath>
#include <cstdio>
#include <initializer_list>

static double worst = 0.0;
static int failures = 0;
static void check(double value, double expected)
{
    double error = std::fabs(value - expected);
    if(error > worst) worst = error;
    if(!std::isfinite(value) || error > 2e-12) ++failures;
}

int main()
{
    CAutoRefine ref;
    ref.m_clsGlobalVars.m_cScheme = AUTOREFINE_COLLOCATION;
    const int xyz[8][3] = {{-1,-1,-1},{1,-1,-1},{1,1,-1},{-1,1,-1},
                          {-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}};
    const int faces[6][4] = {{0,3,2,1},{4,5,6,7},{0,1,5,4},
                            {3,7,6,2},{0,4,7,3},{1,2,6,5}};
    for(double scale : {0.125, 1.0, 8.0}) for(bool outward : {true, false}) {
        CAutoPanel targets[12];
        int nt = 0;
        for(int face=0; face<6; ++face) for(int part=0; part<2; ++part) {
            int indices[3] = {faces[face][0], faces[face][part+1], faces[face][part+2]};
            for(int k=0; k<3; ++k) {
                int v=indices[k];
                targets[nt].m_clsVertex[k].pos(scale*xyz[v][0], scale*xyz[v][1], scale*xyz[v][2]);
            }
            targets[nt].CalcPanelGeomPar();
            targets[nt].m_ucType = AUTOPANEL_IS_DIEL | (outward ? AUTOPANEL_OUTPERM_NORMAL_DIR : 0);
            ++nt;
        }
        // Every unit source on a planar face contributes half its flux. The
        // self-panel principal value is zero; the remaining panels carry 2*pi.
        for(int j=0; j<nt; ++j) {
            double flux = 0.0, pairFlux = 0.0;
            for(int i=0; i<nt; ++i) if(i!=j) {
                double single, first, second;
                if(ref.PotEstimateOpt(&targets[i], &targets[j], single) ||
                   ref.PotEstimateOpt(&targets[i], &targets[j], first, second)) return 2;
                flux += single*targets[i].GetArea()*FOUR_PI_TIMES_E0/(4.0*PI);
                pairFlux += first*targets[i].GetArea()*FOUR_PI_TIMES_E0/(4.0*PI);
                check(first*FOUR_PI_TIMES_E0*scale*scale, single*FOUR_PI_TIMES_E0*scale*scale);
            }
            check(flux, outward ? 0.5 : -0.5);
            check(pairFlux, outward ? 0.5 : -0.5);
        }
        for(bool dielectric : {false,true}) for(double height : {0.0, 1.0-1.0/1024, 1.0+1.0/1024}) {
            CAutoPanel source;
            source.m_clsVertex[0].pos(-scale/2,-scale/4,scale*height);
            source.m_clsVertex[1].pos(scale/4,-scale/4,scale*height);
            source.m_clsVertex[2].pos(-scale/2,scale/4,scale*height);
            source.CalcPanelGeomPar();
            source.m_ucType = AUTOPANEL_OUTPERM_NORMAL_DIR | (dielectric ? AUTOPANEL_IS_DIEL : 0);
            double flux = 0.0, pairFlux = 0.0;
            for(int i=0; i<nt; ++i) {
                double single, first, second;
                if(ref.PotEstimateOpt(&targets[i], &source, single) ||
                   ref.PotEstimateOpt(&targets[i], &source, first, second)) return 2;
                flux += single*targets[i].GetArea()*FOUR_PI_TIMES_E0/(4.0*PI);
                pairFlux += first*targets[i].GetArea()*FOUR_PI_TIMES_E0/(4.0*PI);
            }
            double expected = height<1.0 ? 1.0 : 0.0;
            if(!outward) expected = -expected;
            check(flux,expected); check(pairFlux,expected);
        }
    }
    std::printf("Collocation flux failures: %d; worst error as fraction of total source flux: %.17g\n", failures, worst);
    return failures ? 1 : 0;
}
