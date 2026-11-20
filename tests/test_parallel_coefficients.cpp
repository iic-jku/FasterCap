#include "Autorefine.h"
#include "omp.h"
#include <cmath>
#include <cstdio>
#include <initializer_list>

class CParallelProbe : public CAutoRefine
{
public:
    int Evaluate(CAutoElement *target, CAutoElement *source, double &value)
    {
        CAutoElement **chunk = &source;
        double *coefficients = &value;
        m_pNodes = &target;
        m_pdPanelPtrLinks[0] = &chunk;
        m_dPotCoeffLinks[0] = &coefficients;
        m_ulLinkChunkNum[0] = 1;
        m_ulCurrBlock = 0;
        target->m_ulLinkIndexStart[0] = 0;
        target->m_ulLinkIndexEnd[0] = 1;
        m_ulNumofpotest = 37;
        int status = ComputeLinkCoefficients(0, 1, 0);
        if(m_ulNumofpotest != 37) status = -100;
        m_pNodes = nullptr;
        m_pdPanelPtrLinks[0] = nullptr;
        m_dPotCoeffLinks[0] = nullptr;
        m_ulLinkChunkNum[0] = 0;
        return status;
    }
};

int main()
{
    g_bFCContinue = true;
    CParallelProbe probe;
    probe.m_clsGlobalVars.m_cScheme = AUTOREFINE_COLLOCATION;
    for(int threads : {1, 2, 4}) {
        omp_set_num_threads(threads);
        g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
        CAutoPanel panels[2], group, source;
        for(int i=0; i<2; ++i) {
            panels[i].m_clsVertex[0].pos(i*2,0,0);
            panels[i].m_clsVertex[1].pos(i*2,1,0);
            panels[i].m_clsVertex[2].pos(i*2,0,1);
            panels[i].CalcPanelGeomPar();
        }
        double value = -12345, expected;
        if(probe.PotEstimateOpt(&panels[0], &panels[1], expected)) return 1;
        if(probe.Evaluate(&panels[0], &panels[1], value) || value != expected) return 2;
        value = -12345;
        if(probe.Evaluate(&panels[0], &panels[0], value) != FC_INVALID_RESULT || value != -12345) return 3;

        for(auto &panel : panels) panel.m_ucType |= AUTOPANEL_IS_DIEL;
        group.MakeSuperPanel(&panels[0], &panels[1]);
        source.GetCentroid() = group.GetCentroid();
        value = -12345;
        if(probe.Evaluate(&group, &source, value) != FC_INVALID_RESULT || value != -12345) return 4;

        g_ucSolverType = SOLVERGLOBAL_2DSOLVER;
        CAutoSegment segments[2];
        for(int i=0; i<2; ++i) {
            segments[i].m_clsVertex[0].pos(i*2,0);
            segments[i].m_clsVertex[1].pos(i*2,1);
            segments[i].CalcSegmentGeomPar();
        }
        if(probe.PotEstimateOpt(&segments[0], &segments[1], expected)) return 5;
        if(probe.Evaluate(&segments[0], &segments[1], value) || value != expected) return 6;
        value = -12345;
        if(probe.Evaluate(&segments[0], &segments[0], value) != FC_INVALID_RESULT || value != -12345) return 7;
    }
    std::puts("Parallel coefficients preserve values, reject failed links and defer statistics");
    return 0;
}
