#include "SolveCapacitance.h"
#include <cmath>
#include <cstdio>
#include <limits>

class CBlockProbe : public CSolveCap
{
public:
    CBlockProbe()
    {
        m_pdBlockPrecond = new double[3][SOLVE_MAX_BLOCK_PRECOND_NUM]();
        m_pucBlockPrecondDim = new unsigned char[3]();
        m_pclsSupPotMtx = new double[2][SOLVE_MAX_SUPER_PRECOND_NUM]();
        m_pclsTmpCapMtx = new double[2][SOLVE_MAX_SUPER_PRECOND_NUM]();
        m_clsBlockPrecondElements = new CAutoElement*[2];
    }

    int CheckApplication()
    {
        int failures = 0;
        m_clsMulthier.m_ulPanelNum[0] = 3;
        m_pucBlockPrecondDim[0] = 2;
        m_pucBlockPrecondDim[2] = 1;
        m_pdBlockPrecond[0][0] = 2; m_pdBlockPrecond[0][1] = 1;
        m_pdBlockPrecond[1][0] = -1; m_pdBlockPrecond[1][1] = 3;
        m_pdBlockPrecond[2][0] = 4;
        const double inputs[] = {1,2,3,-2,1,5};
        const double expected[] = {4,5,12,-3,5,20};
        for(bool complex : {false, true}) {
            m_clsGlobalVars.m_ucHasCmplxPerm = complex ? AUTOREFINE_CPLX_PERM : AUTOREFINE_REAL_PERM;
            CLin_Vector q(complex ? 6 : 3), result(q.size());
            for(unsigned long i=0; i<q.size(); ++i) q[i] = inputs[i];
            for(double oldValue : {std::numeric_limits<double>::quiet_NaN(), 123.0}) {
                result = oldValue;
                ComputePrecondVectFast(&result, &q, AUTOREFINE_PRECOND_BLOCK);
                for(unsigned long i=0; i<q.size(); ++i) {
                    if(!std::isfinite(result[i]) || result[i] != expected[i]) {
                        std::printf("FAILED block application complex=%d entry=%lu: %g expected %g\n",
                                    complex, i, result[i], expected[i]);
                        ++failures;
                    }
                }
            }
        }
        return failures;
    }

    int CheckMaterial()
    {
        g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
        m_clsGlobalVars.m_ucHasCmplxPerm = AUTOREFINE_REAL_PERM;
        CAutoPanel panel;
        panel.m_clsVertex[0].pos(0,0,0);
        panel.m_clsVertex[1].pos(1,0,0);
        panel.m_clsVertex[2].pos(0,1,0);
        panel.CalcPanelGeomPar();
        panel.m_ucType |= AUTOPANEL_IS_DIEL;
        panel.m_lIndex[0] = 0;
        CAutoConductor correct, stale;
        correct.m_dInperm[0] = 10; correct.m_dOutperm[0] = 1;
        stale.m_dInperm[0] = 100; stale.m_dOutperm[0] = 1;
        m_clsMulthier.SetCurrentConductor(&correct);
        double diagonal, imaginary;
        m_clsMulthier.SelfPotential(&panel, &diagonal, &imaginary);
        m_clsMulthier.m_clsSelfPotCoeff.newsize(1);
        m_clsMulthier.m_clsSelfPotCoeff[0] = diagonal;
        m_clsMulthier.SetCurrentConductor(&stale);
        m_pCurrCond = &correct;
        m_clsBlockPrecondElements[0] = &panel;
        m_uiBlockPreNum = 1;
        m_ulBlockPreBaseNum = 0;
        ComputeBlockPrecond();
        double identity = m_pdBlockPrecond[0][0] * diagonal;
        if(!std::isfinite(identity) || std::fabs(identity-1.0) > 1e-12) {
            std::printf("FAILED material diagonal: inverse times assembled diagonal = %.17g\n", identity);
            return 1;
        }
        return 0;
    }

    int CheckSegments()
    {
        g_ucSolverType = SOLVERGLOBAL_2DSOLVER;
        CAutoSegment segments[2];
        double matrix[2][2], imaginary;
        m_clsMulthier.m_clsSelfPotCoeff.newsize(2);
        for(int i=0; i<2; ++i) {
            segments[i].m_clsVertex[0].pos(0,3*i);
            segments[i].m_clsVertex[1].pos(1,3*i);
            segments[i].CalcSegmentGeomPar();
            segments[i].m_lIndex[0] = i;
            m_clsMulthier.SelfPotential(&segments[i], &matrix[i][i], &imaginary);
            m_clsMulthier.m_clsSelfPotCoeff[i] = matrix[i][i];
            m_clsBlockPrecondElements[i] = &segments[i];
        }
        m_clsMulthier.PotEstimateOpt(&segments[0], &segments[1], matrix[0][1], matrix[1][0], AUTOREFINE_PRECOND_BLOCK);
        m_uiBlockPreNum = 2;
        m_ulBlockPreBaseNum = 0;
        ComputeBlockPrecond();
        int failures = 0;
        for(int i=0; i<2; ++i) {
            for(int j=0; j<2; ++j) {
                double identity = matrix[i][0]*m_pdBlockPrecond[0][j] + matrix[i][1]*m_pdBlockPrecond[1][j];
                if(!std::isfinite(identity) || std::fabs(identity-(i == j ? 1.0 : 0.0)) > 1e-12) {
                    std::printf("FAILED 2D block inverse (%d,%d): %.17g\n", i, j, identity);
                    ++failures;
                }
            }
        }
        return failures;
    }
};

int main()
{
    g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
    CBlockProbe probe;
    int failures = probe.CheckApplication();
    failures += probe.CheckMaterial();
    failures += probe.CheckSegments();
    return failures ? 1 : 0;
}
