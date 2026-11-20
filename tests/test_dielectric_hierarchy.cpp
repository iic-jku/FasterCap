// A uniform electric field must project onto each physical leaf normal,
// independently of how those panels are grouped or of a previous product.
#include "MultiplyHierarchical.h"
#include <cmath>
#include <cstdio>

class CFieldProbe : public CMultHier
{
public:
    void Project(CAutoPanel *root, const C3DVector &field, CLin_Vector &result)
    {
        CAutoConductor conductor;
        conductor.m_uTopElement.m_pTopPanel = root;
        m_stlConductors.push_back(&conductor);
        CLin_Vector charge(result.size());
        charge = 0.0;
        m_clsSelfPotCoeff.newsize(result.size());
        m_clsSelfPotCoeff = 0.0;
        m_clsChargeVect = charge;
        m_clsPotVect = result;
        // This also tests clearing every field accumulator between products.
        ComputePanelCharges_fast();
        root->m_clsField = field;
        ComputeLeafPotentials_fast();
        m_stlConductors.clear(); // All tree elements and this conductor are stack-owned.
    }
};

int main()
{
    g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
    CFieldProbe probe;
    probe.m_clsGlobalVars.m_cScheme = AUTOREFINE_COLLOCATION;
    int failures = 0;
    for(int grouping=0; grouping<2; ++grouping) {
        CAutoPanel leaves[4], left, right, root;
        const double vertices[4][3][3] = {
            {{0,0,0},{0,1,0},{0,0,1}},
            {{2,0,0},{2,0,1},{2,1,0}},
            {{0,0,2},{0,0,3},{1,0,2}},
            {{0,2,0},{1,2,0},{0,3,0}}
        };
        for(int i=0; i<4; ++i) {
            for(int k=0; k<3; ++k)
                leaves[i].m_clsVertex[k].pos(vertices[i][k][0], vertices[i][k][1], vertices[i][k][2]);
            leaves[i].CalcPanelGeomPar();
            leaves[i].m_ucType |= AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR;
        }
        // Grouping 0 includes two opposing normals: its mean normal is zero.
        left.MakeSuperPanel(&leaves[0], &leaves[grouping ? 2 : 1]);
        right.MakeSuperPanel(&leaves[grouping ? 1 : 2], &leaves[3]);
        root.MakeSuperPanel(&left, &right);
        const C3DVector fields[] = {C3DVector(1.25,-2,3.5), C3DVector(-3,0.5,2), C3DVector(0,0,0)};
        for(const C3DVector &field : fields) {
            CLin_Vector result(4);
            probe.Project(&root, field, result);
            for(int i=0; i<4; ++i) {
                double expected = DotProd(field, leaves[i].GetDielNormal());
                double obtained = leaves[i].m_dPotential;
                if(!std::isfinite(obtained) || std::fabs(obtained-expected) > 1e-12) {
                    std::printf("FAILED group %d leaf %d: %.17g expected %.17g\n",
                                grouping, i, obtained, expected);
                    ++failures;
                }
            }
        }
    }
    return failures ? 1 : 0;
}
