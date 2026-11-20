// A uniform electric field must project onto each physical leaf normal,
// independently of how those panels are grouped or of a previous product.
#include "MultiplyHierarchical.h"
#include <cmath>
#include <cstdio>

class CFieldProbe : public CMultHier
{
public:
    int Interact(CAutoPanel *source, CAutoPanel *target, CLin_Vector &charge, CLin_Vector &result)
    {
        CAutoConductor sources, targets;
        sources.m_uTopElement.m_pTopPanel = source;
        targets.m_uTopElement.m_pTopPanel = target;
        m_stlConductors.push_back(&sources);
        m_stlConductors.push_back(&targets);
        m_clsSelfPotCoeff.newsize(charge.size());
        m_clsSelfPotCoeff = 0.0;
        m_clsChargeVect = charge;
        m_clsPotVect = result;
        CAutoElement *node = target, *link = source;
        CAutoElement **links = &link;
        C3DVector distance = C3DVector(target->GetCentroid()) - C3DVector(source->GetCentroid());
        double r = Mod(distance), coefficient = 1.0 / (FOUR_PI_TIMES_E0*r*r*r);
        double *coefficients = &coefficient;
        m_pNodes = &node;
        m_pdPanelPtrLinks[0] = &links;
        m_dPotCoeffLinks[0] = &coefficients;
        m_ulNodeNum[0] = m_ulLinksNum[0] = m_ulLinkChunkNum[0] = m_ulBlocksNum = 1;
        target->m_ulLinkIndexStart[0] = 0;
        target->m_ulLinkIndexEnd[0] = 1;
        ComputePanelCharges_fast();
        int status = ComputePanelPotentials_2fast();
        ComputeLeafPotentials_fast();
        // The artificial one-link hierarchy and its panels are stack-owned.
        m_pNodes = nullptr;
        m_pdPanelPtrLinks[0] = nullptr;
        m_dPotCoeffLinks[0] = nullptr;
        m_ulNodeNum[0] = m_ulLinksNum[0] = m_ulLinkChunkNum[0] = 0;
        m_stlConductors.clear();
        return status;
    }

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
        root->m_clsFieldExpansion.m_clsField = field;
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
    // Harmonic potential x^3-3xy^2 gives field (3x^2-3y^2,-6xy,0).
    // Its mean x component on this triangle is 1; centroid sampling gives 5/3.
    // Translate through two arbitrary centers to exercise the downward shifts.
    CDielectricFieldExpansion polynomial, intermediate, leaf;
    polynomial.Clear();
    intermediate.Clear();
    leaf.Clear();
    polynomial.m_dHessian[DielectricExpansion::Index(0,0,0)] = 6;
    polynomial.m_dHessian[DielectricExpansion::Index(0,1,1)] = -6;
    C3DVector middle(0.25,-0.5,0.75), center(1,2.0/3.0,2.0/3.0);
    intermediate.AddParent(polynomial, middle);
    leaf.AddParent(intermediate, center-middle);
    C3DVector_float triangle[3] = {C3DVector_float(1,0,0),
                                  C3DVector_float(1,2,0), C3DVector_float(1,0,2)};
    double mean = leaf.MeanNormalField(triangle, center, C3DVector(1,0,0));
    if(std::fabs(mean-1.0) > 1e-12) {
        std::printf("FAILED quadratic triangle flux: %.17g expected 1\n", mean);
        ++failures;
    }
    // Stored float centroids are rounded, even when the vertices are exact.
    // Integration must include the first moment about that stored center.
    for(double offset : {0.0, 1048576.0}) {
        CAutoPanel translated;
        for(int i=0; i<3; ++i) {
            translated.m_clsVertex[i].pos(offset+triangle[i].x,
                                          offset+triangle[i].y, offset+triangle[i].z);
        }
        translated.CalcPanelGeomPar();
        C3DVector storedCenter = translated.GetCentroid();
        leaf.Clear();
        leaf.AddParent(polynomial, storedCenter-C3DVector(offset,offset,offset));
        mean = leaf.MeanNormalField(translated.m_clsVertex, storedCenter, C3DVector(1,0,0));
        if(std::fabs(mean-1.0) > 1e-12) {
            std::printf("FAILED rounded-center triangle flux at %g: %.17g expected 1\n", offset, mean);
            ++failures;
        }
    }
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
    // Independent oracle: exact triangle flux from the two physical source
    // centroids. Halving both clusters' size/distance ratio must reduce the
    // truncation error, including a neutral source whose monopole is zero.
    for(double secondCharge : {-0.7, -1.0}) {
        double previousError = 1.0;
        for(double distance : {8.0, 16.0, 32.0}) {
            CAutoPanel sources[2], source, targets[2], target;
            for(int i=0; i<2; ++i) {
                double y = i ? 0.5 : -0.5;
                sources[i].m_clsVertex[0].pos(0,y,0);
                sources[i].m_clsVertex[1].pos(0,y+0.25,0);
                sources[i].m_clsVertex[2].pos(0,y,0.25);
                sources[i].CalcPanelGeomPar();
                targets[i].m_clsVertex[0].pos(distance,y,-0.25);
                targets[i].m_clsVertex[1].pos(distance+0.5,y+0.5,0);
                targets[i].m_clsVertex[2].pos(distance,y,0.5);
                targets[i].CalcPanelGeomPar();
                targets[i].m_ucType |= AUTOPANEL_IS_DIEL | AUTOPANEL_OUTPERM_NORMAL_DIR;
            }
            source.MakeSuperPanel(&sources[0], &sources[1]);
            target.MakeSuperPanel(&targets[0], &targets[1]);
            CLin_Vector charge(4), result(4);
            charge = 0.0;
            charge[0] = 1.0;
            charge[1] = secondCharge;
            if(probe.Interact(&source, &target, charge, result)) return 1;
            double error = 0.0, scale = 0.0;
            for(int i=0; i<2; ++i) {
                double expected = 0.0;
                for(int j=0; j<2; ++j) {
                    double coefficient;
                    probe.PotEstimateOpt(&targets[i], &sources[j], coefficient);
                    expected += charge[j] * coefficient;
                }
                error = std::fmax(error, std::fabs(result[2+i] - expected));
                scale = std::fmax(scale, std::fabs(expected));
            }
            error /= scale;
            if(!std::isfinite(error) || error > previousError*0.3) {
                std::printf("FAILED coupled field q2=%g distance=%g error=%g previous=%g\n",
                            secondCharge, distance, error, previousError);
                ++failures;
            }
            previousError = error;
            charge = 0.0;
            if(probe.Interact(&source, &target, charge, result)) return 1;
            if(result[2] != 0.0 || result[3] != 0.0) {
                std::printf("FAILED field moments retained across products\n");
                ++failures;
            }
        }
    }
    return failures ? 1 : 0;
}
