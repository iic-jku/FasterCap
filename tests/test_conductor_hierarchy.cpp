// Physical conductor descendants must determine near-field link selection.
#include "MultiplyHierarchical.h"
#include <cmath>
#include <cstdio>

class CConductorProbe : public CMultHier
{
public:
	int Interact(CAutoPanel *group, CAutoPanel *single, CLin_Vector &charge,
		CLin_Vector &result, unsigned long &count)
	{
		CAutoConductor grouped, isolated;
		grouped.m_uTopElement.m_pTopPanel = group;
		isolated.m_uTopElement.m_pTopPanel = single;
		m_stlConductors.push_back(&grouped);
		m_stlConductors.push_back(&isolated);
		CAutoElement *nodes[] = {group, group->GetLeftChild(), group->GetRightChild(), single};
		CAutoElement *storage[4] = {}, **links = storage;
		double values[4] = {}, *coefficients = values;
		m_clsGlobalVars.m_cScheme = AUTOREFINE_COLLOCATION;
		m_clsGlobalVars.m_dEps = 1e10; // Isolate geometric admissibility.
		m_dMaxSide = m_dMaxArea = 1.0;
		m_bComputeLinks = false;
		m_iLevel = m_iMaxLevel = 0;
		m_ulLinksNum[0] = 0;
		for(CAutoElement *node : nodes)
			node->m_ulLinkIndexStart[0] = node->m_ulLinkIndexEnd[0] = 0;
		group->CalcExpansionBounds();
		single->CalcExpansionBounds();
		int status = RefineMutual(group, single);
		count = m_ulLinksNum[0];
		unsigned long start = 0;
		for(CAutoElement *node : nodes) {
			unsigned long nodeCount = node->m_ulLinkIndexEnd[0];
			node->m_ulLinkIndexStart[0] = node->m_ulLinkIndexEnd[0] = start;
			start += nodeCount;
		}
		m_pNodes = nodes;
		m_pdPanelPtrLinks[0] = &links;
		m_dPotCoeffLinks[0] = &coefficients;
		m_ulNodeNum[0] = 4;
		m_ulLinkChunkNum[0] = m_ulBlocksNum = 1;
		m_ulCurrBlock = 0;
		m_bComputeLinks = true;
		if(status == FC_NORMAL_END && count <= 4) status = RefineMutual(group, single);
		else status = FC_INVALID_RESULT;
		if(status == FC_NORMAL_END) status = ComputeLinkCoefficients(0, 4, 0);
		m_clsSelfPotCoeff.newsize(3);
		m_clsSelfPotCoeff = 0.0;
		m_clsChargeVect = charge;
		m_clsPotVect = result;
		if(status == FC_NORMAL_END) {
			ComputePanelCharges_fast();
			status = ComputePanelPotentials_2fast();
			ComputeLeafPotentials_fast();
		}
		// Nodes, conductors and the deliberately tiny link block are stack-owned.
		m_pNodes = nullptr;
		m_pdPanelPtrLinks[0] = nullptr;
		m_dPotCoeffLinks[0] = nullptr;
		m_ulNodeNum[0] = m_ulLinksNum[0] = m_ulLinkChunkNum[0] = 0;
		m_stlConductors.clear();
		return status;
	}
};

int main()
{
	g_bFCContinue = true;
	g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
	CConductorProbe probe;
	int failures = 0;
	for(double distance : {3.0, 64.0}) {
		CAutoPanel patches[2], group, single;
		CAutoPanel *leaves[] = {&patches[0], &patches[1], &single};
		for(int i=0; i<3; ++i) {
			double x = i==2 ? distance : i==0 ? -4.0 : 4.0;
			leaves[i]->m_clsVertex[0].pos(x,0,0);
			leaves[i]->m_clsVertex[1].pos(x,0.125,0);
			leaves[i]->m_clsVertex[2].pos(x,0,0.125);
			leaves[i]->CalcPanelGeomPar();
		}
		group.MakeSuperPanel(&patches[0], &patches[1]);
		const double charges[][3] = {{1,1,0}, {1,-1,0}, {0,0,1}};
		for(const auto &values : charges) {
			CLin_Vector charge(3), result(3);
			for(int i=0; i<3; ++i) charge[i] = values[i];
			unsigned long links;
			if(probe.Interact(&group, &single, charge, result, links)) return 2;
			unsigned long expectedLinks = distance==3.0 ? 4 : 2;
			if(links != expectedLinks) {
				std::printf("FAILED conductor extent at %g: %lu links, expected %lu\n",
					distance, links, expectedLinks);
				++failures;
			}
			// Near interactions must agree with their physical leaf kernels. Both
			// source grouping (including neutral charge) and target grouping matter.
			if(distance==3.0) {
				double expected[3] = {};
				for(int i=0; i<2; ++i) {
					double forward, reverse;
					if(probe.PotEstimateOpt(&single, &patches[i], forward) ||
						probe.PotEstimateOpt(&patches[i], &single, reverse)) return 2;
					expected[2] += charge[i]*forward;
					expected[i] += charge[2]*reverse;
				}
				for(int i=0; i<3; ++i) {
					double error = std::fabs(result[i]-expected[i])*FOUR_PI_TIMES_E0;
					if(!std::isfinite(error) || error > 1e-12) {
						std::printf("FAILED conductor q=(%g,%g,%g) leaf %d: normalized error %.17g\n",
							values[0], values[1], values[2], i, error);
						++failures;
					}
				}
			}
		}
	}
	return failures ? 1 : 0;
}
