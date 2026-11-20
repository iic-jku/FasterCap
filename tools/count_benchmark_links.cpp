// Count the bounded benchmark's fixed-mesh interactions before allocating links.
#include "MultiplyHierarchical.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>

int main(int argc, char **argv)
{
	if(argc != 3) return 64;
	char *end;
	double ratio = std::strtod(argv[2], &end);
	if(*end || end == argv[2] || !std::isfinite(ratio) || ratio <= 0.0) return 64;
	g_bFCContinue = true;
	g_ucSolverType = SOLVERGLOBAL_3DSOLVER;
	CMultHier matrix;
	CAutoRefGlobalVars settings;
	settings.m_sFileIn = argv[1];
	settings.m_dMeshEps = 1e32;
	settings.m_dEps = ratio * settings.m_dMeshEps;
	settings.m_dOutOfCoreRatio = 0;
	if(!std::isfinite(settings.m_dEps) || settings.m_dEps == 0.0) return 64;
	int status = matrix.ReadFastCapFile(&settings);
	if(status != FC_NORMAL_END) return status;
	status = matrix.BuildSuperHierarchy();
	if(status != FC_NORMAL_END) return status;
	status = matrix.AutoRefinePanels(settings);
	if(status != FC_NORMAL_END) return status;
	std::printf("BENCHMARK panels %lu nodes %lu links %lu\n",
		matrix.GetPanelNum(0), matrix.GetNodeNum(0), matrix.GetLinksNum(0));
	return matrix.GetPanelNum(0) > 16000 || matrix.GetLinksNum(0) > 8000000 ? 3 : 0;
}
