#include "SolveCapacitance.h"
#include <cstdio>

int main()
{
    CSolveCap solver;
    const unsigned char extras[] = {AUTOREFINE_PRECOND_JACOBI, AUTOREFINE_PRECOND_SUPER,
                                    AUTOREFINE_PRECOND_JACOBI | AUTOREFINE_PRECOND_SUPER};
    for(unsigned char extra : extras) {
        CAutoRefGlobalVars settings;
        settings.m_ucPrecondType = AUTOREFINE_PRECOND_BLOCK | extra;
        // No input exists: invalid combinations must be rejected before input I/O.
        if(solver.Run(settings) != FC_COMMAND_LINE_ERROR) {
            std::fprintf(stderr, "shared solver accepted unsupported preconditioners\n");
            return 1;
        }
    }
    return 0;
}
