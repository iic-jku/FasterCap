# Parallel coefficient diagnostics

Coefficient workers compute values without changing shared estimate counters or
warning flags. Successful link construction counts each directed coefficient once,
including dielectric radial coefficients, in addition to existing serial estimates.
Counts therefore become independent of worker count and storage-block partitioning.
Older builds could undercount through races and omit grouped dielectric coefficients.

Finite near-distance warnings remain visible. Workers record the lowest-index
affected pair in each loaded block, and the serial caller reports it once without
recomputing the coefficient. The number of warning messages can change with block
partitioning; warnings do not affect numerical results or estimate counts.

A failed or nonfinite coefficient returns `FC_INVALID_RESULT` before that block is
saved or used by the solver. The failing value is not published, and no failed
calculation is repeated merely to print diagnostics. Existing numerical kernels and
valid coefficient values are unchanged.
