# Repeatable GMRES scalar products

GMRES orthogonalization uses an ordered compensated dot product. Equal vectors
follow the same reduction order regardless of their buffer alignment or OpenMP
worker count. Four fixed lanes retain independent work without an alignment-dependent
prologue; the lane sums and corrections are combined in a fixed order. A scaled norm remains responsible for avoiding norm underflow and
overflow; the dot-product change does not replace that separate check.

BLAS is still used elsewhere. Its dot-product reduction order can depend on vector
alignment even with one worker, so fixing thread counts alone does not guarantee
repeatable GMRES decisions. Compensated accumulation also retains small terms that
ordinary summation loses during cancellation.

Numerical results and GMRES iteration counts can change from previous builds.
This controls the orthogonalization reduction, not differences between compilers,
geometry kernels, transcendental functions or other BLAS operations. Cross-platform
physical comparisons remain necessary.
