# Dielectric collocation and hierarchical fields

The experimental 3D collocation scheme integrates the normal field over each
dielectric target triangle and uses a source-centroid quadrature point. The
two mutual-interaction interfaces use the same calculation. The matrix is not
symmetrized after solving.

Grouped dielectric interactions use a coupled Cartesian expansion through total
degree two in source and target displacement. For the electric-field kernel `f`,
its gradient `G`, and its Hessian `H`, the expansion is:

- Field: `Q*f - G*p + H:S/2`
- Field gradient: `Q*G - H.p`
- Field Hessian: `Q*H`

Here `Q` is total source charge, `p` its first moment, and `S` its symmetric second
moment about the group centroid. Physical source leaves have zero moments because
qualocation represents them by centroid points. Groups gather the moments of their
actual children; the synthetic superpanel triangle contributes no source moment.
A neutral group can still have nonzero moments and must not be skipped.

Each multiplication clears the transient state, gathers source moments, accumulates
linked fields, and translates local expansions down to physical target triangles.
Leaf evaluation integrates the quadratic field over the triangle before projecting
onto its oriented normal. Its covariance about the centroid is
`sum((vertex-centroid)*(vertex-centroid)^T)/12`.

- A parent normal cannot replace the different or opposite normals of its leaves
- Translating only target fields or correcting only source moments omits coupled terms
- A leaf target's scalar link already integrates the source monopole exactly; only source-moment corrections are added there
- Direct leaf-to-leaf interactions, self terms, conductor equations, Galerkin and 2D equations keep their existing kernels
- Real and imaginary products each initialize their own moment and field accumulation

`CAutoRefine::UsesDielectricField()` defines where this representation applies.
Non-leaf dielectric targets store the scalar radial factor
`1 / (4*pi*epsilon0*r^3)`. Link generation and multiplication must agree on that
condition. `CAutoRefine::PotEstimateOpt()` continues to return scalar estimates for
preconditioners and direct interactions.

The transient state stores 28 doubles per 3D panel, including internal tree nodes:
224 bytes, or 200 bytes more than the earlier constant-vector implementation.
Symmetric tensors store only their distinct components. Source moments are gathered
only for 3D collocation models containing a dielectric. Each directed link still
stores one double, and the out-of-core block layout is unchanged. One worker owns
each node's link sum; downward propagation occurs after every block has contributed.

The bounded regressions cover uniform fields, exact integration of a known harmonic
quadratic, far-field convergence with neutral and nonneutral sources, repeated
products, and fixed-mesh moderate- and high-contrast comparisons against direct
interactions. Source and target moments must be improved together: isolated
higher-order corrections can substantially worsen dielectric charge conservation.

Higher expansion order does not establish an accuracy bound for every accepted
link. Loose grouping, uneven meshes, automatic refinement and very large dielectric
contrast still need investigation. A small linear-system residual or reciprocal
matrix alone does not establish physical accuracy or release readiness.
