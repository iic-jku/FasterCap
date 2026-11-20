# Dielectric collocation and hierarchical fields

The 3D collocation scheme integrates the normal field over each
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
onto its oriented normal. The stored float center can differ from the exact triangle
centroid, so both first and second moments about that center contribute. For vertex
displacements `u_i` and `s = sum(u_i)`, these moments are `s/3` and
`(s*s^T + sum(u_i*u_i^T))/12`.

- A parent normal cannot replace the different or opposite normals of its leaves
- Translating only target fields or correcting only source moments omits coupled terms
- A leaf target's scalar link already integrates the source monopole exactly; only source-moment corrections are added there
- Direct dielectric leaf pairs use the target-integrated kernel; self terms, conductor equations, Galerkin and 2D equations keep their existing kernels
- Real and imaginary products each initialize their own moment and field accumulation

`CAutoRefine::UsesDielectricField()` defines where this representation applies.
Non-leaf dielectric targets store the scalar radial factor
`1 / (4*pi*epsilon0*r^3)`. Link generation and multiplication must agree on that
condition. `CAutoRefine::PotEstimateOpt()` continues to return scalar estimates for
preconditioners and direct interactions.

The transient moments and fields store 28 doubles per 3D panel, including internal
tree nodes. Two additional doubles bound the physical source and target extents:
240 bytes in total, or 216 bytes more than the earlier constant-vector implementation.
Symmetric tensors store only their distinct components. Source moments are gathered
only for 3D collocation models containing a dielectric. Physical source leaves skip
moment contractions because their moments are exactly zero. Neutral groups still
contribute their nonzero moments. Explicit symmetric components avoid indexed-vector
branches in the linked-pair kernel. Each directed link still
stores one double, and the out-of-core block layout is unchanged. One worker owns
each node's link sum; downward propagation occurs after every block has contributed.

After mesh refinement, source bounds enclose every descendant leaf centroid and
target bounds enclose every descendant triangle vertex. Child bounds are translated
about the same stored centers used by the expansion. The artificial triangle whose
area represents a superpanel cannot supply these bounds: its area does not constrain
the distance between disconnected or unevenly distributed descendants.

A compressed dielectric interaction requires the sum of its source and target radii
to be at most one quarter of the center distance, in addition to the existing
interaction checks. The retained quadratic expansion omits cubic terms; reducing
the displacement ratio from one half to one quarter reduces their leading size by
a factor of eight and keeps the expansion away from its convergence boundary.
This is a geometric convergence margin, not a bound on the conditioned capacitance
error. Failed checks descend the existing hierarchy; physical leaf pairs keep their
direct kernel. Conductor-only, Galerkin and 2D link selection is unchanged.

The bounded regressions cover uniform fields, exact integration of a known harmonic
quadratic, far-field convergence with neutral and nonneutral sources, repeated
products, and fixed-mesh moderate- and high-contrast comparisons against direct
interactions. Source and target moments must be improved together: isolated
higher-order corrections can substantially worsen dielectric charge conservation.

The correction is enabled by default for 3D dielectric collocation because the old
group-normal representation does not represent its physical leaves. Retaining that
representation as a compatibility option would retain the defect. Numerical results,
selected links, memory use and runtime can change. Tighter grouping can be necessary
for high dielectric contrast; `-d` remains an interaction heuristic, not a relative
capacitance-error guarantee. In particular, loose conductor-potential grouping can
dominate the remaining matrix error even with exact dielectric interactions.

Automatic charge corrections test algebraic convergence separately. Mesh approximation,
hierarchical compression and very large dielectric contrast can still limit physical
accuracy. Compare fixed-mesh grouped and direct solves on bounded models, or analytic
values across successive mesh refinements. A small linear-system residual or reciprocal
matrix alone does not establish physical accuracy.
