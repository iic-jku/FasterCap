# Dielectric collocation and hierarchical fields

The experimental 3D collocation scheme integrates the normal field over each
dielectric target triangle and uses a source-centroid quadrature point. The
two mutual-interaction interfaces use the same calculation. The matrix is not
symmetrized after solving.

Grouped dielectric targets carry an electric-field vector. Their linked scalar
coefficient is `1 / (4*pi*epsilon0*r^3)`, where `r` is the distance between target
and source group centroids. Multiplication forms the field by multiplying this
factor by source charge and the source-to-target displacement.

- Each multiplication clears the transient field on every dielectric node
- A parent passes its vector field to its children
- Each leaf projects the accumulated vector onto its own oriented normal
- Leaf interactions and self terms remain scalar contributions to the same equation
- Real and imaginary products each initialize their own field accumulation

Projecting at a parent and copying that scalar to its children loses information:
children can have different or opposite normals, even when their average normal
vanishes. Vector propagation preserves the action of a uniform field independently
of the tree grouping.

`CAutoRefine::UsesDielectricField()` defines where this representation applies.
Link generation and multiplication must use the same condition, including the
non-leaf test for the radial coefficient. `CAutoRefine::PotEstimateOpt()` continues
to return scalar normal-field or potential estimates for preconditioners and
direct interactions.

The transient state adds three doubles per 3D panel. Each directed link still
stores one double, and the out-of-core block layout is unchanged. A node's link
sum is owned by one worker. Downward propagation is sequential, after all blocks
have contributed to the field.

This is a constant-field approximation within each target group and a total-charge
approximation within each source group. Spatial field variation and cancellation
within a source group remain approximation errors. High dielectric contrast can
amplify them substantially. A small linear-system residual or reciprocal matrix
does not establish the accuracy of these approximations.

The bounded regressions require exact uniform-field projection, including opposite
normals and alternate groupings. Moderate-contrast matrix tests compare the same
fixed mesh with grouped and all-direct interactions. They are not a high-contrast
accuracy guarantee or evidence that the default scheme is ready for release.
