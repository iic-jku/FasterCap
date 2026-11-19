# Dielectric flux integration

With `-g`, finite triangular dielectric interactions integrate the target field analytically
and use the existing distance-dependent quadrature over the source. The target-normal
integral is an oriented solid angle. Tangential components use boundary line integrals,
so averaged superpanel normals remain supported. Both panel and vertex entry points share
this calculation.

For a closed target surface with outward normals, the total flux from a unit source is
`1/eps0` inside, zero outside and the principal value `1/(2*eps0)` on a smooth face.
These identities hold for polyhedra: spherical geometry approximation cannot explain a
large violation. Integrating the target analytically preserves the identities at each
source quadrature point when the same source quadrature is used for every target face.
Exactly coplanar source triangles use the principal value before quadrature rounding can
move their points to either side of the field jump.

Coplanarity uses a filtered determinant with an exact fallback for the existing binary32
vertex coordinates. Error-free product tails and expansion summation distinguish a zero
determinant from a small nonzero one, without a distance cutoff. An ordinary rounded dot
product is insufficient: even a triangle tested against itself can otherwise acquire a
false field jump. This fallback assumes the input vertices remain binary32; changing the
vertex representation requires rechecking its arithmetic contract.

The former near interaction subtracted potentials displaced by `h = distance/20` and
sampled their derivative on the target. Near an interface this smeared the field jump,
and target quadrature introduced a separate flux error. Reducing `h` did not eliminate
that second error. The retained `h` argument of `CPotential::MutualDHalfNumerical()` is
ignored for source compatibility.

This changes numerical results in the 3D Galerkin scheme. It does not enforce symmetry on
the output matrix. Source quadrature, distance-dependent rule changes, far monopoles,
hierarchical approximation, mesh geometry and the linear solve still contribute error.
The complete operator is therefore not promised to conserve flux to rounding or produce
a universally symmetric matrix. The default collocation scheme and the 2D solver retain
their existing formulas.

The bounded regressions use twelve target triangles for exact flux identities and fixed
nested-sphere meshes of 960 or 3840 panels for physical capacitance. They cover reversed
normals, scaling, translation, skewed faces, arbitrary field directions, lossy permittivity,
offset interfaces, thin coats and permittivity contrast 100. Coplanar principal values and
the one-sided jump across a representable gap of `2^-60` are checked separately. The
accuracy tests require no adaptive growth or out-of-core storage. Local consistency checks
also use 1280 panels and 1,637,120 links to compare threads and storage modes.

Fixed meshes use `-m1e32`; `-d` must then be reduced accordingly because the interaction
tolerance is the product of `-m` and `-d`. Tight linear-solver tolerance alone does not
control source quadrature, hierarchy or collocation error. Neither this integration change
nor a low reciprocity defect establishes convergence to the continuum solution.

Analytical integration is useful near panels, while numerical cubature remains useful for
far interactions; see [Glück and Hilk, 2016](https://arxiv.org/abs/1606.03743). The local
experiments and their limitations are recorded separately in project memory. The filtered
exact predicate follows the expansion-arithmetic approach described by
[Shewchuk, 1997](https://www.cs.cmu.edu/~quake/robust.html).
