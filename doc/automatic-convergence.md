# Automatic algebraic convergence

Automatic mode checks the charge solve before using its capacitance matrix to decide
whether to refine the mesh. A small residual relative to the excitation is insufficient:
a weak dielectric mode can have a tiny residual and still carry a large charge error.

For each excitation, the solver first runs GMRES with the existing automatic residual
tolerance. It then recomputes `r = b - A*x`, solves `A*delta = r` with a tolerance
relative to this new right-hand side, and updates `x`. Reusing the original excitation
as the normalization would miss a small residual containing a large forward error.

A correction must be small in both measures:

- The Euclidean norm of all charge changes relative to the updated charge vector
- A bound on the change in the capacitance row relative to that row's norm

For the second measure, each panel's correction is multiplied by its surrounding
permittivity. Absolute magnitudes are summed per conductor before taking the row norm,
so opposite panel contributions cannot conceal a large correction. Real and complex
systems use the same conversion as final capacitance extraction. In 2D, the raw
reference-conductor column participates in this check.

Corrections must settle below one quarter of the requested `-a` tolerance, keeping
algebraic changes smaller than the mesh comparison target. The true residual of the
updated solution must also meet the original GMRES target. An exactly zero residual
is accepted directly. At most eight correction solves are allowed per excitation;
exhaustion or an unconverged inner solve stops automatic mode with
`FC_CONVERGENCE_FAILED`. Operator errors, nonfinite arithmetic and cancellation retain
their existing error paths. A completed failed matrix remains identified as failed in
`--report` output and is not used to accept mesh convergence.

The existing Krylov storage and preconditioner are reused. Two additional charge-sized
vectors hold the residual and correction; links and mesh storage do not grow during
these solves. Automatic numerical results, selected mesh, runtime and diagnostic output
can change. Manual `-t` behavior retains its residual-based contract.

These are algebraic convergence estimates, not a rigorous bound on physical capacitance
error. Mesh approximation, hierarchical compression, inaccurate input geometry and
conditioning can still limit the result. Neither reciprocity nor a small residual alone
proves accuracy. The tests distinguish a hidden weak mode from planar conductor shielding
and compare automatic results with a tightly solved version of the identical mesh.
