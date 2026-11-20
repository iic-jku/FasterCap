# Controlled platform comparison

The original plate from issue #10 is captured on each supported CI build with
`tests/check_platform_comparison.py`. The generated fixture preserves commit
`9d7d5b8`'s panel ordering, materials and reference point. Its coordinate numbers
are micrometers; FasterCap interprets them as meters, so multiply the reported
capacitances by `1e-6` for the physical values in farads. No analytic plate
capacitance or stored reference matrix is asserted.

The capture separates three questions:

- `fixed_mesh` divides the original planar faces into strips of aspect at most ten,
  yielding 1,404 explicit triangles and at most 1,969,812 directed links. It has
  the same physical surfaces but a deliberately
  different discretization from the quadrilateral input. The initially investigated
  36-triangle mesh was rejected because its thin panels gave wrong matrix signs. It uses
  `-m1e32 -d1e-50 -t1e-10 -pj -f0`. The geometry hash fixes the input mesh. Even this tiny positive interaction
  threshold can retain grouped dielectric interactions whose refinement estimate
  is zero, so this is a tight-hierarchy comparison, not an all-direct operator.
  Link counts and all effective settings must match before its tight cross-build
  comparison can pass. On the f7a957b macOS capture this coarse mesh still has
  about 4.6% reciprocal-pair disagreement and a negative row sum; it is arithmetic
  comparison evidence, not a physically converged capacitance result
- `adaptive` uses the original `-a0.05 -d0.5 -m0.5 -ap` options, in-core. Every
  completed matrix, effective setting and panel/link count is retained. Quadrilateral
  input processing subdivides skinny faces before mesh refinement: on the f7a957b
  macOS build, even `-m1e32` starts with 2,764 triangles, not 36
- `manual_replay` repeats the final reported mesh and interaction settings with
  `-t1e-10 -pj`, recording its count match and matrix difference from the automatic
  result. Passing `-t` alongside `-a` would not perform this experiment: automatic
  mode overrides the GMRES tolerance. Matching counts and settings is evidence of
  matching work, not a proof of identical hierarchical topology. Unknown automatic
  settings or a replay count mismatch fail explicitly; algorithms that additionally
  enrich the mesh require replay from their exact exported refined geometry

All captures force one OpenMP and one BLAS worker through their environment.
The existing `solver_repeatability` regression separately checks exact matrices
within a build across repeated one-, two- and four-thread runs. Cross-build
comparison uses full-precision JSON reports and all platform pairs, never rounded
console values. It rejects missing/failed results, dirty or differing repository
revisions, per-discretization input hashes, options, matrix shapes and conductor
ordering. The report
retains actual compiler, configuration, wxWidgets, linked BLAS, source revisions
and runtime host/package metadata. In particular Windows currently builds
`Release`, while the Unix configurations use `RelWithDebInfo`.

The predeclared investigation budgets are relative to the larger magnitude of each
corresponding complex entry:

- `1e-7` for the fixed-mesh solve and manual replays with matching reported work,
  allowing amplification above the `1e-10` algebraic target while detecting the
  historical percent-level discrepancy
- `0.005` for automatic results or manual replays with different reported work,
  one tenth of the original `-a0.05` target and below the historical roughly 3%
  divergence

These are regression budgets, not conditioning-based error estimates or physical
accuracy guarantees. Every comparison records actual entrywise and Frobenius
changes, matching work/path flags and iteration counts. A failure requires an
investigation; changing the tolerance to admit a new result does not explain it.
Physical mesh convergence is tracked separately in issue #78.

Each solver process has a 60-second timeout. Each report/log is capped at 2 MiB
with a failed capture if it grows beyond that limit. Reports above 16,000 panels or eight
million links are rejected. These report checks are post-run limits, not memory
allocation guards; local probes must also use `tools/run_guarded.py` and perform
count-only preflight where applicable. CI retains only small reports, logs and
geometry files in `platform-report-*` artifacts. It does not duplicate binaries.

For a local fixed-mesh first probe, use `capture --fixed-only --fastercap ...
--output-dir ...`. This intentionally incomplete capture cannot satisfy the
four-platform comparison. For complete captures, omit `--fixed-only`, then use
`compare <capture.json> ... --output=comparison.json`. A controlled dirty local
comparison requires explicit `--allow-dirty`; CI does not enable it.
