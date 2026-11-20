# Bounded performance comparisons

`FasterCap -i` adds wall timings for upward charge/moment gathering, linked
interactions, downward field propagation and out-of-core link loading. Link time
includes loading; subtract the separately reported loading time to estimate compute
time. Counts and accumulators reset for each refined solve. The timings include
hierarchical products used by preconditioning as well as GMRES and apply to each
real product in a complex solve. Without `-i`, these phase timers do not read the clock.

Existing timing labels and the numerical JSON schema are unchanged. Input,
hierarchy-building, discretization and preconditioner timings use CPU time; matrix
construction, solve time and the new hierarchy timings use elapsed wall time.
Do not add CPU timings to wall timings or add nested timings twice.

For a headless build, run the benchmark through the resource guard, using one shared
scratch directory and one build:

```sh
python3 tools/run_guarded.py --root /path/to/scratch --seconds 180 -- \
  python3 tools/benchmark_solver.py \
  --fastercap /path/to/build/FasterCap \
  --counter /path/to/build/tests/count_benchmark_links \
  --workdir /path/to/scratch --output /path/to/scratch/baseline.json
```

The default runs small conductor, dielectric, uneven and lossy fixed meshes with
1/2/4 OpenMP threads and one BLAS worker, interleaving three repetitions. It requires
identical full-precision matrices, panel/link counts, printed GMRES step counts and
hierarchical product counts across runs. Potential-estimate diagnostics are recorded
separately because older builds have a shared-counter race. Numerical quality and
analytic errors are recorded, not converted into machine-dependent speed thresholds.
This benchmark does not replace the physical regression suite.

`--size medium` and `--size large` select depth-3 and depth-4 sphere meshes; the
uneven fixture always retains its small depth-2 mesh. `--cases` selects fixtures,
and `--repeats 1..5` controls repetitions. The helper counts interactions through
the production refiner, without allocating links, and rejects more than 16,000
panels or eight million links. Each subprocess has a 30-second limit. The fixed
`-m1e32 -d1e-34 -t1e-9 -pj -f0` settings isolate thread scheduling from adaptation.
They are diagnostic settings, not general accuracy recommendations.

Generated geometry and temporary reports are removed after the run. The summary
records build provenance, executable hash, wall times, solver work, matrix quality
and solver-reported allocated memory. The outer guard receipt supplies sampled peak
process-group RSS and resource-limit outcomes. Solver allocation accounting is not
RSS, and repeated timings are not guarantees for another machine.

Use matching compiler flags, build types and dependencies when comparing branches.
Keep numerical acceptance unchanged during performance work and separately verify
genuine multi-block in-core/out-of-core equality for changes to multiplication or
link handling. The timer regression checks behavior and accounting, never speed.
