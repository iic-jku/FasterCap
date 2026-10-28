# Project conventions

## Repositories

- This is `iic-jku/FasterCap` (remote `jku`, formerly `martinjankoehler/FasterCap`), the
  maintained fork of `ediloren/FasterCap` (remote `upstream`, 6.0.7, `b42179a`, dormant).
  Its issues are written `iic-jku/FasterCap#N` outside this repository.
- CMake builds the sibling checkouts `../LinAlgebra` and `../Geometry` (`add_subdirectory`),
  and the sources include `LinAlgebra/…` and `Geometry/…` from `..`. CI checks out the
  `master` of `iic-jku/LinAlgebra` and `iic-jku/Geometry`, forks of the ediloren
  repositories with their own AGENTS.md.
- Locally (`~/Source/FastFieldSolvers`), the working checkouts are `FasterCap_MJK`,
  `LinAlgebra_MJK` and `Geometry`; `FasterCap`, `LinAlgebra`, `LinAlgebra.orig` and
  `LinAlgebra_orig` are references only. Build from `work/`, whose symlinks `FasterCap`,
  `LinAlgebra` and `Geometry` point to the working checkouts: a build from `FasterCap_MJK`
  itself would use the reference `LinAlgebra`.
- Remotes use SSH URLs (`git@github.com:…`).
- The repository stays in ediloren's fork network while PRs to upstream are open, so GitHub
  proposes `ediloren/FasterCap` as the base of a new PR: pick `iic-jku/FasterCap` `master`
  (`gh pr create --repo iic-jku/FasterCap --base master`).
- `klayout-pex` uses the released `.deb`/`.tar.gz` assets: a change of the numbers or of the
  output format is a behavior change for it.

## Branches

- `master` is `upstream/master` plus merge commits (`--no-ff`) of feature branches only;
  `master-backup-<date>` branches keep earlier histories. Feature branches are never deleted.
- Changes offered upstream (`ediloren/FasterCap#4`–`#7`: `wxwidgets-update`,
  `cmake-openmp-fix`, `out-of-core-unix-support`, `blas-support`) and other upstreamable
  ones (`hot-loop-optimizations`, `pdf-doc`) live on branches based on `upstream/master`.
  Fix them with new commits on the branch, never by rewriting it (the PRs show the history),
  then merge the branch into `master` again.
- Fork-only changes (CI, packaging, AGENTS.md) go on branches based on `master`, e.g.
  `ci-and-agents`.
- Tests go on `regression-tests`, based on `master`: every fix gets a test there that fails
  without it. `master` merges `regression-tests` together with the fixes.

## Scope

- One concern per branch; report unrelated findings instead of fixing them.
- Say in the commit body when behavior changes: command line options, output format, exit
  and return codes (`FC_*`), numerical results.
- Fix root causes in the shared code (`Solver/`, `LinAlgebra`, `Geometry`), not with special
  cases. Platform `#ifdef`s only where the OS API differs, each with a fallback.

## Domain

- Results must make physical sense: the capacitance matrix is symmetric within the solver
  tolerance, its diagonal is positive, the off-diagonal terms are negative, and row sums are
  not negative.
- In-core and out-of-core runs give identical matrices. Force out-of-core with `-f<ratio>`
  on an input with more than one link block (more than 2^20 links, see
  `AUTOREFINE_LINK_CHUNK_SIZE`); check both after every change to link storage, loading or
  multiplication.
- Results differ slightly between platforms, compilers, BLAS libraries and thread counts:
  compare against analytic values, the in-core run, or an upstream build on the same
  machine, never against numbers from another platform.

## Code

- Keep each file's line endings (most sources are CRLF), indentation (tabs) and naming
  (`m_ul…`, `m_d…`, `C…` classes).
- In hot loops (`CMultHier::ComputePanelPotentials_2fast()`, `CAutoRefine::AutoRefineLinks()`)
  hoist only values that are invariant: `LoadLinks()` changes `m_ulCurrBlock` inside the
  block loop.
- `ASSERT` is `wxASSERT`, which wxWidgets 3 also evaluates in release builds: keep it out of
  hot loops.
- Errors return the documented error value and take the existing error path (e.g. `-1` from
  `portableGetFreeMemory()` means "unknown", so the run continues in-core), never `0` or a
  partial result.
- CMake: use imported targets (`OpenMP::OpenMP_CXX`, `BLAS::BLAS`), no hard-coded `-l`
  flags and no paths outside the source tree; BLAS detection belongs in `LinAlgebra`.
- Docs describe contracts, not file layout.
- Don't duplicate code, build logic included: extend the shared code rather than adding a
  variant.

## Build and tests

- macOS: `export OpenMP_ROOT=$(brew --prefix)/opt/libomp`, then
  `cmake -S ../work/FasterCap -B ../work/build -DCMAKE_BUILD_TYPE=RelWithDebInfo
  -DFASTFIELDSOLVERS_HEADLESS=ON
  -DwxWidgets_CONFIG_EXECUTABLE=$(brew --prefix)/bin/wx-config` and
  `cmake --build ../work/build -j`. Linux: the same with `/usr/bin/wx-config`.
- Spell build types exactly (`Release`, `RelWithDebInfo`): an unknown `CMAKE_BUILD_TYPE`
  silently builds without `-O`. Benchmarks compare builds of the same type.
- Run `ctest --test-dir ../work/build --output-on-failure` before calling work done; report
  failures with their output. The tests (`tests/`, Python 3 standard library only) generate
  their geometries in `tests/geometries.py` and check each matrix with
  `tests/check_capacitance.py`.
- Each bug fix gets a test that fails without it. A test that needs a workaround means the
  code is wrong: fix the code.
- Prefer geometries with analytic values or symmetries to stored numbers. Each tolerance
  comes with its reason in `tests/CMakeLists.txt`.
- FasterCap doesn't read PATRAN `.neu` files: FastCap2 examples using them, like
  `1x1bus.lst`, give `nan`. Don't use them as a reference.
- No loosened tolerances or re-baselined numbers without a physical explanation.
- Packaging changes: install the `.deb` in a clean `ubuntu:22.04` and `ubuntu:24.04`
  container and run `FasterCap -bv` there.

## Git

- Martin reviews in detail and commits himself: leave changes uncommitted in the working
  copy. For a change that should become several commits, stage it and propose the split.
  Proposed commit messages are welcome as a starting point. Don't commit or push unless
  asked; never commit scratch files.
- Martin keeps WIP stashes: never chain `git stash push`/`pop`. Test staged changes via
  `git checkout-index -a --prefix=<dir>/` or a temporary worktree.
- Commit subject: `Area: verb phrase`, 6–12 words, no articles. Area is a file, class or
  function, e.g. `ci-cd.yml: include distribution version in asset package filenames`,
  `CMultHier::ComputePanelPotentials_2fast(): remove assertion from hotspot`. Usually no
  body; one or two lines for a non-obvious why or a behavior change.
- Commit and PR texts contain only what the reader can't get from the issue, diff or CI:
  the non-obvious why, rejected alternatives, behavior changes, risks. PR: issue reference
  plus 1–4 sentences, no Summary/Changes/Testing headings.
- Propose a PR title with every PR text, in the style of the commit subjects; for several
  commits, name what they have in common rather than repeating the first subject.
- Every PR for a bug or feature references an issue: `Fixes #N` or `Part of #N`; without an
  issue, file one first. CI, tooling and developer-docs PRs need none. Upstream issues are
  written `ediloren/FasterCap#N`.
- A PR is ready when its CI is green on its final commits: after each push, check the CI
  and report every red job with its cause.
- After rewriting a pushed branch, say that it must be pushed before the merge.
- Merge only green PRs: `master` stays green, as a red baseline hides new failures.
- Issues show the smallest input that reproduces the problem and end with the versions
  used: FasterCap version and commit, OS and CPU architecture, compiler, wxWidgets, BLAS
  library, build type, `OMP_NUM_THREADS`, e.g. `FasterCap 6.0.9 (1955bf7), macOS 26.6.2
  arm64, Apple clang 21, wxWidgets 3.3.3, Accelerate, RelWithDebInfo, 11 threads`.
- Aim for the sweet spot, not minimal length: too terse makes the reader guess the why;
  repetition and self-justification make the reader parse and check sentences that add
  nothing.
- In commit bodies, PR texts, issues, comments and docs, prefer bullet lists to long
  sentences chained with semicolons: one fact per bullet, under a short lead-in such as
  `New options:`. No period at the end of a bullet.
- Mark identifiers as code with backticks: files, classes, functions, macros, CMake
  variables, command line tools and their options, e.g. `Solver/Autorefine.cpp`,
  `CAutoRefine::LoadLinks()`, `CMAKE_BUILD_TYPE`, `FasterCap -b -f1500`. Tick every mention,
  not only the first.
- No AI attribution: no `Co-Authored-By` trailer, no "Generated with …" line.
- Before handing over, re-read your diff against `master` and check it against this file.
