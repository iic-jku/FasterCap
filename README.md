# FasterCap

FasterCap is a three- and two-dimensional capacitance extraction program.

This is the fork maintained by the [Institute for Integrated Circuits and Quantum Computing (IICQC)](https://iic.jku.at/) of Johannes Kepler University Linz (IIC-JKU). It continues FasterCap 6.0.7 by Enrico Di Lorenzo, FastFieldSolvers S.R.L., whose repository [ediloren/FasterCap](https://github.com/ediloren/FasterCap) has had no commits since 6.0.7 in 2018. The fork is not a FastFieldSolvers release: report problems in [its issues](https://github.com/iic-jku/FasterCap/issues).

[`klayout-pex`](https://github.com/iic-jku/klayout-pex), the parasitic extraction for `KLayout`, runs `FasterCap` from this fork's releases.

## Changes from 6.0.7

[`History.txt`](History.txt) lists them per release. In short:

- fixes of wrong results: out-of-core runs, dielectric interfaces in 3D and 2D, non-planar `Q` panels, the block preconditioner
- runs that fail say so: exit code 75 when GMRES doesn't converge, 65 for NaN or infinite results, 74 for input `FasterCap` can't read
- out-of-core runs on Linux and macOS, BLAS (OpenBLAS, Accelerate), wxWidgets 3.2 and 3.3, headless builds
- release assets for Ubuntu 22.04 and 24.04, Windows and macOS, built and tested in CI
- regression tests and unit tests, run with `ctest`

## Installing

The [releases](https://github.com/iic-jku/FasterCap/releases) hold console builds without GUI:

- Ubuntu 22.04 and 24.04, x86-64: a `.deb` package, e.g. `sudo apt install ./FasterCap_6.2.0_amd64_ubuntu-24.04.deb`, or a `.tar.gz`
- macOS, arm64: a `.tar.gz` with the executable, which needs Homebrew's `libomp` and `wxwidgets` (`brew install libomp wxwidgets`)
- Windows, x64: a `.tar.gz` with `FasterCap.exe` and the OpenBLAS and wxWidgets DLLs, which needs the Microsoft Visual C++ 2015-2022 Redistributable (x64), also for OpenMP

## Running

```
FasterCap -b <input file> [options]
```

For example, `FasterCap -b cubes.lst -a0.01` refines the geometry until the capacitance matrix changes by less than 1 %. `FasterCap -b?` lists the options, and the help in [`doc/pdf`](doc/pdf) describes the input files, the options and the output.

Exit codes:

| code | meaning |
|---|---|
| 0 | the matrix is printed |
| 1 | other error |
| 64 | wrong command line |
| 65 | the matrix has NaN or infinite entries; it is not printed |
| 66 | the input file can't be opened |
| 71 | out of memory |
| 74 | the input can't be read, e.g. a PATRAN neutral file or a self-intersecting `Q` panel |
| 75 | GMRES didn't converge; the matrix of the last iteration is printed |
| 97 | the run can't go out-of-core |
| 98 | unexpected exception |
| 125 | stopped, e.g. by Ctrl-C or SIGTERM |

## Building

`FasterCap` needs the [`LinAlgebra`](https://github.com/iic-jku/LinAlgebra) and [`Geometry`](https://github.com/iic-jku/Geometry) sources of this fork next to it:

```
git clone https://github.com/iic-jku/FasterCap.git
git clone https://github.com/iic-jku/LinAlgebra.git
git clone https://github.com/iic-jku/Geometry.git
```

Requirements: CMake 3.12 or later, a C++ compiler with OpenMP, a BLAS library and wxWidgets 3.0 or later.

- Ubuntu 24.04: `sudo apt install cmake g++ libomp-dev libopenblas-dev wx3.2-headers libwxgtk3.2-dev` (22.04: `wx3.0-headers libwxbase3.0-dev`)
- macOS: `brew install cmake libomp wxwidgets`, then `export OpenMP_ROOT=$(brew --prefix)/opt/libomp`; macOS uses the Accelerate BLAS
- Windows: Visual Studio with the wxWidgets 3.2 and OpenBLAS binaries; [`.github/workflows/ci-cd.yml`](.github/workflows/ci-cd.yml) shows the steps

A console build without GUI:

```
cmake -S FasterCap -B build -DCMAKE_BUILD_TYPE=Release -DFASTFIELDSOLVERS_HEADLESS=ON
cmake --build build -j
```

On macOS, add `-DwxWidgets_CONFIG_EXECUTABLE=$(brew --prefix)/bin/wx-config` to the first command. Without `-DFASTFIELDSOLVERS_HEADLESS=ON`, the build includes the GUI and needs the wxWidgets core and html libraries as well.

## Testing

```
ctest --test-dir build --output-on-failure
```

The tests in [`tests/`](tests) need Python 3. Most run `FasterCap` on generated geometries and check each matrix: exit code, signs, symmetry, and analytic values where they are known. Unit tests check the potential integrals and the triangulation; they are built in console builds only.

## License

LGPL 2.1 or any later version, see [`LICENCE.txt`](LICENCE.txt).

Copyright (c) 2019 FastFieldSolvers S.R.L., http://www.fastfieldsolvers.com. The fork's changes are under the same license.
