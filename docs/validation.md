# Validation of the initial public release

## Reproduce the checks

Run `make test` after building the solver. Python 3 and the same Fortran
compiler/LAPACK/BLAS toolchain used for the executable are required.
The test script uses only Python's standard library and temporary directories.
It compiles the module portion of the original source unchanged to exercise
the public numerical API independently of the command-line input restrictions.

The suite performs four checks:

1. Polarized zero-temperature CLI state, `Nc=1`, at two fields.
2. The same state with `Nc=3`, checking per-spin normalization and periodic
   stacking.
3. Tilted isolated-ion state through the module API, checking energy and
   all dipole/quadrupole components on both sublattices.
4. A deliberately insufficient iteration budget, checking that the program
   terminates unsuccessfully and explicitly reports lack of convergence.

These checks passed using GNU Fortran 14.2.0, LAPACK/BLAS 3.12.1,
Python 3.12.14, and x86-64 Linux. Both optimized builds (`-O2` and `-O3`)
and a separate debug build
(`-O0 -g -std=f2008 -fopenmp -fcheck=all -ffpe-trap=invalid,zero,overflow`)
passed. The tests took approximately one second, including compilation of
the small API driver, in the preparation environment; runtime is hardware
dependent.

The comparison tolerance is an absolute `2e-8` in the units of each checked
quantity. Fixed-point residuals must be below `1e-9`. The tests use one
OpenMP thread and one BLAS thread.

## Polarized-state reference

The supplied `examples/saturation.nml` uses
$J_0=0.1$, $J_1=0.02$, and $J_2=0.01$ meV, $\Delta=0.5$, zero crystal
fields and quadrupolar exchanges, zero temperature, and a field along $z$.

For the fully polarized state the exact energy and moment per spin are

$$
f=J^2(2J_0\Delta+J_1+J_2)-g\mu_B BJ,\qquad
M_z/\mu_B=gJ,\qquad J=\frac52,\quad g=\frac67.
\tag{1}
$$

Each spin contributes two in-plane bonds, one nearest-plane bond, and one
next-nearest-plane bond to the energy when each physical bond is counted
once. This makes the reference sensitive to coordination and mean-field
double counting.

| Field (T) | Analytic energy per spin (meV) | Analytic moment (Bohr magnetons/spin) |
| ---: | ---: | ---: |
| 100 | -11.5911752985714 | 2.14285714285714 |
| 110 | -12.8315428284286 | 2.14285714285714 |

The executable reproduced these energies within approximately `2e-14` meV
per spin in the release preparation runs. Both sublattices had
$J_z=2.5$; transverse dipoles and both quadrupoles vanished to numerical
precision. Self-consistency residuals were below `1e-14`.

These reference values are also provided in
[expected_values.csv](../examples/expected_values.csv).
The test derives the analytic result directly and additionally verifies
the magnetization obtained from the energy's field derivative.

## Tilted isolated-ion reference

For this API-level test all exchanges and crystal fields vanish.
The field is $10$ T, at $45^\circ$ from $z$, with equal $x$ and $y$
components. For a polarized spin coherent state along unit vector
$\mathbf n$, the exact expectation values are

$$
\langle\mathbf J\rangle=J\mathbf n,\qquad
\langle Q_x\rangle=J(J-\tfrac12)n_xn_z,\qquad
\langle Q_y\rangle=J(J-\tfrac12)n_yn_z,\qquad
f=-g\mu_B BJ.
\tag{2}
$$

The reference values are:

| Quantity | Expected value |
| --- | ---: |
| Energy per spin (meV) | -1.24036752985714 |
| Jx | 1.25 |
| Jy | 1.25 |
| Jz | 1.76776695296637 |
| Qx | 1.76776695296637 |
| Qy | 1.76776695296637 |

The complex $y$ component and quadrupole normalization are explicitly
tested on both sublattices.

## Scope and known limits

These tests cover zero-temperature limiting cases, basic continuation,
CSV observables, and reporting of convergence failure. They do not establish
the correctness of every search setting or parameter regime. Finite-temperature
dense diagonalization, Anderson acceleration, competing low-field branches,
and material phase boundaries are outside this release's regression coverage.

GNU Fortran's optimized build reports possible-uninitialized diagnostics for
the conditionally allocated Anderson-history arrays in the preserved solver
source. No bounds or floating-point exception was reported in the debug
runs of the covered cases. The checks use `anderson_depth=0`.

Numerical source SHA-256 for this release:

`9b9899a4390869479fd448b7a983e84cd377d3455b49b00b15b347bc492b9d6d`

The source is byte-identical to the selected existing solver. Build files,
documentation, the small example, analytic checks, and licensing were prepared
for this public release.
