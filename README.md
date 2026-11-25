# Plaquette CMFT

A Fortran 2008 solver for four-site `J = 5/2` plaquettes with self-consistent
dipolar and quadrupolar fields in stacked square lattices.

The local six-state multiplet is retained at each site, giving a 1296-state
plaquette Hilbert space. The program combines an interacting four-site cluster
with mean fields for the remaining bonds. It searches over user-selected
stacking periods, starting configurations, and magnetic fields.

This is an initial public release of a research solver. The included examples
check well-defined limiting cases; they are not a fitted phase diagram for a
specific material. The numerical source is preserved from the maintainer's
existing research implementation. This release adds a small runnable example,
analytic regression checks, and documentation.

## Build

Requirements:

- A Fortran 2008 compiler; GNU Fortran is the tested compiler.
- LAPACK and BLAS libraries.
- GNU Make.
- Python 3 for the regression checks; only its standard library is used.

On Debian or Ubuntu, the required packages are `gfortran`, `make`,
`liblapack-dev`, `libblas-dev`, and `python3`.

```sh
make
make test
make example
```

Compiler and library settings can be overridden on the command line:

```sh
make FC=gfortran FFLAGS="-O2 -std=f2008 -Wall -Wextra -fopenmp" LDLIBS="-llapack -lblas"
```

OpenMP parallelizes calculations for different planes. For small examples, use
one thread to avoid overhead. When using a threaded BLAS, limit its threads
separately to avoid nested parallelism. The example and test commands set both
`OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.

## Run the supplied example

`make example` runs the model in `examples/saturation.nml` and places its four
CSV output files in `runs/saturation/`. The example is a zero-temperature,
high-field polarized state with a known energy and magnetization. It uses two
field values and one stacking period.

For your own calculation, provide a Fortran `&model` namelist:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./plaquette_cmf path/to/input.nml
```

Output paths are relative to the working directory, and an output prefix's
parent directory must already exist. Files with the same output prefix are
overwritten. Running without an argument looks for `cmf.nml` in the working
directory; this release supplies its example explicitly under `examples/`.

## Model and units

The Hamiltonian contains nearest-neighbor in-plane XXZ exchange, isotropic
first- and second-neighbor interplane exchange, tetragonal Stevens crystal
fields, Zeeman coupling, and the two quadrupoles
`Qx = (Jx Jz + Jz Jx)/2` and `Qy = (Jy Jz + Jz Jy)/2`.

Energies, exchange constants, and crystal-field coefficients are in **meV**;
temperature is in **K**; magnetic field is in **T**. The code uses
`g = 6/7`, `mu_B = 0.05788381806 meV/T`, and
`k_B = 0.08617333262 meV/K`.

See [model conventions](docs/model.md) for the Hamiltonian, bond counting,
field geometry, order parameters, and free-energy normalization.

## Principal inputs

| Input | Meaning |
| --- | --- |
| `j0, j1, j2` | In-plane, nearest-plane, and next-nearest-plane dipolar exchange, in meV. |
| `k0, k1, k2` | Corresponding two-component quadrupolar exchanges, in meV. |
| `delta` | In-plane XXZ anisotropy; the command-line program requires `0 < delta < 1`. |
| `b20, b40, b44` | Crystal-field coefficients in the Stevens convention documented below. |
| `temperature` | Temperature in K; zero invokes the Lanczos ground-state calculation. |
| `theta_deg, inplane` | Tilt from the positive z axis; `inplane='x'` or `'diag'`. |
| `b_start, b_stop, b_step` | Increasing field grid, in T. |
| `n_periods, periods` | Up to eight candidate stacking periods, each from 1 to 80 planes. |
| `n_starts, start_ids` | Number and identifiers of dipolar starting configurations. |
| `n_quad_starts, quad_start_ids` | Number and identifiers of quadrupolar starting configurations. |
| `mixing, tolerance, max_iter` | Self-consistency controls. |
| `random_starts, random_seed, random_scale` | Reproducible perturbations of starting configurations. |
| `seed_stride, keep_branches` | Frequency of fresh searches and number of branches continued. |
| `anderson_depth, dual_mixing` | Anderson acceleration and additional damped-iteration starts. |
| `reverse_sweep` | Also sweep the field downward from the upper endpoint. |
| `output_prefix` | Prefix for the four CSV output files. |

The command-line interface requires all three dipolar exchanges to be positive.
These are interface restrictions in this release; the module-level isolated-ion
test uses the numerical API directly. The source includes a finite-temperature
dense-diagonalization path, which is considerably more expensive than the
zero-temperature calculation.

## Outputs

For an output prefix `run`:

| File | Contents |
| --- | --- |
| `run_summary.csv` | Lowest free energy among converged candidates at each field; order parameters, moment, residual, and search counts. |
| `run_branches.csv` | All converged attempted candidates, including repeated or equivalent states; it is not a count of distinct physical phases. |
| `run_moments.csv` | Dimensionless expectation values of the three angular-momentum components for each plane and sublattice. |
| `run_quadrupoles.csv` | Dimensionless quadrupole expectation values for each plane and sublattice. |

The summary energy is per spin, in meV. Summary magnetic moments are per spin,
in Bohr magnetons. A vanishing order-parameter amplitude makes the corresponding
reported Fourier wavevector physically uninformative.

## Validation and limitations

Run `make test` for analytic checks of saturated-state energy and moment,
stacking-period normalization, and tilted isolated-ion spin and quadrupole
expectation values. The test script exits unsuccessfully if a check fails.
Details, tolerances, and the tested environment are in
[validation notes](docs/validation.md).

Self-consistency convergence establishes a stationary solution within the
chosen ansatz. Comparing a finite set of periods and starting configurations
does not establish a global thermodynamic minimum. Increase the search scope
and check residuals, competing solutions, and field-sweep dependence for
research calculations. A failed search or diagonalization is a numerical
failure and is not evidence that a phase is absent.

The mean-field environment permits two sublattices per plane and selected
periodic stacking patterns. It does not span arbitrary in-plane order.
The validation supplied here is limited to the explicitly documented tests;
it is not a validation of every parameter regime or of quantitative predictions
for a material.

## License and citation

Copyright 2026 Eduardo Miranda. Released under the [MIT License](LICENSE).
LAPACK, BLAS, and the compiler are external dependencies and retain their own
licenses.

For scientific use, cite this repository and the exact commit or version used.
Machine-readable citation metadata is provided in [CITATION.cff](CITATION.cff).
Issues and pull requests are welcome; include the input namelist, compiler
version, and a small reproducing case when reporting a numerical problem.
