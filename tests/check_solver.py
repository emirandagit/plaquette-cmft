#!/usr/bin/env python3
"""Small analytic regression checks for the plaquette CMFT solver.

Uses Python's standard library plus the Fortran compiler and LAPACK/BLAS
already needed to build the solver. No material-specific data are required.
"""

import argparse
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
J = 2.5
G = 6.0 / 7.0
MU_B = 0.05788381806  # meV/T
TOL = 2.0e-8


def close(actual, expected, label, tolerance=TOL):
    if not math.isfinite(actual) or abs(actual - expected) > tolerance:
        raise AssertionError(f"{label}: got {actual:.16g}, expected {expected:.16g}")


def namelist(period, max_iter=10):
    return f"""&model
  j0=0.1, j1=0.02, j2=0.01, delta=0.5
  k0=0.0, k1=0.0, k2=0.0
  b20=0.0, b40=0.0, b44=0.0
  temperature=0.0, theta_deg=0.0, inplane='x'
  b_start=100.0, b_stop=110.0, b_step=10.0
  n_periods=1, periods={period}
  n_starts=1, start_ids=10
  n_quad_starts=1, quad_start_ids=5
  random_starts=0, keep_branches=1, seed_stride=1
  anderson_depth=0, dual_mixing=.false., reverse_sweep=.false.
  mixing=1.0, tolerance=1.0e-9, max_iter={max_iter}
  output_prefix='check'
/
"""


def run(command, directory, env, expected_success=True):
    result = subprocess.run(command, cwd=directory, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            timeout=120, check=False)
    if expected_success and result.returncode != 0:
        raise RuntimeError(f"Command failed: {shlex.join(map(str, command))}\n"
                           f"{result.stdout}\n{result.stderr}")
    return result


def saturated_checks(binary, work, env):
    """An exact polarized product state checks energy and per-spin normalization."""
    for period in (1, 3):
        case = work / f"saturation_nc{period}"
        case.mkdir()
        (case / "input.nml").write_text(namelist(period), encoding="utf-8")
        run([str(binary), "input.nml"], case, env)
        with (case / "check_summary.csv").open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        if len(rows) != 2:
            raise AssertionError(f"Nc={period}: expected two field rows, got {len(rows)}")
        for row, field in zip(rows, (100.0, 110.0)):
            close(float(row["B_T"]), field, "field")
            if int(row["Nc"]) != period:
                raise AssertionError("The output has the wrong stacking period")
            # Two in-plane bonds and one bond of each interplane range per spin.
            expected_f = J * J * (2 * 0.1 * 0.5 + 0.02 + 0.01) - G * MU_B * field * J
            close(float(row["F_meV_per_spin"]), expected_f, f"Nc={period} energy")
            for key in ("Mx_muB", "My_muB", "dip_amplitude", "quad_amplitude"):
                close(float(row[key]), 0.0, f"Nc={period} {key}")
            close(float(row["Mz_muB"]), G * J, f"Nc={period} magnetization")
            residual = float(row["residual"])
            if not math.isfinite(residual) or residual >= 1.0e-9:
                raise AssertionError(f"Nc={period}: unconverged residual {residual}")
            if int(row["converged_branches"]) < 1:
                raise AssertionError("No converged branch was recorded")
        with (case / "check_moments.csv").open(newline="", encoding="utf-8") as stream:
            moments = list(csv.DictReader(stream))
        if len(moments) != 4 * period:
            raise AssertionError("Wrong number of plane/sublattice moment rows")
        for row in moments:
            close(float(row["Jx"]), 0.0, "Jx")
            close(float(row["Jy"]), 0.0, "Jy")
            close(float(row["Jz"]), J, "Jz")
        with (case / "check_quadrupoles.csv").open(newline="", encoding="utf-8") as stream:
            for row in csv.DictReader(stream):
                close(float(row["Qx"]), 0.0, "Qx")
                close(float(row["Qy"]), 0.0, "Qy")
        slope = (float(rows[1]["F_meV_per_spin"]) - float(rows[0]["F_meV_per_spin"])) / 10.0
        close(-slope / MU_B, G * J, "magnetization from energy slope")
        print(f"PASS: saturated CLI state, Nc={period}, two fields")


DRIVER = """program check_isolated_spin
  use plaquette_cmf
  implicit none
  type(parameters) :: par
  type(solution) :: result
  real(dp) :: bvec(3)
  par%j0=0.0_dp; par%j1=0.0_dp; par%j2=0.0_dp
  par%k0=0.0_dp; par%k1=0.0_dp; par%k2=0.0_dp
  par%b20=0.0_dp; par%b40=0.0_dp; par%b44=0.0_dp
  par%temperature=0.0_dp; par%theta_deg=45.0_dp
  par%inplane='diag'; par%mixing=1.0_dp
  par%tolerance=1.0e-9_dp; par%max_iter=4
  par%anderson_depth=0
  call initialize(par)
  call field_vector(par,10.0_dp,bvec)
  call solve_branch(par,1,10,5,0,bvec,result=result)
  if (.not.result%converged) error stop 'Isolated spin did not converge'
  write(*,'(*(es25.16,1x))') result%f, result%m(:,1,1), &
    result%qm(:,1,1), result%m(:,2,1), result%qm(:,2,1), result%residual
end program check_isolated_spin
"""


def tilted_spin_check(source, args, work, env):
    """An independent spin coherent state exercises complex Jy and quadrupoles."""
    text = source.read_text(encoding="utf-8")
    marker = "program run_plaquette_cmf"
    if text.count(marker) != 1:
        raise RuntimeError("Cannot identify the module/program boundary in the solver source")
    case = work / "isolated_spin"
    case.mkdir()
    # The original module is compiled unchanged; only the CLI program is omitted
    # so the driver can call its public API with zero exchange constants.
    (case / "module.f90").write_text(text.split(marker)[0], encoding="utf-8")
    (case / "driver.f90").write_text(DRIVER, encoding="utf-8")
    command = (shlex.split(args.fc) + shlex.split(args.fflags)
               + ["-o", "isolated_spin", "module.f90", "driver.f90"]
               + shlex.split(args.ldflags) + shlex.split(args.ldlibs))
    run(command, case, env)
    result = run([str(case / "isolated_spin")], case, env)
    values = [float(word) for word in result.stdout.split()]
    if len(values) != 12:
        raise AssertionError(f"Expected 12 driver values, got {len(values)}")
    axis = (0.5, 0.5, 1.0 / math.sqrt(2.0))
    moment = [J * n for n in axis]
    quad = [J * (J - 0.5) * axis[i] * axis[2] for i in (0, 1)]
    expected = [-G * MU_B * 10.0 * J] + moment + quad + moment + quad
    for i, (actual, target) in enumerate(zip(values[:-1], expected)):
        close(actual, target, f"isolated-spin component {i}")
    if not math.isfinite(values[-1]) or values[-1] >= 1.0e-9:
        raise AssertionError("Isolated-spin fixed-point residual is too large")
    print("PASS: tilted isolated-spin energy, both sublattice moments and quadrupoles")


def nonconvergence_check(binary, work, env):
    case = work / "nonconvergence"
    case.mkdir()
    (case / "input.nml").write_text(namelist(1, max_iter=1), encoding="utf-8")
    result = run([str(binary), "input.nml"], case, env, expected_success=False)
    if result.returncode == 0:
        raise AssertionError("A deliberately unconverged run returned success")
    if "No branch converged at a field" not in result.stdout + result.stderr:
        raise AssertionError("Failed run did not report the expected convergence failure")
    print("PASS: an insufficient iteration budget reports failure")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default=str(ROOT / "plaquette_cmf"))
    parser.add_argument("--source", default=str(ROOT / "plaquette_cmf.f90"))
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--fflags", default=os.environ.get("FFLAGS", "-O2 -std=f2008 -fopenmp"))
    parser.add_argument("--ldflags", default=os.environ.get("LDFLAGS", ""))
    parser.add_argument("--ldlibs", default=os.environ.get("LDLIBS", "-llapack -lblas"))
    args = parser.parse_args()
    binary = Path(args.executable).resolve(strict=True)
    source = Path(args.source).resolve(strict=True)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    start = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="plaquette-cmft-check-") as temporary:
        work = Path(temporary)
        saturated_checks(binary, work, env)
        tilted_spin_check(source, args, work, env)
        nonconvergence_check(binary, work, env)
    print(f"All four analytic/behavior checks passed in {time.perf_counter() - start:.2f} s")


if __name__ == "__main__":
    main()
