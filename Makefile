FC = gfortran
FFLAGS = -O3 -std=f2008 -Wall -Wextra -fopenmp
LDFLAGS =
LDLIBS = -llapack -lblas
PYTHON = python3

all: plaquette_cmf

plaquette_cmf: plaquette_cmf.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

test: plaquette_cmf
	$(PYTHON) tests/check_solver.py --executable ./plaquette_cmf --fc "$(FC)" --fflags "$(FFLAGS)" --ldflags "$(LDFLAGS)" --ldlibs "$(LDLIBS)"

example: plaquette_cmf
	mkdir -p runs/saturation
	cd runs/saturation && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../plaquette_cmf ../../examples/saturation.nml

clean:
	rm -f plaquette_cmf *.o *.mod

.PHONY: all test example clean
