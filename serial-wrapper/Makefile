
shtns :=
debug :=
SHT_PREFIX := $(HOME)/cueball/local-software/SHTns
################################################################

SHT_INC := -I$(SHT_PREFIX)/include
SHT_LIB := -L$(SHT_PREFIX)/lib -lfftw3_omp -lfftw3 -lm -lshtns_omp

DGEMM_INC := -I$(MKLROOT)/include
DGEMM_LIB := -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core

include_flags := $(SHT_INC) $(DGEMM_INC)
library_flags := $(SHT_LIB) $(DGEMM_LIB)

f90_comp := ifort

ifdef debug
f90_flags := -FR -r8 -O0 -traceback -g -check all -noinline -debug all -warn all -shared-intel -qopenmp
else
f90_flags := -FR -r8 -O3 -shared-intel -qopenmp
endif

f90_flags += -DINTEL_COMPILER -fpp
ifdef shtns
f90_flags += -DUSE_SHTns
endif

f90_compile := $(f90_flags)

all: xmain

Legendre_Polynomials.o: src/Legendre_Polynomials.F90
cheb_grids.o: src/cheb_grids.f90
Structures.o: src/Structures.F90
probin.o: src/probin.f90
ifdef shtns
Legendre_Transforms_SHTns.o: Legendre_Polynomials.o Structures.o src/Legendre_Transforms_SHTns.F90
else
Legendre_Transforms.o: Legendre_Polynomials.o Structures.o src/Legendre_Transforms.F90
endif
ifdef shtns
test_suite.o: Structures.o cheb_grids.o Legendre_Polynomials.o Legendre_Transforms_SHTns.o src/test_suite.f90
else
test_suite.o: Structures.o cheb_grids.o Legendre_Polynomials.o Legendre_Transforms.o src/test_suite.f90
endif
ifdef shtns
main.o: probin.o test_suite.o Legendre_Polynomials.o cheb_grids.o Legendre_Transform_SHTns.o src/main.f90
else
main.o: probin.o test_suite.o Legendre_Polynomials.o cheb_grids.o src/main.f90
endif

objects := probin.o cheb_grids.o Structures.o Legendre_Polynomials.o
ifdef shtns
objects += Legendre_Transforms_SHTns.o
else
objects += Legendre_Transforms.o
endif
objects += test_suite.o main.o

xmain: $(objects)
	$(f90_comp) $(f90_compile) $(library_flags) $(include_flags) $(objects) -o xmain
%.o: src/%.f90
	$(f90_comp) $(f90_compile) -c $< -o $@ $(include_flags)
%.o: src/%.F90
	$(f90_comp) $(f90_compile) -c $< -o $@ $(include_flags)

.PHONY: clean distclean realclean purge
clean:
	rm -f ./*.o
	rm -f ./*.mod
distclean: clean
	rm -f xmain
realclean: distclean
purge : distclean

