
src_dirs := src ../../utilities
skip_modules := iso_c_binding iso_fortran_env

# this is actually used
build_probin := t

# this is not strictly used, but codes in ../../utilities do use it
# and I'm too lazy to remove them from the GPackage.mak file, because
# they are used in other test/example cases.
build_info := t

exe_base := main

verbose := t
debug :=

f90_comp := gfortran

f90_flags :=
OMP := t

xtr_f90_flags :=
SHT_PREFIX := $(HOME)/cueball/local-software/SHTns
include_flags := -I$(SHT_PREFIX)/include
library_flags := -L$(SHT_PREFIX)/lib -lfftw3_omp -lfftw3 -lm -lshtns_omp

MAKEFILE_HOME_DIR := $(HOME)/Programs/FortranDep/Makefiles
include $(MAKEFILE_HOME_DIR)/GMain.mak

