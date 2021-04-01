
src_dirs := src
skip_modules := iso_c_binding iso_fortran_env

build_probin :=

build_info :=

exe_base := sht

debug :=

f90_compiler := intel

f90_flags :=
OMP := t

xtr_f90_flags :=
SHT_PREFIX := $(HOME)/cueball/local-software/SHTns/icc-built
include_flags := -I$(SHT_PREFIX)/include
library_flags := -L$(SHT_PREFIX)/lib -lfftw3_omp -lfftw3 -lm -lshtns_omp

MAKEFILE_HOME_DIR := $(HOME)/Programs/FortranDep/Makefiles
include $(MAKEFILE_HOME_DIR)/GMain.mak

