
# user options
# ------------------
#   sht          -- build using the SHTns library
#   debug        -- compile with debugging flags
#   f90_compiler -- what compiler to use [intel or gfortran]
shtns :=
debug :=
f90_compiler := gfortran

#################################################################

# what to build, what to ignore
src_dirs := src ../utilities
skip_modules := iso_c_binding iso_fortran_env
build_probin := t
build_info := t

exe_base := main

OMP := t

# compiler-dependent preprocessor stuff
ifeq ($(f90_compiler), intel)
  xtr_f90_flags := -DINTEL_COMPILER -fpp
  pp_macros := INTEL_COMPILER=1
else
  xtr_f90_flags := -DGNU_COMPILER -cpp
  pp_macros := GNU_COMPILER=1
endif

#SHTns stuff
ifdef shtns
  xtr_f90_flags += -DUSE_SHTns
  pp_macros += USE_SHTns=1
endif
SHT_PREFIX := $(HOME)/cueball/local-software/SHTns
SHT_INC := -I$(SHT_PREFIX)/include
SHT_LIB := -L$(SHT_PREFIX)/lib -lfftw3_omp -lfftw3 -lm -lshtns_omp

DGEMM_INC := -I$(MKLROOT)/include
DGEMM_LIB := -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core

include_flags := $(SHT_INC) $(DGEMM_INC)
library_flags := $(SHT_LIB) $(DGEMM_LIB)

# main makefile stuff
MAKEFILE_HOME_DIR := $(HOME)/Programs/FortranDep/Makefiles
include $(MAKEFILE_HOME_DIR)/GMain.mak

