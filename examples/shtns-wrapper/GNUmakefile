#
# General Makefile
#

# Define main makefile directory
makefile_directory := $(HOME)/Programs/Build-System

# Debugging flag
debug :=

# OpenMP/MPI flag
OMP := t
MPI :=

# Print verbosely
verbose := t

# Executable basename
exe_base := sht

# Fortran compiler (choices are "gfortran" or "intel" so far)
f90_comp := gfortran

# Use the Intel Math Kernel Library (this adjusts the intel compiler flags)
Intel_MKL := 

# If using a locally built MPI, specify the location such that 
# the compiler command is $(MPIHOME)/bin/mpif90
MPIHOME := $(HOME)/Local-Installs/openmpi

# List of modules to skip when running through the dependency checker.
# If file my_code.f90 contains:
#       use omp
#       use mpi
#       use my_freq
# and the skip-modules variable is set to:
#       skip_modules := omp mpi
# then the checker will not search for a file that contains a module 
# named "omp" or a file containing module named "mpi" but it will search 
# for a file containing module "my_freq". This allows only "part" of the 
# file to be processed by the dependency checker
skip_modules := iso_c_binding iso_fortran_env

# Path to main source code directory
src_dirs := src ../../utilities

SHT_PREFIX := $(HOME)/cueball/local-software/SHTns
SHT_LIBS := -lfftw3_omp -lfftw3 -lm 

# Other Include directories (the -I flag is specified here)
xtr_inc_dirs := -I$(SHT_PREFIX)/include

# Libraries to load (the -L & -l flags are specified here)
#xtr_libs := -L<lib dir> -l<lib>
xtr_libs := -L$(SHT_PREFIX)/lib $(SHT_LIBS) -lshtns_omp

# Include main build info
include $(makefile_directory)/GMain.mak

