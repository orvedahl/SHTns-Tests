
# change this to be any non-empty value to enable debugging flags
debug :=
compiler := GNU #[either INTEL or GNU]

#######################################################################
# where to find source files
source_dir := src

#----------------------------------------------------------------------
# where to put built modules/objects
tdir := _build
mdir := $(tdir)/m
odir := $(tdir)/o

#----------------------------------------------------------------------
# some basic macros
RM := -rm
bold   = `tput bold`
normal = `tput sgr0`
red    = `tput setaf 1`
green  = `tput setaf 2`

#----------------------------------------------------------------------
# DGEMM include/library flags
DGEMM_INC := -I$(MKLROOT)/include
DGEMM_LIB := -L$(MKLROOT)/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core

include_flags := $(DGEMM_INC)
library_flags := $(DGEMM_LIB)

#----------------------------------------------------------------------
# compiler FLAGS
ifeq ($(findstring INTEL, $(compiler)), INTEL)
  f90_comp := ifort
  comp_suf := .intel

  ifdef debug
    f90_flags = -FR -r8 -O0 -traceback -g -check all -noinline -debug all -warn all
  else
    f90_flags := -FR -r8 -O3
  endif

  f90_flags += -DINTEL_COMPILER -fpp -shared-intel -qopenmp -module $(mdir) -I$(mdir)
else
  ifeq ($(findstring GNU, $(compiler)), GNU)
    f90_comp := gfortran
    comp_suf := .gfortran

    ifdef debug
      dbg1 = -g -fno-range-check -O0 -fbounds-check -fbacktrace -Wuninitialized -Wunused
      dbg2 = -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=nan
      f90_flags := $(dbg1) $(dbg2)
    else
      f90_flags := -O3
    endif

    f90_flags += -DGNU_COMPILER -cpp -ffixed-line-length-132 -fopenmp -J$(mdir) -I$(mdir)
  else
    $(error "Compiler must be either INTEL or GNU")
  endif
endif

ifdef debug
  debug_suf := .debug
endif

suf := $(shell uname)$(comp_suf)$(debug_suf)

exe := xmain.$(suf).exe

#----------------------------------------------------------------------
# what to build and where to look for source files
all: $(exe)

vpath %.f90 $(source_dir)
vpath %.F90 $(source_dir)

#----------------------------------------------------------------------
# dependencies
include makefile_dependencies

#----------------------------------------------------------------------
# link the executable
$(exe): $(objects)
	@echo "${bold}Linking $@ ...${normal}"
	$(f90_comp) $(f90_flags) -o $(exe) $(objects) $(library_flags) $(include_flags)
	@echo
	@echo "${bold}${green}SUCCESS${normal}"
	@echo

#----------------------------------------------------------------------
# compile object files
$(odir)/%.o: %.f90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	@echo "${bold}Building $< ...${normal}"
	$(f90_comp) $(f90_flags) -c $< -o $@ $(include_flags)

$(odir)/%.o: %.F90
	@if [ ! -d $(odir) ]; then mkdir -p $(odir); fi
	@if [ ! -d $(mdir) ]; then mkdir -p $(mdir); fi
	@echo "${bold}Building $< ...${normal}"
	$(f90_comp) $(f90_flags) -c $< -o $@ $(include_flags)

#----------------------------------------------------------------------
# clean up
.PHONY: clean distclean realclean purge
clean:
	$(RM) -f ./*.o $(odir)/*.o
	$(RM) -f ./*.mod $(mdir)/*.mod
realclean: clean
	$(RM) -f $(exe)
distclean: realclean
	$(RM) -rf $(tdir)
	$(RM) -f ./*.exe
purge: distclean

#----------------------------------------------------------------------
# debugging tool: "make print-objects" will print the value of $(objects)
print-%: ; @echo $* is $($*)

