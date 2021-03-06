#
#============
# Important:
#============
#	be sure you generate a file called "makefile_SHT_prefix"
#	with the following single line in it:
#
#		SHT_PREFIX := /path/to/prefix/location/for/SHTns

# change this to be any non-empty value to enable debugging flags
debug :=

# change this to be any non-empty value to build the SHTns version
shtns :=

# either INTEL or GNU
compiler := INTEL

# run with openmp support --- seems to be necessary???
omp := t

# use custom library/include commands, defined in "makefile_custom",
# very useful if building on supercomputers with odd requirements
#     to override the automatically generated include_flags/library_flags
#     variables, generate a file called "makefile_custom" with the
#     following definitions:
#         ...
#         force_library := -L<lib_dir> -llibrary_1 ...
#         ...
#         force_include := -I<inc_dir>
#         ...
# change this to be any non-empty value to use custom library/includes
custom :=

# access to the SHT_PREFIX variable
include makefile_SHT_prefix

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
# include/library flags
include_flags :=
library_flags :=

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

  f90_flags += -DINTEL_COMPILER -fpp -shared-intel -module $(mdir) -I$(mdir)
  ifdef omp
    f90_flags += -qopenmp
  endif

  # MKL is only used when SHTns is not
  ifndef shtns
    include_flags += -I$(MKLROOT)/include
    library_flags += -mkl
  endif
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

    f90_flags += -DGNU_COMPILER -cpp -ffixed-line-length-132 -J$(mdir) -I$(mdir)
  ifdef omp
    f90_flags += -fopenmp
  endif

    # MKL is only used when SHTns is not
    ifndef shtns
      include_flags += -I$(MKLROOT)/include
      libs := -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
      library_flags += -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed $(libs)
    endif
  else
    $(error "Compiler must be either INTEL or GNU")
  endif
endif

ifdef shtns
  f90_flags += -DUSE_SHTns
  include_flags += -I$(SHT_PREFIX)/include
  ifdef omp
    library_flags += -L$(SHT_PREFIX)/lib -lfftw3_omp -lfftw3 -lm -lshtns_omp
  else
    library_flags += -L$(SHT_PREFIX)/lib -lfftw3 -lm -lshtns
  endif
  SHT_suf := .SHTns
endif

# use custom defined library/include lines instead of automatically generated ones
ifdef custom
  -include makefile_custom

  ifdef force_library
    library_flags := $(force_library)
  endif

  ifdef force_include
    include_flags := $(force_include)
  endif
endif

ifdef debug
  debug_suf := .debug
endif

suf := $(shell uname)$(comp_suf)$(SHT_suf)$(debug_suf)

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
	$(RM) -rf $(tdir)
	$(RM) -f $(exe)
distclean: realclean
	$(RM) -f ./*.exe
purge: distclean
rebuild: clean
	@$(MAKE) --no-print-directory all

#----------------------------------------------------------------------
# debugging tool: "make print-objects" will print the value of $(objects)
print-%: ; @echo $* is $($*)

